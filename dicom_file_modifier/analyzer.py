#!/usr/bin/env python3
"""
Analyzer Analyse von Zielgebieten und Risikoorganen
aus DICOM RT Structure Set Dateien.

Berechnet pro Struktur:
  - Volumen (cm³)
  - Schwerpunkt (x, y, z in mm)
  - Bounding Box
  - Formmetriken: Sphärizität, Kompaktheit, Elongation

Berechnet zwischen Strukturen:
  - Minimaler Abstand (mm)
  - Hausdorff-Abstand (mm)
  - Schwerpunkt-Abstand (mm)

Verwendung:
  python -m dicom_file_modifier.analyzer rtstruct.dcm [--targets PTV,CTV,GTV] [--oars Parotis,Rueckenmark]
  python -m dicom_file_modifier.analyzer rtstruct.dcm --list    # Nur Strukturnamen auflisten
  python -m dicom_file_modifier.analyzer --self-test            # Synthetische Konsistenztests (Löcher, z-Lücken)

Benötigte Packages:
  pip install pydicom numpy scipy shapely matplotlib
"""

from __future__ import annotations

import argparse
import contextlib
import contextvars
import json
import re
import sys
from functools import reduce
from pathlib import Path
from typing import Optional

import numpy as np
import pydicom
from scipy.spatial import ConvexHull
from scipy.spatial.distance import directed_hausdorff, pdist
from shapely.geometry import Polygon
from shapely.ops import unary_union
from shapely.validation import make_valid

from .dicom_utils import find_point_markers, get_rs_frame_of_references
from .issues import Issue

# Fester Zufallsgenerator: Distanz-/Subsample-Operationen sollen reproduzierbar
# sein (vorher unverseedetes np.random.choice -> nicht-deterministische QA-Werte).
# ``analyze_rtstruct`` setzt je Lauf einen frischen Generator (Seed 0): ein
# langlebiger Prozess (Worker/GUI) liefert dann dieselben Zahlen wie ein
# frischer CLI-Aufruf.  Ausserhalb eines Laufs gilt ein Generator je Kontext.
_RNG_VAR: contextvars.ContextVar = contextvars.ContextVar("analyzer_rng", default=None)
# Waehrend ``analyze_rtstruct``: Sammelliste fuer Befunde und die aktuelle ROI
_ISSUES_VAR: contextvars.ContextVar = contextvars.ContextVar("analyzer_issues", default=None)
_ROI_VAR: contextvars.ContextVar = contextvars.ContextVar("analyzer_roi", default="")


def _rng() -> np.random.Generator:
    rng = _RNG_VAR.get()
    if rng is None:
        rng = np.random.default_rng(0)
        _RNG_VAR.set(rng)
    return rng


def _note(level: str, code: str, message: str, exc: Optional[BaseException] = None) -> None:
    """Befund des laufenden ``analyze_rtstruct`` sammeln (ausserhalb: verworfen)."""
    sink = _ISSUES_VAR.get()
    if sink is None:
        return
    roi = _ROI_VAR.get()
    sink.append(Issue(level, code, f"{roi}: {message}" if roi else message,
                      detail=f"{type(exc).__name__}: {exc}" if exc is not None else ""))


@contextlib.contextmanager
def _roi_context(name: str):
    token = _ROI_VAR.set(name)
    try:
        yield
    finally:
        _ROI_VAR.reset(token)


def parse_name_list(spec: Optional[str]) -> Optional[list]:
    """``"PTV, GTV"`` -> ``["PTV", "GTV"]`` (getrimmt, leere Eintraege entfallen); leer -> None."""
    names = [t.strip() for t in (spec or "").split(",") if t.strip()]
    return names or None


# ---------------------------------------------------------------------------
# 1. DICOM RTSTRUCT einlesen
# ---------------------------------------------------------------------------

def load_rtstruct(filepath: str) -> pydicom.Dataset:
    """Lädt eine DICOM RTSTRUCT Datei und prüft die Modalität."""
    ds = pydicom.dcmread(filepath)
    if ds.Modality != "RTSTRUCT":
        raise ValueError(f"Datei ist keine RTSTRUCT (Modalität: {ds.Modality})")
    return ds


def get_structure_names(ds: pydicom.Dataset) -> dict[int, str]:
    """Gibt ein Dictionary {ROI-Nummer: Name} zurück."""
    return {
        roi.ROINumber: roi.ROIName
        for roi in ds.StructureSetROISequence
    }


def get_structure_type(ds: pydicom.Dataset) -> dict[int, str]:
    """Gibt ein Dictionary {ROI-Nummer: RT ROI Interpreted Type} zurück."""
    type_map = {}
    if hasattr(ds, "RTROIObservationsSequence"):
        for obs in ds.RTROIObservationsSequence:
            roi_num = obs.ReferencedROINumber
            rt_type = getattr(obs, "RTROIInterpretedType", "UNKNOWN")
            type_map[roi_num] = rt_type
    return type_map


# ---------------------------------------------------------------------------
# 1b. Strukturklassifikation (kategorisiert ROIs für Auswertung + Plots)
# ---------------------------------------------------------------------------

# Kategorien (Reihenfolge = Sortier-/Plot-Reihenfolge)
CAT_TARGET = "TARGET"          # echte GTV/PTV/CTV/ITV
CAT_OAR_SERIAL = "OAR_SERIAL"  # serielle (Maximaldosis-kritische) Risikoorgane
CAT_OAR_PARALLEL = "OAR_PARALLEL"  # parallele (Volumeneffekt-)Risikoorgane
CAT_HELPER = "HELPER"          # Hilfs-/Planungs-/Optimierungs-/Vereinigungsstrukturen
CAT_EXTERNAL = "EXTERNAL"      # Außenkontur/Body
CAT_MARKER = "MARKER"          # POINT-Marker (Fiducials, Isozentren)

CATEGORY_ORDER = [CAT_TARGET, CAT_OAR_SERIAL, CAT_OAR_PARALLEL,
                  CAT_HELPER, CAT_EXTERNAL, CAT_MARKER]

# Serielle OARs: Schlüsselwörter (case-insensitive, Teilstring)
_SERIAL_OAR_KEYWORDS = (
    "hirnstamm", "brainstem", "rueckenmark", "rückenmark", "spinal", "myelon",
    "sehnerv", "optic", "chiasma", "chiasm", "hypophyse", "pituitary",
)
# Hilfsstruktur-Erkennung über den Namen (Vereinigungen, Dosis-Shells, Opt-
# Strukturen). Bewusst eng gehalten: ein echtes Organ wie "Hirn gesamt" darf
# NICHT als Hilfsstruktur gelten. Die Konvention dieses Datensatzes präfixt
# Hilfsstrukturen mit "h_" bzw. "opt"; zusätzlich generische Shell-/Margin-Namen.
_HELPER_NAME_RE = re.compile(
    r"(^h_|opt[\s_]*system|^opt_|\bring\b|\bshell\b|\+\s*\d+\s*mm"
    r"|^iso_?\d|_x_iso\d|_minus_iso\d|^iso\d+_minus_)",   # dose_indices-ROIs
    re.IGNORECASE,
)
# RT-Typen, die nie klinische Strukturen sind (TPS-Isodosen, Optimierungs-/
# Kontrollstrukturen): Eclipse exportiert Isodosen-Strukturen als CONTROL.
_HELPER_RT_TYPES = ("CONTROL", "DOSE_REGION")
# Echte Zielvolumina: Name beginnt mit GTV/PTV/CTV/ITV (gefolgt von _ oder Ziffer)
_TARGET_NAME_RE = re.compile(r"^(gtv|ptv|ctv|itv)[ _0-9]", re.IGNORECASE)
# Läsions-Schlüssel zum Paaren von GTV mit seinem PTV (z.B. "GTV_1"/"PTV_1" -> "1")
_LESION_KEY_RE = re.compile(r"^(?:gtv|ptv|ctv|itv)_(.+)$", re.IGNORECASE)


def classify_structure(name: str, rt_type: str, geom_types: set[str]) -> str:
    """
    Ordnet eine ROI genau einer Kategorie zu.

    Reihenfolge der Regeln ist wichtig: Marker/External zuerst, dann Hilfs-
    strukturen *per Name* (da das DICOM `RTROIInterpretedType` Vereinigungs-
    PTVs wie ``h_PTV_gesamt`` fälschlich als ``GTV`` und Opt-Strukturen als
    ``ORGAN`` taggt), dann echte Targets, zuletzt OAR-Aufteilung seriell/parallel.
    """
    rt = (rt_type or "").upper()
    nm = name or ""

    # 1) POINT-Marker / Fiducials
    if "POINT" in geom_types or rt == "MARKER":
        return CAT_MARKER
    # 2) Außenkontur / Body
    if rt == "EXTERNAL" or re.search(r"aussenkontur|außenkontur|external|\bbody\b|koerper|körper",
                                     nm, re.IGNORECASE):
        return CAT_EXTERNAL
    # 3) Hilfs-/Planungsstrukturen (Name-basiert, überschreibt fehlerhaftes RT-Type;
    #    CONTROL/DOSE_REGION = TPS-Isodosen und Kontrollstrukturen)
    if _HELPER_NAME_RE.search(nm) or rt in _HELPER_RT_TYPES:
        return CAT_HELPER
    # 4) Echte Zielvolumina
    if _TARGET_NAME_RE.match(nm) or rt in ("PTV", "CTV", "GTV", "ITV", "TV"):
        return CAT_TARGET
    # 5) Risikoorgane: seriell vs. parallel
    low = nm.lower()
    if any(k in low for k in _SERIAL_OAR_KEYWORDS):
        return CAT_OAR_SERIAL
    return CAT_OAR_PARALLEL


def lesion_key(name: str) -> Optional[str]:
    """Extrahiert den Läsions-Schlüssel eines Targets (Teil nach GTV_/PTV_)."""
    m = _LESION_KEY_RE.match(name or "")
    return m.group(1) if m else None


# ---------------------------------------------------------------------------
# 2. Konturen extrahieren
# ---------------------------------------------------------------------------

def extract_contours(ds: pydicom.Dataset, roi_number: int) -> list[np.ndarray]:
    """
    Extrahiert die Konturen einer Struktur als Liste von Nx3 Arrays.
    Jedes Array enthält die (x, y, z) Koordinaten einer Kontur-Schicht.
    """
    contours = []
    for roi_contour in ds.ROIContourSequence:
        if roi_contour.ReferencedROINumber != roi_number:
            continue
        if not hasattr(roi_contour, "ContourSequence"):
            continue
        for contour in roi_contour.ContourSequence:
            pts = np.array(contour.ContourData).reshape(-1, 3)
            contours.append(pts)
    return contours


def contours_to_points(contours: list[np.ndarray]) -> np.ndarray:
    """Fasst alle Konturpunkte zu einem einzigen Nx3 Array zusammen."""
    if not contours:
        return np.empty((0, 3))
    return np.vstack(contours)


# ---------------------------------------------------------------------------
# 3. Schichten, XOR-Geometrie und Volumen
# ---------------------------------------------------------------------------

# Konturen, deren z-Werte sich um höchstens diese Toleranz unterscheiden, liegen
# auf derselben Schicht (DICOM-z ist pro Ebene exakt; die Toleranz fängt nur
# Rundungsrauschen ab).
_Z_MERGE_TOL_MM = 0.05
# Ein z-Abstand > _GAP_FACTOR x nominale Schichtdicke gilt als Lücke
# (fehlende Schicht oder räumlich getrennte Komponente) und wird nie überbrückt.
_GAP_FACTOR = 1.5


def _polygonal(geom):
    """Nur die flächigen Teile einer Geometrie (``make_valid`` kann Linien liefern)."""
    if geom.geom_type in ("Polygon", "MultiPolygon"):
        return geom
    parts = [g for g in getattr(geom, "geoms", [])
             if g.geom_type in ("Polygon", "MultiPolygon")]
    return unary_union(parts) if parts else Polygon()


def _contour_polygon(pts: np.ndarray):
    """Gültiges Shapely-(Multi)Polygon einer Kontur (x, y = erste zwei Spalten)
    oder ``None`` bei < 3 Punkten bzw. verschwindender Fläche.

    Ungültige Ringe (Selbstschnitt-Spitzen, wie sie Eclipse gelegentlich
    schreibt) werden mit ``make_valid`` repariert; das erhält -- anders als
    ``buffer(0)`` -- beide Lappen einer Schleife.  Die Reparatur MUSS vor jeder
    Overlay-Operation erfolgen, sonst wirft GEOS eine TopologyException.
    """
    if len(pts) < 3:
        return None
    try:
        poly = Polygon(np.asarray(pts, dtype=float)[:, :2])
        if not poly.is_valid:
            try:
                poly = _polygonal(make_valid(poly))
            except Exception as e:
                _note("info", "ANA.CONTOUR_REPAIR_FALLBACK",
                      "ungueltige Kontur mit buffer(0) statt make_valid repariert", e)
                poly = poly.buffer(0)
        if poly.is_empty or poly.area <= 0:
            return None
        return poly
    except Exception as e:
        _note("warning", "ANA.CONTOUR_DROPPED", "Kontur nicht auswertbar und uebergangen", e)
        return None


def polygon_area(pts_2d: np.ndarray) -> float:
    """Fläche eines 2D-Polygons (Shoelace via Shapely; 0.0 bei ungültig)."""
    poly = _contour_polygon(np.asarray(pts_2d, dtype=float))
    return float(poly.area) if poly is not None else 0.0


def _group_slices(contours: list[np.ndarray],
                  tol: float = _Z_MERGE_TOL_MM) -> list[tuple[float, list[np.ndarray]]]:
    """Gruppiert Konturen nach Schicht, z aufsteigend: ``[(z, [pts, ...]), ...]``.

    Konturen mit z-Abstand <= ``tol`` zur laufenden Gruppe gehören zur selben
    Schicht; z der Gruppe ist der Mittelwert.  Ersetzt das frühere Runden auf
    3 bzw. 4 Dezimalen (Volumen vs. Raster), das inkonsistent war.
    """
    order = sorted(range(len(contours)), key=lambda i: float(contours[i][0, 2]))
    groups = []                      # [([z, ...], [pts, ...]), ...]
    for i in order:
        z = float(contours[i][0, 2])
        if groups and z - groups[-1][0][-1] <= tol:
            groups[-1][0].append(z)
            groups[-1][1].append(contours[i])
        else:
            groups.append(([z], [contours[i]]))
    return [(float(np.mean(zs)), pts) for zs, pts in groups]


def _nominal_slice_spacing(zs) -> tuple[Optional[float], int]:
    """Nominale Schichtdicke und Anzahl der z-Lücken aus Schicht-z-Werten.

    Der *Mittelwert* der z-Abstände ist (z_max - z_min)/(n - 1) und wird bei
    Lücken (fehlende Schicht, getrennte Läsionen in einer ROI) aufgebläht.
    Stattdessen: Median der Abstände, die höchstens ``_GAP_FACTOR`` x den
    kleinsten Abstand betragen (robust gegen Lücken und Rundungsrauschen).
    Abstände > ``_GAP_FACTOR`` x Nominalwert werden als Lücken gezählt und
    NICHT überbrückt.  Liefert ``(None, 0)`` bei weniger als zwei Schichten.
    """
    zs = np.sort(np.asarray(list(zs), dtype=float))
    if len(zs) < 2:
        return None, 0
    d = np.diff(zs)
    d_min = float(d.min())
    if d_min <= 0:
        return None, 0
    dz = float(np.median(d[d <= _GAP_FACTOR * d_min]))
    n_gaps = int(np.count_nonzero(d > _GAP_FACTOR * dz))
    return dz, n_gaps


def _slice_geometry(polys):
    """XOR (Even-Odd) aller Polygone einer Schicht.

    Das ist die DICOM-Semantik von ``CLOSED_PLANAR_XOR`` und die Eclipse-
    Konvention für mehrere Konturen einer ROI auf einer Ebene: eine
    verschachtelte Innenkontur ist ein Loch, eine Insel im Loch zählt wieder,
    getrennte Inseln addieren sich; reihenfolgeunabhängig.  Nebenwirkung:
    identische Doppelkonturen löschen sich aus (wie im TPS).
    """
    if len(polys) == 1:
        return polys[0]
    return reduce(lambda a, b: a.symmetric_difference(b), polys)


def _slice_geometries(contours: list[np.ndarray]) -> list[tuple[float, object]]:
    """Pro Schicht die XOR-Geometrie: ``[(z, geom), ...]``, nur Flächen > 0."""
    out = []
    for z, pts_list in _group_slices(contours):
        polys = [p for p in (_contour_polygon(pts) for pts in pts_list) if p is not None]
        if not polys:
            continue
        try:
            geom = _slice_geometry(polys)
        except Exception as e:
            _note("warning", "ANA.SLICE_SKIPPED",
                  f"Schicht z={z:.2f} mm nicht kombinierbar (XOR) und uebergangen", e)
            continue
        if geom.is_empty or geom.area <= 0:      # z.B. Selbst-XOR -> leer
            continue
        out.append((z, geom))
    return out


def compute_volume(contours: list[np.ndarray], geoms=None) -> float:
    """
    Berechnet das Volumen in cm³ als Schichtstapel:
    V = dz_nom * Σ_k A_k,  A_k = XOR-Fläche aller Konturen der Schicht k.

    dz_nom ist die nominale Schichtdicke (``_nominal_slice_spacing``);
    z-Lücken werden nicht überbrückt.  ``geoms`` (aus ``_slice_geometries``)
    kann übergeben werden, um die XOR-Geometrien nicht doppelt zu berechnen.
    """
    if len(contours) < 2:
        return 0.0
    dz, _ = _nominal_slice_spacing(z for z, _ in _group_slices(contours))
    if dz is None:
        return 0.0
    if geoms is None:
        geoms = _slice_geometries(contours)
    total_area = sum(g.area for _, g in geoms)
    return total_area * dz / 1000.0  # mm³ -> cm³


# ---------------------------------------------------------------------------
# 4. Schwerpunkt
# ---------------------------------------------------------------------------

def compute_centroid(contours: list[np.ndarray], geoms=None) -> np.ndarray:
    """
    Flächengewichteter Schwerpunkt (x, y, z) in mm.

    Gewicht jeder Schicht ist ihre XOR-Fläche, die In-Plane-Position der
    Flächenschwerpunkt der XOR-Region (Shapely ``centroid`` = Flächenmoment,
    NICHT das Eckpunkt-Mittel, das zu dicht besetzten Randabschnitten hin
    verzerrt wäre).  Ein Loch verschiebt den Schwerpunkt damit korrekt von
    sich weg, statt ihn -- wie bei positiver Gewichtung der Innenkontur --
    zu sich hin zu ziehen.
    """
    if not contours:
        return np.array([0.0, 0.0, 0.0])
    if geoms is None:
        geoms = _slice_geometries(contours)

    weighted_sum = np.zeros(3)
    total_weight = 0.0
    for z, g in geoms:
        c = g.centroid
        if c.is_empty:
            continue
        weighted_sum += g.area * np.array([c.x, c.y, z])
        total_weight += g.area

    if total_weight == 0:
        return contours_to_points(contours).mean(axis=0)
    return weighted_sum / total_weight


# ---------------------------------------------------------------------------
# 5. Formanalyse
# ---------------------------------------------------------------------------

def rasterize_contours(contours: list[np.ndarray], xc: np.ndarray,
                       yc: np.ndarray, zc: np.ndarray) -> np.ndarray:
    """
    Rastert Konturen per XOR auf ein gegebenes Gitter (Voxelmittelpunkte
    ``xc``, ``yc``, ``zc`` in mm, jeweils aufsteigend) -> bool-Maske (nz, ny, nx).

    Jede Kontur wird der nächsten z-Ebene zugeordnet (Toleranz: halber
    Ebenenabstand, sonst übersprungen) und per XOR mit der Ebene kombiniert,
    d.h. Löcher/Inseln exakt wie in ``_slice_geometry``.  Bewusst KEIN
    Matplotlib-Compound-Path: dessen ``contains_points`` füllt verschachtelte
    Konturen (getestet).  ``contains_points`` läuft nur über das Bounding-Box-
    Teilgitter jeder Kontur (halbiert die Laufzeit, Pflicht für feine Raster).

    Wiederverwendbar für beliebige Gitter, z.B. ein feines Dosisraster mit
    z auf den Dosisebenen.
    """
    from matplotlib.path import Path as MplPath

    xc = np.asarray(xc, dtype=float)
    yc = np.asarray(yc, dtype=float)
    zc = np.asarray(zc, dtype=float)
    z_tol = 0.5 * float(np.min(np.diff(zc))) if len(zc) > 1 else np.inf
    mask = np.zeros((len(zc), len(yc), len(xc)), dtype=bool)

    for pts in contours:
        if len(pts) < 3:
            continue
        z = float(pts[0, 2])
        k = int(np.argmin(np.abs(zc - z)))
        if abs(zc[k] - z) > z_tol + 1e-6:
            continue                        # Kontur liegt auf keiner Gitterebene
        i0 = int(np.searchsorted(xc, pts[:, 0].min()))
        i1 = int(np.searchsorted(xc, pts[:, 0].max(), side="right"))
        j0 = int(np.searchsorted(yc, pts[:, 1].min()))
        j1 = int(np.searchsorted(yc, pts[:, 1].max(), side="right"))
        if i1 <= i0 or j1 <= j0:
            continue
        gx, gy = np.meshgrid(xc[i0:i1], yc[j0:j1])          # (j1-j0, i1-i0)
        try:
            inside = MplPath(pts[:, :2]).contains_points(
                np.column_stack([gx.ravel(), gy.ravel()]))
        except Exception as e:
            _note("warning", "ANA.RASTER_CONTOUR_SKIPPED", f"Kontur bei z={z:.2f} mm beim Rastern uebergangen", e)
            continue
        mask[k, j0:j1, i0:i1] ^= inside.reshape(j1 - j0, i1 - i0)
    return mask


def _rasterize_structure(contours: list[np.ndarray],
                         target_dim: int = 96):
    """
    Rastert die gestapelten Konturen in eine binäre 3D-Voxelmaske.

    Liefert ``(mask, spacing, centers)`` mit ``mask`` (nz, ny, nx) bool,
    ``spacing`` = (sz, sy, sx) in mm und ``centers`` = (xc, yc, zc) der
    Voxelmittelpunkte in LPS-mm, oder ``None`` bei zu wenig Daten.  Die
    z-Ebenen liegen exakt auf den Konturebenen (nominale Schichtdicke, siehe
    ``_nominal_slice_spacing``); z-Lücken ergeben leere Ebenen, so dass
    ``n_components`` getrennte Teile korrekt zählt.  Eine konsistente Maske
    ist die Grundlage für mathematisch *beschränkte* Sphärizität/Solidität.
    """
    try:
        import matplotlib.path  # noqa: F401  (nur Verfügbarkeit prüfen)
    except Exception as e:
        _note("warning", "ANA.RASTER_UNAVAILABLE",
              "matplotlib.path fehlt; Formmetriken nur aus der konvexen Huelle", e)
        return None

    all_pts = contours_to_points(contours)
    if len(all_pts) < 4:
        return None

    x0, y0, z0 = all_pts.min(axis=0)
    x1, y1, z1 = all_pts.max(axis=0)
    span_x, span_y, span_z = x1 - x0, y1 - y0, z1 - z0
    if span_x <= 0 or span_y <= 0:
        return None

    dz, _ = _nominal_slice_spacing(z for z, _ in _group_slices(contours))
    s_z = dz if dz is not None else 1.0

    s_xy = max(span_x, span_y) / target_dim
    s_xy = float(np.clip(s_xy, 0.3, 5.0))
    pad = 2

    nx = int(np.ceil(span_x / s_xy)) + 1 + 2 * pad
    ny = int(np.ceil(span_y / s_xy)) + 1 + 2 * pad
    nz = int(round(span_z / s_z)) + 1 + 2 * pad
    if nx * ny * nz > 6_000_000:          # Sicherheitskappe gegen Speicher-Spikes
        # nz ist durch die Schichtdicke fest -> nur in-plane vergröbern
        # (Quadratwurzel, nicht Kubikwurzel).
        s_xy *= (nx * ny * nz / 6_000_000) ** 0.5
        nx = int(np.ceil(span_x / s_xy)) + 1 + 2 * pad
        ny = int(np.ceil(span_y / s_xy)) + 1 + 2 * pad

    ox, oy = x0 - pad * s_xy, y0 - pad * s_xy
    oz = z0 - (pad + 0.5) * s_z            # -> zc[pad] == z0 exakt
    xc = ox + (np.arange(nx) + 0.5) * s_xy
    yc = oy + (np.arange(ny) + 0.5) * s_xy
    zc = oz + (np.arange(nz) + 0.5) * s_z

    mask = rasterize_contours(contours, xc, yc, zc)
    if not mask.any():
        return None
    return mask, (s_z, s_xy, s_xy), (xc, yc, zc)


def _voxel_shape_metrics(contours: list[np.ndarray],
                         all_pts: np.ndarray) -> Optional[dict]:
    """Beschränkte Sphärizität/Solidität aus einer konsistenten Voxelmaske.

    Sphärizität = π^(1/3)(6V)^(2/3) / A   (V, A beide aus derselben Maske)
    Solidität   = V_Maske / V_konvexe-Hülle  (beide voxelbasiert -> ≤ 1)
    Beide sind hier mathematisch auf (0, 1] beschränkt; das alte Verfahren
    mischte planimetrisches V mit Hüllen-A/V und lieferte unphysikalische
    Werte > 1.
    """
    raster = _rasterize_structure(contours)
    if raster is None:
        return None
    mask, (s_z, s_y, s_x), _axes = raster
    vox_vol = s_x * s_y * s_z
    v_mask = float(mask.sum()) * vox_vol
    if v_mask <= 0:
        return None

    out = {"volume_voxel_cm3": round(v_mask / 1000.0, 3)}

    # Zusammenhangskomponenten (Erkennung von Vereinigungs-/Mehrkomponenten-ROIs)
    try:
        from scipy import ndimage
        out["n_components"] = int(ndimage.label(mask)[1])
    except Exception as e:
        _note("warning", "ANA.COMPONENTS_FAILED", "Zusammenhangskomponenten nicht bestimmbar (als 1 gezaehlt)", e)
        out["n_components"] = 1

    # Oberfläche via Marching Cubes (gleiche Maske wie V) -> Sphärizität ≤ 1
    sph = None
    try:
        from skimage import measure as skmeasure
        verts, faces, _, _ = skmeasure.marching_cubes(
            mask.astype(np.float32), level=0.5, spacing=(s_z, s_y, s_x)
        )
        area = float(skmeasure.mesh_surface_area(verts, faces))
        if area > 0:
            sph = (np.pi ** (1 / 3) * (6 * v_mask) ** (2 / 3)) / area
    except Exception as e:
        _note("warning", "ANA.SURFACE_FAILED", "Oberflaeche (Marching Cubes) nicht bestimmbar; keine Sphaerizitaet", e)
        sph = None
    if sph is not None:
        out["sphericity"] = round(float(min(sph, 1.0)), 4)

    # Solidität = V_Maske / V_konvexe-Hülle.  Beide sind echte Volumina (die
    # gefüllte Voxelmaske, NICHT das aufgeblähte planimetrische Volumen), und
    # das Konturgebiet liegt in der konvexen Hülle -> Verhältnis ≤ 1 (bis auf
    # Sub-Voxel-Diskretisierung, daher geklippt).  Die exakte Hüllen-Volumen-
    # berechnung vermeidet den teuren Halbraumtest über das gesamte Gitter.
    try:
        v_hull = float(ConvexHull(all_pts).volume)
        if v_hull > 0:
            out["solidity"] = round(float(min(v_mask / v_hull, 1.0)), 4)
    except Exception as e:
        _note("warning", "ANA.HULL_FAILED", "konvexe Huelle nicht bestimmbar; keine Soliditaet", e)

    return out


def compute_shape_metrics(contours: list[np.ndarray], volume_cm3: float) -> dict:
    """
    Berechnet Formmetriken (alle dimensionslosen Ratios sind beschränkt ≤ 1):
    - Sphärizität: π^(1/3)(6V)^(2/3)/A aus einer konsistenten Voxelmaske
                   (1.0 = perfekte Kugel)
    - Solidität:   V / V_konvexe-Hülle (vorher fälschlich "Kompaktheit" genannt;
                   1.0 = konvex, < 1 = konkav/lückenhaft)
    - Elongation:  Verhältnis der Hauptachsenlängen (PCA), ≥ 1
    - Bounding Box, Äquivalentdurchmesser, max. 3D-Durchmesser, #Komponenten
    """
    all_pts = contours_to_points(contours)
    metrics = {
        "sphericity": 0.0,
        "solidity": 0.0,
        "elongation": 0.0,
        "bbox_mm": (0, 0, 0, 0, 0, 0),
        "bbox_size_mm": (0, 0, 0),
        "equivalent_diameter_mm": 0.0,
        "max_diameter_mm": 0.0,
        "volume_voxel_cm3": 0.0,
        "n_components": 0,
        "shape_valid": False,
    }

    if len(all_pts) < 4 or volume_cm3 <= 0:
        return metrics

    # Bounding Box (achsenparallel)
    mins = all_pts.min(axis=0)
    maxs = all_pts.max(axis=0)
    metrics["bbox_mm"] = tuple(np.round(np.concatenate([mins, maxs]), 2))
    metrics["bbox_size_mm"] = tuple(np.round(maxs - mins, 2))

    volume_mm3 = volume_cm3 * 1000.0
    # Äquivalent-Kugeldurchmesser aus dem (planimetrischen) Volumen
    metrics["equivalent_diameter_mm"] = round((6.0 * volume_mm3 / np.pi) ** (1 / 3), 2)

    # Max. 3D-Durchmesser (größter paarweiser Abstand = Durchmesser der Hülle)
    try:
        hull = ConvexHull(all_pts)
        hv = all_pts[hull.vertices]
        metrics["max_diameter_mm"] = round(float(pdist(hv).max()), 2)
    except Exception as e:
        _note("warning", "ANA.HULL_FAILED", "konvexe Huelle nicht bestimmbar; kein max. Durchmesser", e)

    # Beschränkte Sphärizität/Solidität + #Komponenten aus konsistenter Voxelmaske
    voxel = _voxel_shape_metrics(contours, all_pts)
    if voxel is not None:
        for k in ("sphericity", "solidity", "volume_voxel_cm3", "n_components"):
            if k in voxel:
                metrics[k] = voxel[k]
        metrics["shape_valid"] = (
            "sphericity" in voxel and "solidity" in voxel
            and voxel.get("n_components", 1) == 1
        )
    else:
        # Fallback: hüllenkonsistente Sphärizität (≤ 1), geklippte Solidität
        try:
            hull = ConvexHull(all_pts)
            if hull.area > 0:
                sph = (np.pi ** (1 / 3) * (6 * hull.volume) ** (2 / 3)) / hull.area
                metrics["sphericity"] = round(float(min(sph, 1.0)), 4)
            if hull.volume > 0:
                metrics["solidity"] = round(float(min(volume_mm3 / hull.volume, 1.0)), 4)
            metrics["n_components"] = 1
        except Exception as e:
            _note("warning", "ANA.HULL_FAILED", "konvexe Huelle nicht bestimmbar; keine Ersatz-Formmetriken", e)

    # Elongation via PCA (Verhältnis größte/kleinste Hauptachse)
    try:
        centered = all_pts - all_pts.mean(axis=0)
        cov = np.cov(centered.T)
        eigenvalues = np.sort(np.linalg.eigvalsh(cov))[::-1]
        if eigenvalues[-1] > 0:
            metrics["elongation"] = round(
                np.sqrt(eigenvalues[0] / eigenvalues[-1]), 4
            )
    except Exception as e:
        _note("warning", "ANA.ELONGATION_FAILED", "Elongation (Hauptachsen) nicht bestimmbar", e)

    return metrics


# ---------------------------------------------------------------------------
# 6. Abstandsberechnungen
# ---------------------------------------------------------------------------

def _cap_points(pts: np.ndarray, cap: int) -> np.ndarray:
    """Deterministisches Ausdünnen NUR bei sehr großen Wolken (fester Seed).

    Die alte Implementierung subsamplete jede Wolke > paar-tausend Punkte mit
    *unverseedetem* np.random.choice. Das machte (a) den Minimalabstand
    nicht-reproduzierbar und (b) verzerrte ihn nach OBEN (das Weglassen von
    Punkten kann das nächste Paar entfernen, nie ein näheres erzeugen) – also
    eine zu große, *unsichere* Abstandsangabe für die OAR-Schonung. Jetzt
    exakt; nur bei > cap Punkten wird mit festem Seed ausgedünnt.
    """
    if len(pts) > cap:
        _note("info", "ANA.POINTS_CAPPED",
              f"Punktwolke fuer Abstaende von {len(pts)} auf {cap} Punkte ausgeduennt (angenaehert)")
        idx = _rng().choice(len(pts), cap, replace=False)
        return pts[idx]
    return pts


def min_distance(pts_a: np.ndarray, pts_b: np.ndarray,
                 cap: int = 50000) -> float:
    """Exakter minimaler Abstand zwischen zwei Punktwolken in mm.

    Volle KD-Baum-Abfrage (kein verlustbehaftetes Subsampling); nur jenseits
    von ``cap`` Punkten wird deterministisch ausgedünnt.
    """
    if len(pts_a) == 0 or len(pts_b) == 0:
        return float("inf")

    pts_a = _cap_points(pts_a, cap)
    pts_b = _cap_points(pts_b, cap)

    from scipy.spatial import cKDTree
    tree = cKDTree(pts_b)
    dists, _ = tree.query(pts_a, k=1)
    return float(np.min(dists))


def hausdorff_distance(pts_a: np.ndarray, pts_b: np.ndarray,
                       cap: int = 50000) -> float:
    """Exakter (symmetrischer Maximum-)Hausdorff-Abstand in mm.

    Hinweis: der rohe Maximum-Hausdorff ist ausreißerdominiert; robustere
    Varianten (HD95, ASSD) sind klinischer Standard für die Übereinstimmung
    derselben Struktur (siehe HD95/ASSD-Metriken). Hier deterministisch.
    """
    if len(pts_a) == 0 or len(pts_b) == 0:
        return float("inf")

    pts_a = _cap_points(pts_a, cap)
    pts_b = _cap_points(pts_b, cap)

    d1 = directed_hausdorff(pts_a, pts_b)[0]
    d2 = directed_hausdorff(pts_b, pts_a)[0]
    return float(max(d1, d2))


def centroid_distance(c1: np.ndarray, c2: np.ndarray) -> float:
    """Euklidischer Abstand zwischen zwei Schwerpunkten in mm."""
    return float(np.linalg.norm(c1 - c2))


def pair_distances(pts_a: np.ndarray, pts_b: np.ndarray,
                   cap: int = 50000) -> dict:
    """Alle Punktwolken-Abstandsmaße aus *einem* Paar KD-Baum-Abfragen.

    Liefert ``min`` (nächste Annäherung), ``hausdorff`` (Maximum, ausreißer-
    empfindlich), ``hd95`` (95. Perzentil = robust) und ``assd`` (mittlerer
    symmetrischer Oberflächenabstand) – alle deterministisch in mm.
    """
    from scipy.spatial import cKDTree
    empty = {"min": float("inf"), "hausdorff": float("inf"),
             "hd95": float("inf"), "assd": float("inf")}
    if len(pts_a) == 0 or len(pts_b) == 0:
        return empty
    a = _cap_points(pts_a, cap)
    b = _cap_points(pts_b, cap)
    da, _ = cKDTree(b).query(a, k=1)   # für jeden a-Punkt der nächste in b
    db, _ = cKDTree(a).query(b, k=1)
    return {
        "min": float(min(da.min(), db.min())),
        "hausdorff": float(max(da.max(), db.max())),
        "hd95": float(max(np.percentile(da, 95), np.percentile(db, 95))),
        "assd": float((da.sum() + db.sum()) / (len(da) + len(db))),
    }


# ---------------------------------------------------------------------------
# 7. Gesamtanalyse
# ---------------------------------------------------------------------------

def get_structure_geom_types(ds: pydicom.Dataset) -> dict[int, set]:
    """Gibt {ROI-Nummer: Menge der ContourGeometricType-Werte} zurück."""
    out: dict[int, set] = {}
    if not hasattr(ds, "ROIContourSequence"):
        return out
    for rc in ds.ROIContourSequence:
        num = int(getattr(rc, "ReferencedROINumber", -1))
        g = set()
        for c in getattr(rc, "ContourSequence", []):
            g.add(str(getattr(c, "ContourGeometricType", "")))
        out[num] = g
    return out


def analyze_structure(ds: pydicom.Dataset, roi_number: int, roi_name: str,
                      category: str = "", oar_subtype: Optional[str] = None) -> dict:
    """Vollständige Analyse einer einzelnen Struktur."""
    contours = extract_contours(ds, roi_number)
    all_pts = contours_to_points(contours)
    slices = _group_slices(contours)
    dz, n_gaps = _nominal_slice_spacing(z for z, _ in slices)
    geoms = _slice_geometries(contours)          # XOR pro Schicht, einmal berechnet
    volume = compute_volume(contours, geoms=geoms)
    centroid = compute_centroid(contours, geoms=geoms)
    shape = compute_shape_metrics(contours, volume)

    return {
        "roi_number": roi_number,
        "name": roi_name,
        "category": category,
        "oar_subtype": oar_subtype,          # "serial" | "parallel" | None
        "lesion_key": lesion_key(roi_name),  # zum Paaren von GTV mit PTV
        "num_contours": len(contours),
        "num_slices": len(slices),
        "slice_spacing_mm": None if dz is None else round(dz, 3),
        "n_gaps": n_gaps,                    # z-Lücken (nie überbrückt)
        "num_points": len(all_pts),
        "volume_cm3": round(volume, 3),
        "centroid_mm": tuple(np.round(centroid, 2)),
        "shape": shape,
        "contours": contours,       # für Abstandsberechnung
        "all_points": all_pts,
    }


def _load(source) -> pydicom.Dataset:
    return source if isinstance(source, pydicom.Dataset) else load_rtstruct(str(source))


def analyze_rtstruct(source, target_names: Optional[list[str]] = None,
                     oar_names: Optional[list[str]] = None) -> tuple:
    """
    Stiller Kern von ``run_analysis`` (druckt nichts).

    ``source`` ist ein Pfad oder ein geladenes RTSTRUCT.  Rueckgabe
    ``(results, info)``: ``results`` genau wie ``run_analysis`` (Grundlage des
    JSON), ``info`` mit ``names``, ``types``, ``categories``, den gewaehlten
    ROI-Nummern (``target_nums``, ``oar_nums``, ``helper_nums``), der
    Analysereihenfolge ``analyzed`` (``[(abschnitt, ergebnis), ...]``, auch bei
    doppelten ROI-Namen vollstaendig) und ``issues``: still behandelte
    Geometriefehler und ausgeduennte Punktwolken als ``Issue``; die Werte
    selbst aendern sich dadurch nicht.  Die Zufallszahlen starten je Aufruf
    neu (Seed 0), wie in einem frischen Prozess.
    """
    ds = _load(source)
    names = get_structure_names(ds)
    types = get_structure_type(ds)
    issues: list = []
    tokens = (_RNG_VAR.set(np.random.default_rng(0)), _ISSUES_VAR.set(issues))
    try:
        # Strukturen klassifizieren (kategoriebasiert; Name-Override per --targets/--oars)
        geom_types = get_structure_geom_types(ds)
        categories = {
            num: classify_structure(nm, types.get(num, ""), geom_types.get(num, set()))
            for num, nm in names.items()
        }

        def _name_matches(nm, patterns):
            return any(p.lower() in nm.lower() for p in patterns)

        if target_names:
            target_nums = {n for n, nm in names.items()
                           if _name_matches(nm, target_names)
                           and categories[n] not in (CAT_MARKER, CAT_EXTERNAL)}
        else:
            target_nums = {n for n, c in categories.items() if c == CAT_TARGET}

        if oar_names:
            oar_nums = {n for n, nm in names.items()
                        if _name_matches(nm, oar_names)
                        and categories[n] not in (CAT_MARKER, CAT_EXTERNAL)}
        else:
            oar_nums = {n for n, c in categories.items()
                        if c in (CAT_OAR_SERIAL, CAT_OAR_PARALLEL)}

        helper_nums = {n for n, c in categories.items()
                       if c == CAT_HELPER and n not in target_nums and n not in oar_nums}

        def _subtype(num):
            c = categories[num]
            return "serial" if c == CAT_OAR_SERIAL else ("parallel" if c == CAT_OAR_PARALLEL else None)

        results = {"targets": {}, "oars": {}, "helpers": {}, "distances": [], "meta": {}}
        analyzed = []
        for section, nums in (("targets", target_nums), ("oars", oar_nums), ("helpers", helper_nums)):
            for roi_num in sorted(nums):
                roi_name = names[roi_num]
                with _roi_context(roi_name):
                    if section == "oars":
                        r = analyze_structure(ds, roi_num, roi_name, category=categories[roi_num],
                                              oar_subtype=_subtype(roi_num))
                    else:
                        r = analyze_structure(ds, roi_num, roi_name,
                                              category=categories[roi_num] if section == "targets"
                                              else CAT_HELPER)
                results[section][roi_name] = r
                analyzed.append((section, r))

        # Meta-Informationen (Kategorie-Zählungen, Marker/External nur vermerken)
        from collections import Counter
        results["meta"] = {
            "category_counts": dict(Counter(categories.values())),
            "external_names": [names[n] for n, c in categories.items() if c == CAT_EXTERNAL],
            "marker_count": sum(1 for c in categories.values() if c == CAT_MARKER),
        }

        # --------------------------------------------------------------
        # Abstände: klinisch relevante Paare (Target↔OAR + GTV↔PTV derselben
        # Läsion) statt aller C(n,2)-Kombinationen inkl. Containment-Artefakte.
        # --------------------------------------------------------------
        target_items = list(results["targets"].items())
        oar_items = list(results["oars"].items())

        pairs = []  # (name_a, ra, name_b, rb, pair_type)
        for tn, tr in target_items:
            for on, orr in oar_items:
                pairs.append((tn, tr, on, orr, "target-oar"))
        # GTV↔PTV derselben Läsion (Margin-Check)
        by_key: dict = {}
        for tn, tr in target_items:
            key = tr.get("lesion_key")
            if key:
                by_key.setdefault(key, {})[tn.split("_")[0].upper()] = (tn, tr)
        for key, d in by_key.items():
            if "GTV" in d and "PTV" in d:
                (an, ar), (bn, br) = d["GTV"], d["PTV"]
                pairs.append((an, ar, bn, br, "gtv-ptv"))

        for name_a, ra, name_b, rb, ptype in pairs:
            with _roi_context(f"{name_a} <-> {name_b}"):
                d = pair_distances(ra["all_points"], rb["all_points"])
            results["distances"].append({
                "structure_a": name_a,
                "structure_b": name_b,
                "category_a": ra.get("category"),
                "category_b": rb.get("category"),
                "oar_subtype": rb.get("oar_subtype"),
                "pair_type": ptype,
                "min_distance_mm": round(d["min"], 2),
                "hd95_mm": round(d["hd95"], 2),
                "hausdorff_distance_mm": round(d["hausdorff"], 2),
                "assd_mm": round(d["assd"], 2),
                "centroid_distance_mm": round(centroid_distance(
                    np.array(ra["centroid_mm"]), np.array(rb["centroid_mm"])), 2),
            })
        results["distances"].sort(key=lambda e: e["min_distance_mm"])
    finally:
        _RNG_VAR.reset(tokens[0])
        _ISSUES_VAR.reset(tokens[1])

    info = {"names": names, "types": types, "categories": categories,
            "target_nums": target_nums, "oar_nums": oar_nums, "helper_nums": helper_nums,
            "analyzed": analyzed, "issues": issues}
    return results, info


def run_analysis(filepath: str,
                 target_names: Optional[list[str]] = None,
                 oar_names: Optional[list[str]] = None,
                 list_only: bool = False) -> dict:
    """
    Hauptfunktion: Lädt RTSTRUCT, analysiert Strukturen, berechnet Abstände.

    Druckschicht um ``analyze_rtstruct`` (gleiche Ausgabe wie bisher).

    Parameters
    ----------
    filepath : Pfad zur RTSTRUCT DICOM Datei
    target_names : Liste der Zielgebiet-Namen (z.B. ["PTV", "CTV", "GTV"])
                   Wenn None, werden alle als "TV" typisierten Strukturen verwendet.
    oar_names : Liste der Risikoorgan-Namen
                Wenn None, werden alle als "OAR" typisierten verwendet.
    list_only : Nur Strukturnamen auflisten

    Returns
    -------
    Dictionary mit allen Analyseergebnissen
    """
    ds = load_rtstruct(filepath)
    names = get_structure_names(ds)
    types = get_structure_type(ds)

    print(f"\n{'=' * 60}")
    print(f"RTSTRUCT Analyse: {Path(filepath).name}")
    print(f"Patient: {getattr(ds, 'PatientName', 'N/A')}")
    print(f"Studie:  {getattr(ds, 'StudyDescription', 'N/A')}")
    print(f"Anzahl Strukturen: {len(names)}")
    print(f"{'=' * 60}")

    # Alle Strukturen auflisten
    print("\nVerfügbare Strukturen:")
    print(f"{'Nr':<6} {'Name':<30} {'Typ':<15}")
    print("-" * 51)
    for roi_num, roi_name in sorted(names.items()):
        rt_type = types.get(roi_num, "—")
        print(f"{roi_num:<6} {roi_name:<30} {rt_type:<15}")

    if list_only:
        return {"structures": names, "types": types}

    results, info = analyze_rtstruct(ds, target_names, oar_names)

    if not info["target_nums"]:
        print("\n(!) Keine Zielgebiete gefunden. Verwende --targets um Namen anzugeben.")
    if not info["oar_nums"]:
        print("(!) Keine Risikoorgane gefunden. Verwende --oars um Namen anzugeben.")

    headers = {
        "targets": "ZIELGEBIETE",
        "oars": "RISIKOORGANE",
        "helpers": "HILFS-/PLANUNGSSTRUKTUREN  (Formmetriken nur eingeschränkt aussagekräftig)",
    }
    for section in ("targets", "oars", "helpers"):
        if section == "helpers" and not info["helper_nums"]:
            continue
        print(f"\n{'=' * 60}")
        print(headers[section])
        print(f"{'=' * 60}")
        for sec, r in info["analyzed"]:
            if sec == section:
                _print_structure(r)

    if results["distances"]:
        print(f"\n{'=' * 60}")
        print("ABSTÄNDE  (Target<->OAR und GTV<->PTV, aufsteigend nach Min-Abstand)")
        print(f"{'=' * 60}")
        print(f"{'Struktur A':<24} {'Struktur B':<20} {'Min':>7} {'HD95':>7} "
              f"{'Haus':>7} {'ASSD':>7} {'Zentr':>7}")
        print("-" * 84)
        for e in results["distances"][:30]:
            print(f"{e['structure_a'][:23]:<24} {e['structure_b'][:19]:<20} "
                  f"{e['min_distance_mm']:>7.2f} {e['hd95_mm']:>7.2f} "
                  f"{e['hausdorff_distance_mm']:>7.2f} {e['assd_mm']:>7.2f} "
                  f"{e['centroid_distance_mm']:>7.2f}")

    return results


def inspect_rtstruct(source, volumes: bool = True) -> dict:
    """
    Schnelle Uebersicht fuer eine Oberflaeche, ohne Formmetriken und Abstaende
    und ohne Ausgabe: ROI-Tabelle (Nummer, Name, DICOM-Typ, Kategorie,
    Konturtypen, Konturzahl, planimetrisches Volumen, Anzeigefarbe, Marker),
    POINT-Marker mit Position, referenzierte FrameOfReferenceUIDs sowie Label,
    Name und Datum des Structure Sets.  ``volumes=False`` spart die
    Volumenberechnung.
    """
    ds = _load(source)
    names = get_structure_names(ds)
    types = get_structure_type(ds)
    geoms = get_structure_geom_types(ds)
    colors = {}
    for rc in ds.get("ROIContourSequence", []):
        col = rc.get("ROIDisplayColor")
        if col is not None and len(col) == 3:
            colors[int(rc.ReferencedROINumber)] = [int(v) for v in col]
    rois = []
    for num, name in sorted(names.items()):
        g = geoms.get(num, set())
        category = classify_structure(name, types.get(num, ""), g)
        contours = extract_contours(ds, num)
        rois.append({
            "number": int(num), "name": name, "rt_type": types.get(num, ""), "category": category,
            "geometric_types": sorted(t for t in g if t),
            "n_contours": len(contours),
            "volume_cm3": round(compute_volume(contours), 3) if volumes else None,
            "color": colors.get(num),
            "is_marker": category == CAT_MARKER,
        })
    return {
        "file": None if isinstance(source, pydicom.Dataset) else Path(source).name,
        "structure_set_label": str(ds.get("StructureSetLabel", "")),
        "structure_set_name": str(ds.get("StructureSetName", "")),
        "structure_set_date": str(ds.get("StructureSetDate", "")),
        "frame_of_reference_uids": sorted(get_rs_frame_of_references(ds)),
        "rois": rois,
        "markers": [{"name": n, "position_mm": [float(v) for v in p]} for n, p in find_point_markers(ds)],
    }


def _print_structure(r: dict):
    """Gibt die Analyseergebnisse einer Struktur formatiert aus."""
    s = r["shape"]
    print(f"\n  > {r['name']} (ROI #{r['roi_number']})")
    print(f"    Konturen: {r['num_contours']} auf "
          f"{r.get('num_slices', r['num_contours'])} Schichten, "
          f"{r['num_points']} Punkte")
    dz = r.get("slice_spacing_mm")
    gaps = r.get("n_gaps", 0) or 0
    dz_txt = f"{dz:.2f} mm" if dz is not None else "n/a"
    gap_txt = f"   (!) z-Lücken: {gaps}" if gaps else ""
    print(f"    Schichtabstand: {dz_txt}{gap_txt}")
    print(f"    Volumen:        {r['volume_cm3']:.3f} cm³")
    print(f"    Schwerpunkt:    x={r['centroid_mm'][0]:.1f}, "
          f"y={r['centroid_mm'][1]:.1f}, z={r['centroid_mm'][2]:.1f} mm")
    print(f"    Bounding Box:   {s['bbox_size_mm'][0]:.1f} × "
          f"{s['bbox_size_mm'][1]:.1f} × {s['bbox_size_mm'][2]:.1f} mm")
    print(f"    Äquiv.-Durchm.: {s['equivalent_diameter_mm']:.1f} mm   "
          f"(max. 3D-Durchm.: {s['max_diameter_mm']:.1f} mm)")
    valid = "" if s.get("shape_valid", False) else "  (!) (Mehrkomponenten/ungueltig)"
    print(f"    Sphärizität:    {s['sphericity']:.4f}  (1.0 = Kugel){valid}")
    print(f"    Solidität:      {s['solidity']:.4f}  "
          f"(Vol / konvexe Hülle, 1.0 = konvex)")
    print(f"    Elongation:     {s['elongation']:.4f}  "
          f"(Hauptachsen-Verhältnis)")
    if s.get("n_components", 1) and s["n_components"] > 1:
        print(f"    Komponenten:    {s['n_components']}  (Vereinigung -> "
              f"Hüllen-Formmetriken nicht aussagekräftig)")


# ---------------------------------------------------------------------------
# 7b. Self-Test (synthetische Geometrie, kein DICOM nötig)
# ---------------------------------------------------------------------------

def _synthetic_circle(z: float, r: float, cx: float = 0.0, cy: float = 0.0,
                      n: int = 64) -> np.ndarray:
    """Regelmäßiges n-Eck (Kreisapproximation) als (n,3)-Kontur bei Höhe z."""
    t = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    return np.column_stack([cx + r * np.cos(t), cy + r * np.sin(t),
                            np.full(n, float(z))])


def _synthetic_keyhole(z: float, r_out: float, r_in: float, n: int = 64,
                       delta: float = 1e-3) -> np.ndarray:
    """Ring als EINE Kontur in Keyhole-Technik: Außenbogen gegen den Uhrzeiger-
    sinn, schmaler Kanal (Breite ~2*delta), Innenbogen im Uhrzeigersinn zurück."""
    t = np.linspace(delta, 2.0 * np.pi - delta, n)
    outer = np.column_stack([r_out * np.cos(t), r_out * np.sin(t)])
    inner = np.column_stack([r_in * np.cos(t[::-1]), r_in * np.sin(t[::-1])])
    xy = np.vstack([outer, inner])
    return np.column_stack([xy, np.full(len(xy), float(z))])


def _run_self_test() -> int:
    """Konsistenztests für XOR-Löcher und z-Lücken.  Rückgabe 0 = PASS, 1 = FAIL."""
    results = []

    def check(name, ok, detail=""):
        results.append(bool(ok))
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))

    def close(a, b, rel=1e-9):
        return abs(a - b) <= rel * max(abs(a), abs(b), 1e-12)

    def sphere(z_center, R=10.0):
        return [_synthetic_circle(z_center + z, np.sqrt(R * R - z * z))
                for z in np.arange(-9.0, 9.5, 1.0)]

    def slice_info(contours):
        sl = _group_slices(contours)
        dz, g = _nominal_slice_spacing(z for z, _ in sl)
        return len(sl), dz, g

    print("\nAnalyzer Self-Test (synthetische Konturen, XOR-Löcher und z-Lücken)")
    print("-" * 70)

    # 1) Nominale Schichtdicke / Lücken
    dz, g = _nominal_slice_spacing([5.0])
    check("Schichtdicke: eine Schicht -> (None, 0)", dz is None and g == 0)
    dz, g = _nominal_slice_spacing([0.0, 3.0])
    check("Schichtdicke: [0,3] -> (3, 0 Lücken)", close(dz, 3.0) and g == 0)
    dz, g = _nominal_slice_spacing([0.0, 1.0, 32.0])
    check("Schichtdicke: [0,1,32] -> (1, 1 Lücke)", close(dz, 1.0) and g == 1,
          f"dz={dz}, gaps={g}")
    dz, g = _nominal_slice_spacing([0.0, 0.999, 2.0])
    check("Schichtdicke: Rundungsrauschen [0,0.999,2] -> ~1, 0 Lücken",
          abs(dz - 1.0) < 2e-3 and g == 0, f"dz={dz}")
    dz, g = _nominal_slice_spacing([0, 1, 2, 3, 4, 7, 10, 13])
    check("Schichtdicke: [0..4,7,10,13] -> (1, 3 Lücken)", close(dz, 1.0) and g == 3,
          f"dz={dz}, gaps={g}")

    # 2) Schichtgruppierung
    grp = _group_slices([_synthetic_circle(0.0, 5), _synthetic_circle(0.02, 5),
                         _synthetic_circle(1.0, 5)])
    check("Schichtgruppierung: z=0/0.02/1.0 -> 2 Schichten, Repräsentant 0.01",
          len(grp) == 2 and abs(grp[0][0] - 0.01) < 1e-9)

    Z = np.arange(20, dtype=float)          # 20 Schichten à 1 mm
    A20 = polygon_area(_synthetic_circle(0, 20)[:, :2])
    A10 = polygon_area(_synthetic_circle(0, 10)[:, :2])
    A4 = polygon_area(_synthetic_circle(0, 4)[:, :2])

    # 3) Volle Scheibe (Referenz)
    disk = [_synthetic_circle(z, 20) for z in Z]
    v_disk = compute_volume(disk)
    s_disk = compute_shape_metrics(disk, v_disk)
    c_disk = compute_centroid(disk)
    check("Scheibe r=20: V = A20*20 exakt", close(v_disk, A20 * 20 / 1000),
          f"V={v_disk:.4f}")
    check("Scheibe: Solidität 1.0, 1 Komponente",
          s_disk["solidity"] == 1.0 and s_disk["n_components"] == 1,
          f"sol={s_disk['solidity']}, nc={s_disk['n_components']}")
    check("Scheibe: Schwerpunkt (0,0,9.5)", np.allclose(c_disk, [0, 0, 9.5], atol=1e-6),
          f"c={np.round(c_disk, 4)}")

    # 4) Verschachtelter Ring (Loch)
    ring = [c for z in Z for c in (_synthetic_circle(z, 20), _synthetic_circle(z, 10))]
    v_ring = compute_volume(ring)
    s_ring = compute_shape_metrics(ring, v_ring)
    v_ring_exp = (A20 - A10) * 20 / 1000
    check("Ring 20/10: V = (A20-A10)*20 exakt (Loch subtrahiert)",
          close(v_ring, v_ring_exp), f"V={v_ring:.4f}, erwartet {v_ring_exp:.4f}")
    check("Ring: |V - pi*300*20/1000| < 0.5 %",
          abs(v_ring - np.pi * 300 * 20 / 1000) < 0.005 * v_ring)
    check("Ring: 1 Komponente, Solidität < 0.85",
          s_ring["n_components"] == 1 and s_ring["solidity"] < 0.85,
          f"sol={s_ring['solidity']}, nc={s_ring['n_components']}")
    check("Ring: Sphärizität < Scheibe - 0.1 (Innenfläche zählt)",
          s_ring["sphericity"] < s_disk["sphericity"] - 0.1,
          f"{s_ring['sphericity']} vs {s_disk['sphericity']}")
    check("Ring: Voxelvolumen innerhalb 1 % des planimetrischen",
          abs(s_ring["volume_voxel_cm3"] - v_ring) < 0.01 * v_ring,
          f"Vvox={s_ring['volume_voxel_cm3']}")
    check("Ring: reihenfolgeunabhängig", close(compute_volume(ring[::-1]), v_ring))

    # 5) Keyhole-Ring (eine Kontur)
    kh = [_synthetic_keyhole(z, 20, 10) for z in Z]
    v_kh = compute_volume(kh)
    s_kh = compute_shape_metrics(kh, v_kh)
    check("Keyhole-Ring: V wie verschachtelter Ring (rel. 1e-3)",
          abs(v_kh - v_ring) < 1e-3 * v_ring, f"V={v_kh:.4f}")
    check("Keyhole-Ring: 1 Komponente, Solidität < 0.85",
          s_kh["n_components"] == 1 and s_kh["solidity"] < 0.85, f"sol={s_kh['solidity']}")

    # 6) Insel im Loch
    isl = [c for z in Z for c in (_synthetic_circle(z, 20), _synthetic_circle(z, 10),
                                  _synthetic_circle(z, 4))]
    v_isl = compute_volume(isl)
    s_isl = compute_shape_metrics(isl, v_isl)
    check("Insel im Loch 20/10/4: V = (A20-A10+A4)*20 exakt",
          close(v_isl, (A20 - A10 + A4) * 20 / 1000), f"V={v_isl:.4f}")
    check("Insel im Loch: 2 Komponenten", s_isl["n_components"] == 2,
          f"nc={s_isl['n_components']}")

    # 7) Kugel
    sph = sphere(0.0)
    areas = sum(polygon_area(c[:, :2]) for c in sph)
    v_sph = compute_volume(sph)
    s_sph = compute_shape_metrics(sph, v_sph)
    _, _, g_sph = slice_info(sph)
    check("Kugel R=10: V = Summe der Flächen exakt", close(v_sph, areas / 1000),
          f"V={v_sph:.4f}")
    check("Kugel: 0 Lücken, 1 Komponente, Sphärizität > 0.8",
          g_sph == 0 and s_sph["n_components"] == 1 and s_sph["sphericity"] > 0.8,
          f"sph={s_sph['sphericity']}")

    # 8) Zwei Kugeln mit 60 mm Abstand in einer ROI
    two = sphere(0.0) + sphere(60.0)
    v_two = compute_volume(two)
    s_two = compute_shape_metrics(two, v_two)
    c_two = compute_centroid(two)
    n_two, dz_two, g_two = slice_info(two)
    check("Zwei Kugeln (60 mm Lücke): V = 2 x Kugel exakt (keine Aufblähung)",
          close(v_two, 2 * v_sph), f"V={v_two:.4f}")
    check("Zwei Kugeln: 2 Komponenten, 1 Lücke, 38 Schichten, dz=1",
          s_two["n_components"] == 2 and g_two == 1 and n_two == 38 and close(dz_two, 1.0),
          f"nc={s_two['n_components']}, gaps={g_two}, slices={n_two}, dz={dz_two}")
    check("Zwei Kugeln: Schwerpunkt z = 30", abs(c_two[2] - 30.0) < 1e-6,
          f"z={c_two[2]:.4f}")

    # 9) Kugel mit fehlender Schicht (keine Überbrückung)
    holey = [c for c in sphere(0.0) if abs(c[0, 2]) > 1e-9]
    areas_h = sum(polygon_area(c[:, :2]) for c in holey)
    v_h = compute_volume(holey)
    s_h = compute_shape_metrics(holey, v_h)
    _, _, g_h = slice_info(holey)
    check("Kugel ohne z=0: V = Summe der 18 Flächen (Lücke nicht überbrückt)",
          close(v_h, areas_h / 1000), f"V={v_h:.4f}, mit Schicht {v_sph:.4f}")
    check("Kugel ohne z=0: 1 Lücke, 2 Komponenten",
          g_h == 1 and s_h["n_components"] == 2, f"gaps={g_h}, nc={s_h['n_components']}")

    # 10) Doppelte identische Kontur -> XOR löscht aus
    dup = [c for z in Z for c in (_synthetic_circle(z, 20), _synthetic_circle(z, 20))]
    v_dup = compute_volume(dup)
    s_dup = compute_shape_metrics(dup, v_dup)
    check("Doppelkontur: V = 0 (XOR-Semantik), shape_valid False",
          v_dup == 0.0 and not s_dup["shape_valid"], f"V={v_dup}")

    n_fail = results.count(False)
    print("-" * 70)
    print(f"Gesamt: {'PASS' if n_fail == 0 else 'FAIL'}  "
          f"({len(results) - n_fail}/{len(results)} Prüfungen bestanden)")
    return 0 if n_fail == 0 else 1


# ---------------------------------------------------------------------------
# 8. CLI
# ---------------------------------------------------------------------------

def main(argv: Optional[list] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Analyse von DICOM RT Structure Sets",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Beispiele:
  %(prog)s rtstruct.dcm --list
  %(prog)s rtstruct.dcm --targets PTV,CTV,GTV --oars Parotis,Rueckenmark,Blase
  %(prog)s rtstruct.dcm   # Auto-Erkennung über DICOM RT ROI Type
  %(prog)s --self-test    # Synthetische Konsistenztests (Löcher, z-Lücken)
        """,
    )
    parser.add_argument("file", nargs="?", default=None,
                        help="Pfad zur RTSTRUCT DICOM Datei")
    parser.add_argument("--list", action="store_true",
                        help="Nur Strukturnamen auflisten")
    parser.add_argument("--targets", type=str, default=None,
                        help="Komma-getrennte Zielgebiet-Namen (z.B. PTV,CTV)")
    parser.add_argument("--oars", type=str, default=None,
                        help="Komma-getrennte Risikoorgan-Namen")
    parser.add_argument("--output", "-o", type=str, default=None,
                        help="Ausgabeverzeichnis für <stem>_analysis.json")
    parser.add_argument("--self-test", action="store_true",
                        help="Synthetische Konsistenztests (Ring/Loch, Keyhole, "
                             "Kugeln, z-Lücken); Exit 0 = pass, 1 = fail")

    args = parser.parse_args(argv)

    if args.self_test:
        return _run_self_test()
    if not args.file:
        parser.error("file ist erforderlich (außer mit --self-test)")

    target_list = parse_name_list(args.targets)
    oar_list = parse_name_list(args.oars)

    results = run_analysis(
        filepath=args.file,
        target_names=target_list,
        oar_names=oar_list,
        list_only=args.list,
    )

    if args.output:
        out_dir = Path(args.output)
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{Path(args.file).stem}_analysis.json"
        with out_path.open("w", encoding="utf-8") as fh:
            json.dump(_results_to_jsonable(results), fh,
                      indent=2, ensure_ascii=False)
        print(f"\nAnalyse-Ergebnisse gespeichert: {out_path}")

    print(f"\n{'=' * 60}")
    print("Analyse abgeschlossen.")
    print(f"{'=' * 60}\n")
    return 0


def _results_to_jsonable(obj):
    """Konvertiert run_analysis-Ergebnisse in JSON-serialisierbare Strukturen.

    Entfernt rohe Punkt-/Konturen-Arrays (zu groß und nicht JSON-tauglich)
    und wandelt NumPy-Skalare/-Arrays sowie Tuples in native Typen um.
    """
    if isinstance(obj, dict):
        return {k: _results_to_jsonable(v) for k, v in obj.items()
                if k not in ("contours", "all_points")}
    if isinstance(obj, (list, tuple)):
        return [_results_to_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    return obj


if __name__ == "__main__":
    sys.exit(main())