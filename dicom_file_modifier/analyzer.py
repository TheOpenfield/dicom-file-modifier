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

Benötigte Packages:
  pip install pydicom numpy scipy shapely matplotlib
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Optional

import numpy as np
import pydicom
from scipy.spatial import ConvexHull
from scipy.spatial.distance import directed_hausdorff, pdist
from shapely.geometry import Polygon

# Fester Zufallsgenerator: Distanz-/Subsample-Operationen sollen reproduzierbar
# sein (vorher unverseedetes np.random.choice -> nicht-deterministische QA-Werte).
_RNG = np.random.default_rng(0)


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
    r"(^h_|opt[\s_]*system|^opt_|\bring\b|\bshell\b|\+\s*\d+\s*mm)",
    re.IGNORECASE,
)
# Echte Zielvolumina: Name beginnt mit GTV/PTV/CTV/ITV (gefolgt von _ oder Ziffer)
_TARGET_NAME_RE = re.compile(r"^(gtv|ptv|ctv|itv)[ _0-9]", re.IGNORECASE)
# Läsions-Schlüssel zum Paaren von GTV mit seinem PTV (z.B. "01M_BM_frontal_li")
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
    # 3) Hilfs-/Planungsstrukturen (Name-basiert, überschreibt fehlerhaftes RT-Type)
    if _HELPER_NAME_RE.search(nm):
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
# 3. Volumenberechnung
# ---------------------------------------------------------------------------

def polygon_area(pts_2d: np.ndarray) -> float:
    """Berechnet die Fläche eines 2D-Polygons mit der Shoelace-Formel."""
    try:
        poly = Polygon(pts_2d)
        if not poly.is_valid:
            poly = poly.buffer(0)
        return abs(poly.area)
    except Exception:
        return 0.0


def compute_volume(contours: list[np.ndarray]) -> float:
    """
    Berechnet das Volumen in cm³.
    Summiert Konturflächen × Schichtabstand.
    """
    if len(contours) < 2:
        return 0.0

    # Schicht-Z-Werte und Flächen sammeln
    slice_data = {}
    for pts in contours:
        z = round(pts[0, 2], 3)
        area = polygon_area(pts[:, :2])
        slice_data.setdefault(z, 0.0)
        slice_data[z] += area  # Mehrere Konturen pro Schicht addieren

    z_values = sorted(slice_data.keys())
    if len(z_values) < 2:
        return 0.0

    # Mittleren Schichtabstand bestimmen
    dz = np.mean(np.diff(z_values))

    total_area = sum(slice_data.values())
    volume_mm3 = total_area * abs(dz)
    return volume_mm3 / 1000.0  # mm³ -> cm³


# ---------------------------------------------------------------------------
# 4. Schwerpunkt
# ---------------------------------------------------------------------------

def _polygon_centroid_xy(pts_2d: np.ndarray) -> Optional[np.ndarray]:
    """Flächenschwerpunkt (Shoelace-Moment) eines 2D-Polygons, oder None.

    Der arithmetische Mittelwert der *Eckpunkte* ist NICHT der Polygon-
    schwerpunkt (er ist zu dicht besetzten Randabschnitten hin verzerrt); wir
    verwenden daher shapely ``Polygon.centroid`` (= Flächenmoment-Integral).
    """
    try:
        poly = Polygon(pts_2d)
        if not poly.is_valid:
            poly = poly.buffer(0)
        if poly.is_empty or poly.area == 0:
            return None
        c = poly.centroid
        return np.array([c.x, c.y])
    except Exception:
        return None


def compute_centroid(contours: list[np.ndarray]) -> np.ndarray:
    """
    Berechnet den flächengewichteten Schwerpunkt (x, y, z) in mm.
    Jede Schicht wird mit ihrer Konturfläche gewichtet; die In-Plane-Position
    jeder Schicht ist der echte Polygon-Flächenschwerpunkt (nicht das
    Eckpunkt-Mittel).
    """
    if not contours:
        return np.array([0.0, 0.0, 0.0])

    weighted_sum = np.zeros(3)
    total_weight = 0.0

    for pts in contours:
        area = polygon_area(pts[:, :2])
        centroid_xy = _polygon_centroid_xy(pts[:, :2])
        if centroid_xy is None or area <= 0:
            continue
        z = pts[0, 2]
        weighted_sum += area * np.array([centroid_xy[0], centroid_xy[1], z])
        total_weight += area

    if total_weight == 0:
        all_pts = contours_to_points(contours)
        return all_pts.mean(axis=0)

    return weighted_sum / total_weight


# ---------------------------------------------------------------------------
# 5. Formanalyse
# ---------------------------------------------------------------------------

def _rasterize_structure(contours: list[np.ndarray],
                         target_dim: int = 96):
    """
    Rastert die gestapelten Konturen in eine binäre 3D-Voxelmaske.

    Liefert ``(mask, spacing, centers)`` mit ``mask`` (nz, ny, nx) bool,
    ``spacing`` = (sz, sy, sx) in mm und ``centers`` als (N,3)-Array der
    Voxelmittelpunkte in LPS-mm (für den Hüllen-Halbraumtest), oder ``None``
    bei zu wenig Daten.  Eine konsistente Maske ist die Grundlage für
    mathematisch *beschränkte* Sphärizität/Solidität (im Gegensatz zur alten,
    inkonsistenten Mischung aus planimetrischem Volumen und konvexer Hülle).
    """
    try:
        from matplotlib.path import Path as MplPath
    except Exception:
        return None

    all_pts = contours_to_points(contours)
    if len(all_pts) < 4:
        return None

    x0, y0, z0 = all_pts.min(axis=0)
    x1, y1, z1 = all_pts.max(axis=0)
    span_x, span_y, span_z = x1 - x0, y1 - y0, z1 - z0
    if span_x <= 0 or span_y <= 0:
        return None

    zs = np.array(sorted({round(float(p[0, 2]), 4) for p in contours}))
    dz = float(np.mean(np.diff(zs))) if len(zs) > 1 else 1.0
    dz = abs(dz) if dz != 0 else 1.0

    s_xy = max(span_x, span_y) / target_dim
    s_xy = float(np.clip(s_xy, 0.3, 5.0))
    s_z = dz
    pad = 2

    nx = int(np.ceil(span_x / s_xy)) + 1 + 2 * pad
    ny = int(np.ceil(span_y / s_xy)) + 1 + 2 * pad
    nz = int(round(span_z / s_z)) + 1 + 2 * pad
    if nx * ny * nz > 6_000_000:          # Sicherheitskappe gegen Speicher-Spikes
        scale = (nx * ny * nz / 6_000_000) ** (1 / 3)
        s_xy *= scale
        nx = int(np.ceil(span_x / s_xy)) + 1 + 2 * pad
        ny = int(np.ceil(span_y / s_xy)) + 1 + 2 * pad

    ox, oy, oz = x0 - pad * s_xy, y0 - pad * s_xy, z0 - pad * s_z
    xc = ox + (np.arange(nx) + 0.5) * s_xy
    yc = oy + (np.arange(ny) + 0.5) * s_xy
    gx, gy = np.meshgrid(xc, yc)          # (ny, nx)
    grid_xy = np.column_stack([gx.ravel(), gy.ravel()])

    mask = np.zeros((nz, ny, nx), dtype=bool)
    for pts in contours:
        z = float(pts[0, 2])
        k = int(round((z - oz) / s_z))
        if not (0 <= k < nz):
            continue
        if len(pts) < 3:
            continue
        try:
            inside = MplPath(pts[:, :2]).contains_points(grid_xy)
        except Exception:
            continue
        mask[k] |= inside.reshape(ny, nx)

    if not mask.any():
        return None

    # Voxelmittelpunkt-Achsen (für den Hüllen-Halbraumtest in _voxel_shape_metrics)
    zc = oz + (np.arange(nz) + 0.5) * s_z
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
    except Exception:
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
    except Exception:
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
    except Exception:
        pass

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
    except Exception:
        pass

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
        except Exception:
            pass

    # Elongation via PCA (Verhältnis größte/kleinste Hauptachse)
    try:
        centered = all_pts - all_pts.mean(axis=0)
        cov = np.cov(centered.T)
        eigenvalues = np.sort(np.linalg.eigvalsh(cov))[::-1]
        if eigenvalues[-1] > 0:
            metrics["elongation"] = round(
                np.sqrt(eigenvalues[0] / eigenvalues[-1]), 4
            )
    except Exception:
        pass

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
        idx = _RNG.choice(len(pts), cap, replace=False)
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
    volume = compute_volume(contours)
    centroid = compute_centroid(contours)
    shape = compute_shape_metrics(contours, volume)

    return {
        "roi_number": roi_number,
        "name": roi_name,
        "category": category,
        "oar_subtype": oar_subtype,          # "serial" | "parallel" | None
        "lesion_key": lesion_key(roi_name),  # zum Paaren von GTV mit PTV
        "num_contours": len(contours),
        "num_points": len(all_pts),
        "volume_cm3": round(volume, 3),
        "centroid_mm": tuple(np.round(centroid, 2)),
        "shape": shape,
        "contours": contours,       # für Abstandsberechnung
        "all_points": all_pts,
    }


def run_analysis(filepath: str,
                 target_names: Optional[list[str]] = None,
                 oar_names: Optional[list[str]] = None,
                 list_only: bool = False) -> dict:
    """
    Hauptfunktion: Lädt RTSTRUCT, analysiert Strukturen, berechnet Abstände.

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

    if not target_nums:
        print("\n(!) Keine Zielgebiete gefunden. Verwende --targets um Namen anzugeben.")
    if not oar_nums:
        print("(!) Keine Risikoorgane gefunden. Verwende --oars um Namen anzugeben.")

    def _subtype(num):
        c = categories[num]
        return "serial" if c == CAT_OAR_SERIAL else ("parallel" if c == CAT_OAR_PARALLEL else None)

    results = {"targets": {}, "oars": {}, "helpers": {}, "distances": [], "meta": {}}

    print(f"\n{'=' * 60}")
    print("ZIELGEBIETE")
    print(f"{'=' * 60}")
    for roi_num in sorted(target_nums):
        roi_name = names[roi_num]
        r = analyze_structure(ds, roi_num, roi_name, category=categories[roi_num])
        results["targets"][roi_name] = r
        _print_structure(r)

    print(f"\n{'=' * 60}")
    print("RISIKOORGANE")
    print(f"{'=' * 60}")
    for roi_num in sorted(oar_nums):
        roi_name = names[roi_num]
        r = analyze_structure(ds, roi_num, roi_name,
                              category=categories[roi_num], oar_subtype=_subtype(roi_num))
        results["oars"][roi_name] = r
        _print_structure(r)

    if helper_nums:
        print(f"\n{'=' * 60}")
        print("HILFS-/PLANUNGSSTRUKTUREN  (Formmetriken nur eingeschränkt aussagekräftig)")
        print(f"{'=' * 60}")
        for roi_num in sorted(helper_nums):
            roi_name = names[roi_num]
            r = analyze_structure(ds, roi_num, roi_name, category=CAT_HELPER)
            results["helpers"][roi_name] = r
            _print_structure(r)

    # Meta-Informationen (Kategorie-Zählungen, Marker/External nur vermerken)
    from collections import Counter
    results["meta"] = {
        "category_counts": dict(Counter(categories.values())),
        "external_names": [names[n] for n, c in categories.items() if c == CAT_EXTERNAL],
        "marker_count": sum(1 for c in categories.values() if c == CAT_MARKER),
    }

    # ------------------------------------------------------------------
    # Abstände: klinisch relevante Paare (Target↔OAR + GTV↔PTV derselben Läsion)
    # statt aller C(n,2)-Kombinationen inkl. Containment-Artefakte.
    # ------------------------------------------------------------------
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

    if pairs:
        print(f"\n{'=' * 60}")
        print("ABSTÄNDE  (Target↔OAR und GTV↔PTV, aufsteigend nach Min-Abstand)")
        print(f"{'=' * 60}")
        for name_a, ra, name_b, rb, ptype in pairs:
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
        print(f"{'Struktur A':<24} {'Struktur B':<20} {'Min':>7} {'HD95':>7} "
              f"{'Haus':>7} {'ASSD':>7} {'Zentr':>7}")
        print("-" * 84)
        for e in results["distances"][:30]:
            print(f"{e['structure_a'][:23]:<24} {e['structure_b'][:19]:<20} "
                  f"{e['min_distance_mm']:>7.2f} {e['hd95_mm']:>7.2f} "
                  f"{e['hausdorff_distance_mm']:>7.2f} {e['assd_mm']:>7.2f} "
                  f"{e['centroid_distance_mm']:>7.2f}")

    return results


def _print_structure(r: dict):
    """Gibt die Analyseergebnisse einer Struktur formatiert aus."""
    s = r["shape"]
    print(f"\n  > {r['name']} (ROI #{r['roi_number']})")
    print(f"    Konturen: {r['num_contours']} Schichten, "
          f"{r['num_points']} Punkte")
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
# 8. CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Analyse von DICOM RT Structure Sets",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Beispiele:
  %(prog)s rtstruct.dcm --list
  %(prog)s rtstruct.dcm --targets PTV,CTV,GTV --oars Parotis,Rueckenmark,Blase
  %(prog)s rtstruct.dcm   # Auto-Erkennung über DICOM RT ROI Type
        """,
    )
    parser.add_argument("file", help="Pfad zur RTSTRUCT DICOM Datei")
    parser.add_argument("--list", action="store_true",
                        help="Nur Strukturnamen auflisten")
    parser.add_argument("--targets", type=str, default=None,
                        help="Komma-getrennte Zielgebiet-Namen (z.B. PTV,CTV)")
    parser.add_argument("--oars", type=str, default=None,
                        help="Komma-getrennte Risikoorgan-Namen")
    parser.add_argument("--output", "-o", type=str, default=None,
                        help="Ausgabeverzeichnis für <stem>_analysis.json")

    args = parser.parse_args()

    target_list = args.targets.split(",") if args.targets else None
    oar_list = args.oars.split(",") if args.oars else None

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
    main()