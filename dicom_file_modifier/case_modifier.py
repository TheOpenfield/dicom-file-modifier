#!/usr/bin/env python3
"""
Case Modifier  Rigid Body Transform fuer CT + RTSTRUCT im Verbund.
====================================================================

Erweitert den bestehenden ``modifier.py`` um eine Case-orientierte Sicht:
Eingabe ist ein Patienten-/Case-Ordner mit fester Struktur

    <case-dir>/
        CT/                       (alle CT-DICOM-Dateien)
        RS*.dcm                   (genau eine RTSTRUCT-Datei)
        [optional] RP*.dcm        (RTPLAN  -- wird NICHT mit-transformiert)
        [optional] RD*.dcm        (RTDOSE  -- wird NICHT mit-transformiert)

Auf diesen Datensatz wird dieselbe rigide Transformation T (Translation
+ intrinsische XYZ-Euler-Rotation, SciPy "XYZ") angewendet, die ``modifier.py`` bereits
fuer das CT bereitstellt.

Stufenweise Implementierung:

  Stage 1 (HIER):  CT wird wie bisher transformiert; RTSTRUCT wird unveraendert
                   in den Output kopiert (mit deutlicher Warnung, dass die
                   Konturen noch nicht mitwandern).  Pre-Flight: Auto-Discovery
                   und Konsistenzpruefung der FrameOfReferenceUID zwischen CT
                   und RTSTRUCT.

  Stage 2+:        Konturpunkte werden mit T transformiert, UID-Verweise im RS
                   auf das neue CT umgeschrieben, Drehpunkt-Marker eingefuegt,
                   variables Rotationszentrum, Aria-sichtbare Metadaten,
                   Robustheitschecks.

Verwendung (Stage 1):
  python -m dicom_file_modifier.case_modifier <case-dir> [Optionen]

Beispiele:
  # Identitaets-Lauf (nur UID-Refresh, gut zum Verifizieren)
  python -m dicom_file_modifier.case_modifier data/<case-id>

  # Reale Transformation (Stage 1 transformiert nur das CT)
  python -m dicom_file_modifier.case_modifier data/<case-id> \
      --tx 10 --ty 0 --tz -5 --rx 0 --ry 0 --rz 15 --output output/run1
"""

from __future__ import annotations

import argparse
import copy
import sys
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

import numpy as np
import pydicom
from pydicom.dataset import Dataset
from pydicom.sequence import Sequence
from pydicom.uid import generate_uid

from . import _runtime
from . import modifier as mod
# Seit P0.6 in modifier; hier weiter importierbar (alte Importpfade)
from .modifier import load_ct_headers, validate_ct_geometry
# Seit P0.3 in dicom_utils; hier weiter importierbar (alte Importpfade)
from .dicom_utils import (_label_with_suffix, _truncate, find_point_markers,
                          get_rs_frame_of_references)
from .issues import Issue, UserInputError

VERIFY_THRESHOLD_MM = 1e-3      # --verify: max. Centroid-Abweichung fuer PASS
PLANE_TOL_MM = 0.01             # resample: Kontur liegt in einer CT-Schichtebene


# ---------------------------------------------------------------------------
# Auto-Discovery
# ---------------------------------------------------------------------------

def discover_case(case_dir: str, rs_override: "str | None" = None, *,
                  return_siblings: bool = False) -> tuple:
    """
    Sucht im ``case_dir`` den ``CT/``-Unterordner und genau eine ``RS*.dcm``-Datei.

    Wirft ``FileNotFoundError`` / ``ValueError`` mit klarer Fehlermeldung,
    falls die Konvention verletzt ist.  Bei mehreren RS-Dateien wird der Aufrufer
    aufgefordert, mit ``--rs <pfad>`` explizit auszuwaehlen.

    Rueckgabe ``(ct_dir, rs_path)``; parallel liegende RP*/RD*-Dateien werden
    dann als Hinweis gedruckt.  Mit ``return_siblings=True`` still und
    ``(ct_dir, rs_path, siblings)`` (Liste der RP*/RD*-Pfade).
    """
    case = Path(case_dir)
    if not case.is_dir():
        raise FileNotFoundError(f"Case-Ordner nicht gefunden: {case_dir!r}")

    ct_dir = case / "CT"
    if not ct_dir.is_dir():
        raise FileNotFoundError(
            f"Erwarteter Unterordner 'CT' fehlt in {case_dir!r} "
            f"(gesucht: {ct_dir})."
        )

    if rs_override is not None:
        rs_path = Path(rs_override)
        if not rs_path.is_file():
            raise FileNotFoundError(f"RS-Datei nicht gefunden: {rs_override!r}")
    else:
        rs_files = sorted(case.glob("RS*.dcm"))
        if len(rs_files) == 0:
            raise FileNotFoundError(
                f"Keine 'RS*.dcm'-Datei in {case_dir!r} gefunden. "
                "Liegt das RTSTRUCT mit anderem Praefix? Dann --rs <pfad> nutzen."
            )
        if len(rs_files) > 1:
            joined = "\n  ".join(str(p) for p in rs_files)
            raise ValueError(
                f"Mehrere 'RS*.dcm'-Kandidaten in {case_dir!r}:\n  {joined}\n"
                "Bitte mit --rs <pfad> explizit auswaehlen."
            )
        rs_path = rs_files[0]

    extras = sorted(list(case.glob("RP*.dcm")) + list(case.glob("RD*.dcm")))
    if return_siblings:
        return ct_dir, rs_path, extras
    if extras:
        print(_siblings_note(extras))
    return ct_dir, rs_path


def _siblings_note(extras: list) -> str:
    names = ", ".join(p.name for p in extras)
    return (f"  Hinweis: Zusaetzliche Plan-/Dosis-Dateien gefunden "
            f"({names}). Diese werden NICHT mit-transformiert.")


# ---------------------------------------------------------------------------
# FrameOfReferenceUID-Validierung
# ---------------------------------------------------------------------------

def get_ct_frame_of_reference(slices: list) -> str:
    """Liest die einheitliche FrameOfReferenceUID aller CT-Slices aus."""
    for_uids = set()
    for s in slices:
        if hasattr(s, "FrameOfReferenceUID"):
            for_uids.add(str(s.FrameOfReferenceUID))
    if not for_uids:
        raise ValueError("CT-Slices besitzen keine FrameOfReferenceUID.")
    if len(for_uids) > 1:
        raise ValueError(
            "CT-Slices haben uneinheitliche FrameOfReferenceUIDs:\n  "
            + "\n  ".join(sorted(for_uids))
        )
    return for_uids.pop()


def validate_for_consistency(ct_for_uid: str, rs_ds: pydicom.Dataset) -> None:
    """Stellt sicher, dass das RTSTRUCT die FoR des CT referenziert."""
    rs_for_uids = get_rs_frame_of_references(rs_ds)
    if ct_for_uid not in rs_for_uids:
        raise ValueError(
            "FrameOfReferenceUID-Mismatch zwischen CT und RTSTRUCT.\n"
            f"  CT-FoR : {ct_for_uid}\n"
            f"  RS-FoRs: {sorted(rs_for_uids) or '(keine)'}\n"
            "Das RTSTRUCT gehoert offenbar nicht zu diesem CT."
        )


# ---------------------------------------------------------------------------
# CT-Geometrie-Validierung (Stage 5)
# ---------------------------------------------------------------------------

def check_contour_clipping(
    new_rs_ds: pydicom.Dataset,
    geom: dict,
    method: str,
) -> list:
    """
    Pruef nach der Transformation, ob Konturpunkte ausserhalb des CT-Voxelgrids
    liegen.  Liefert eine Liste ``[(roi_name, n_outside, n_total, frac), ...]``
    der betroffenen ROIs.

    Hintergrund:
      - **resample**:  Output-CT-Grid behaelt das Original-Affine A.  Konturen
        werden durch T verschoben; ein Punkt p_neu ist nur dann durch das CT
        abgedeckt, wenn ``A^-1 . p_neu`` innerhalb [0, n-1] liegt.  Andernfalls
        gibt es keine Bilddaten zur Kontur (im TPS sieht man die Struktur,
        aber kein darunterliegendes Gewebe).
      - **metadata**:  Output-Affine ist A_neu = T . A.  Da Punkte ebenfalls mit
        T transformiert wurden, gilt ``A_neu^-1 . p_neu = A^-1 . p_orig``,
        also stets im Bild.  Kein Clipping moeglich.

    Wird vom Top-Level-Workflow mit ``method`` aus der Pipeline aufgerufen.
    """
    if method != "resample":
        return []  # in metadata mode geometrisch unmoeglich
    A_inv = np.linalg.inv(geom["affine"])
    nz, ny, nx = geom["shape"]

    name_map = {
        int(r.ROINumber): str(r.ROIName)
        for r in getattr(new_rs_ds, "StructureSetROISequence", [])
    }

    issues: list[tuple[str, int, int, float]] = []
    if not hasattr(new_rs_ds, "ROIContourSequence"):
        return issues

    for rc in new_rs_ds.ROIContourSequence:
        if not hasattr(rc, "ContourSequence"):
            continue
        roi_num = int(getattr(rc, "ReferencedROINumber", -1))
        roi_name = name_map.get(roi_num, f"ROI#{roi_num}")
        # Drehpunkt-Marker selbst muss nicht gepruft werden (1 Punkt am Zentrum)
        if roi_name == "Drehpunkt":
            continue

        all_pts = []
        for c in rc.ContourSequence:
            if not hasattr(c, "ContourData"):
                continue
            pts = np.array(c.ContourData, dtype=np.float64).reshape(-1, 3)
            if pts.size:
                all_pts.append(pts)
        if not all_pts:
            continue

        pts = np.vstack(all_pts)
        n_total = pts.shape[0]
        homog = np.hstack([pts, np.ones((n_total, 1))])
        vox = (A_inv @ homog.T).T[:, :3]

        outside = (
            (vox[:, 0] < 0) | (vox[:, 0] > nz - 1) |
            (vox[:, 1] < 0) | (vox[:, 1] > ny - 1) |
            (vox[:, 2] < 0) | (vox[:, 2] > nx - 1)
        )
        n_outside = int(outside.sum())
        if n_outside > 0:
            issues.append((roi_name, n_outside, n_total, n_outside / n_total))

    return issues


def _verify_rs_centroids(orig_ds: pydicom.Dataset,
                         new_rs_path: str,
                         T: np.ndarray,
                         threshold_mm: float = VERIFY_THRESHOLD_MM) -> dict:
    """
    Liest das frisch geschriebene RTSTRUCT wieder ein, berechnet pro ROI den
    Centroid und vergleicht ihn mit ``T @ centroid_orig``.  Druckt die maximale
    Abweichung pro ROI sowie die Gesamtstatistik.

    Centroids sind unter rigiden Transformationen linear: T(mean(p_i)) =
    mean(T(p_i)).  Eine signifikante Abweichung deutet auf einen Indizierungs-
    oder Reshape-Bug im Transform-Pfad hin.

    Rueckgabe: dict mit ``"max_err_mm"``, ``"checked"``, ``"worst_roi"``,
    ``"threshold_mm"`` und ``"passed"`` (``max_err_mm <= threshold_mm``).
    Nur bei FAIL wird zusaetzlich eine Warnung gedruckt; der CLI-Exit bleibt.
    """
    from .analyzer import load_rtstruct, extract_contours, get_structure_names

    new_ds   = load_rtstruct(new_rs_path)
    names    = get_structure_names(orig_ds)
    max_err  = 0.0
    worst    = None
    checked  = 0

    print("\n--verify  Centroid-Linearitaets-Check:")
    for roi_num, roi_name in names.items():
        c_orig = extract_contours(orig_ds, roi_num)
        c_new  = extract_contours(new_ds, roi_num)
        if not c_orig or not c_new:
            continue
        pts_o = np.vstack(c_orig)
        pts_n = np.vstack(c_new)
        if pts_o.shape != pts_n.shape:
            continue
        cen_o = pts_o.mean(axis=0)
        cen_n = pts_n.mean(axis=0)
        cen_e = (T @ np.append(cen_o, 1.0))[:3]
        err = float(np.linalg.norm(cen_n - cen_e))
        checked += 1
        if err > max_err:
            max_err = err
            worst = roi_name
    print(f"  Geprueft: {checked} ROIs   "
          f"max Abweichung: {max_err:.3e} mm   "
          f"(worst: {worst})")
    passed = bool(max_err <= threshold_mm)
    if not passed:
        print(f"  ! WARNUNG ! Centroid-Abweichung ueber der Schwelle {threshold_mm:g} mm -> FAIL")
    return {"max_err_mm": max_err, "checked": checked, "worst_roi": worst,
            "threshold_mm": float(threshold_mm), "passed": passed}


# Unter Windows in Datei-/Ordnernamen verboten
_WIN_FORBIDDEN = set('<>:"/\\|?*')
_WIN_RESERVED = {"CON", "PRN", "AUX", "NUL", *(f"COM{i}" for i in range(1, 10)),
                 *(f"LPT{i}" for i in range(1, 10))}


def validate_label(label: str, case_id: str) -> None:
    """
    ``--label`` bildet den Ausgabeordner ``<case_id><label>`` und die Datei
    ``RS<label>.dcm``; Zeichen und Namen, die Windows ablehnt, ergeben einen
    ``UserInputError`` (CLI-Exit 2), bevor etwas geschrieben wird.
    """
    bad = sorted({c for c in label if c in _WIN_FORBIDDEN or ord(c) < 32})
    if bad:
        shown = " ".join(repr(c) for c in bad)
        raise UserInputError(Issue(
            "error", "CASE.LABEL_INVALID", f"--label {label!r} enthaelt unzulaessige Zeichen: {shown}",
            hint_de='Buchstaben, Ziffern, "_", "-" und "." verwenden (nicht < > : " / \\ | ? *).',
            field="label"))
    for name in (f"{case_id}{label}", f"RS{label}"):
        if name.split(".")[0].upper() in _WIN_RESERVED or name.endswith((".", " ")):
            raise UserInputError(Issue(
                "error", "CASE.LABEL_INVALID",
                f"--label {label!r} ergibt den unter Windows unzulaessigen Namen {name!r}.",
                hint_de="Reservierte Namen (CON, NUL, COM1 ...) sowie Punkt oder Leerzeichen am Ende vermeiden.",
                field="label"))


# ---------------------------------------------------------------------------
# RTSTRUCT-Transformation (Stage 2)
# ---------------------------------------------------------------------------

def _apply_T_to_flat_coords(flat_data, T: np.ndarray) -> list[str]:
    """
    Wendet T (4x4) auf eine flache DICOM-ContourData-Liste [x1,y1,z1,x2,y2,z2,...]
    an und gibt die transformierte Liste als formattierte Strings zurueck.
    """
    pts = np.array(flat_data, dtype=np.float64).reshape(-1, 3)
    if pts.size == 0:
        return []
    pts_h = np.hstack([pts, np.ones((pts.shape[0], 1))])
    new   = (T @ pts_h.T).T[:, :3]
    return [f"{v:.6f}" for v in new.flatten()]


def _rewrite_referenced_sops(seq, sop_map: dict) -> int:
    """
    Schreibt jedes ``ReferencedSOPInstanceUID`` in der gegebenen Sequenz auf den
    neuen Wert um (alt -> neu via ``sop_map``).  Wirft ``KeyError`` mit klarer
    Meldung, falls eine alte SOP nicht im Map auftaucht; das deutet darauf hin,
    dass das RS auf ein anderes CT verweist.
    """
    n = 0
    for item in seq:
        if hasattr(item, "ReferencedSOPInstanceUID"):
            old = str(item.ReferencedSOPInstanceUID)
            if old not in sop_map:
                raise KeyError(
                    "RTSTRUCT-Verweis auf unbekannte CT-SOPInstanceUID:\n"
                    f"  {old}\n"
                    "Diese SOP gehoerte nicht zum verarbeiteten CT-Verzeichnis."
                )
            item.ReferencedSOPInstanceUID = sop_map[old]
            n += 1
    return n


def _next_roi_number(rs_ds: pydicom.Dataset) -> int:
    """Naechste freie ROINumber im StructureSetROISequence."""
    if not hasattr(rs_ds, "StructureSetROISequence") or not rs_ds.StructureSetROISequence:
        return 1
    return max(int(r.ROINumber) for r in rs_ds.StructureSetROISequence) + 1


def add_drehpunkt_marker(
    rs_ds: pydicom.Dataset,
    position_lps: np.ndarray,
    for_uid: str,
    color: "tuple[int, int, int]" = (255, 255, 0),
) -> int:
    """
    Fuegt eine POINT-Type-ROI mit Namen ``Drehpunkt`` an die gegebene Position
    (LPS, mm) in das uebergebene RTSTRUCT ein.  Im Planungssystem markiert
    dieser Punkt den tatsaechlichen Rotationsmittelpunkt der Transformation
    nach Anwendung der Translation (Rotation laesst das Zentrum invariant).

    Gibt die neue ROINumber zurueck.
    """
    new_roi_num = _next_roi_number(rs_ds)

    # 1. StructureSetROISequence
    roi = Dataset()
    roi.ROINumber = new_roi_num
    roi.ReferencedFrameOfReferenceUID = for_uid
    roi.ROIName = "Drehpunkt"
    roi.ROIGenerationAlgorithm = "MANUAL"
    if not hasattr(rs_ds, "StructureSetROISequence"):
        rs_ds.StructureSetROISequence = Sequence()
    rs_ds.StructureSetROISequence.append(roi)

    # 2. RTROIObservationsSequence
    obs = Dataset()
    if not hasattr(rs_ds, "RTROIObservationsSequence"):
        rs_ds.RTROIObservationsSequence = Sequence()
    if rs_ds.RTROIObservationsSequence:
        obs.ObservationNumber = max(
            int(o.ObservationNumber) for o in rs_ds.RTROIObservationsSequence
        ) + 1
    else:
        obs.ObservationNumber = new_roi_num
    obs.ReferencedROINumber = new_roi_num
    obs.ROIObservationLabel = "Drehpunkt"
    obs.RTROIInterpretedType = "MARKER"
    obs.ROIInterpreter = ""
    rs_ds.RTROIObservationsSequence.append(obs)

    # 3. ROIContourSequence (POINT-Contour mit einem Punkt)
    rc = Dataset()
    rc.ReferencedROINumber = new_roi_num
    rc.ROIDisplayColor = list(color)

    contour = Dataset()
    contour.ContourGeometricType = "POINT"
    contour.NumberOfContourPoints = 1
    contour.ContourData = [f"{v:.6f}" for v in np.asarray(position_lps).reshape(-1)]
    rc.ContourSequence = Sequence([contour])

    if not hasattr(rs_ds, "ROIContourSequence"):
        rs_ds.ROIContourSequence = Sequence()
    rs_ds.ROIContourSequence.append(rc)

    return new_roi_num


def build_transform_description(
    tx: float, ty: float, tz: float,
    rx: float, ry: float, rz: float,
    center_label: str,
    method: str,
    for_strategy: str,
    max_len: int = 64,
) -> str:
    """
    Baut einen kompakten, menschenlesbaren Transform-Beschreibungs-String, der in
    ``StructureSetDescription`` (DICOM VR LO, max. 64 Zeichen) passt.
    """
    s = (
        f"rigid t=({tx:g},{ty:g},{tz:g}) "
        f"r=({rx:g},{ry:g},{rz:g}) "
        f"c={center_label} m={method[:3]} FoR={for_strategy}"
    )
    return _truncate(s, max_len)


def transform_rtstruct(
    rs_ds: pydicom.Dataset,
    T: np.ndarray,
    sop_map: dict,
    new_ct_series_uid: str,
    new_for_uid: "str | None" = None,
    drehpunkt_position: "np.ndarray | None" = None,
    label_suffix: "str | None" = None,
    description: "str | None" = None,
    series_number_offset: int = 1000,
    *,
    rs_sop_uid: "str | None" = None,
    rs_series_uid: "str | None" = None,
) -> pydicom.Dataset:
    """
    Wendet die rigide Transformation T auf das RTSTRUCT an.

    - Jeder ``ContourData``-Punkt wird durch T abgebildet (rein linear pro Punkt
      -> keine Verzerrung der Punktwolke; Volumina / Distanzen bleiben erhalten).
    - Per-Kontur-Verweise auf CT-Slices (``ContourImageSequence``) werden via
      ``sop_map`` auf die neuen SOP-UIDs umgesetzt; ebenso die Top-Level-
      Referenz auf die CT-Serie.
    - Wenn ``new_for_uid`` angegeben ist, werden alle FrameOfReferenceUID-
      Eintraege auf diesen Wert gesetzt; sonst bleibt die alte FoR.
    - Wenn ``drehpunkt_position`` angegeben ist, wird eine zusaetzliche POINT-
      ROI ``Drehpunkt`` an dieser Position eingefuegt (im transformierten
      Koordinatensystem).
    - Frische ``SOPInstanceUID`` und ``SeriesInstanceUID`` fuer das RS selbst
      (``rs_sop_uid`` / ``rs_series_uid``, wenn vorab vergeben).
    - ``InstanceCreationDate/Time`` wird gesetzt.

    Gibt das modifizierte RS-Dataset zurueck (Original wird nicht veraendert).
    """
    new_ds = copy.deepcopy(rs_ds)

    # 1. Konturpunkte transformieren + per-Kontur-CT-Refs umschreiben
    if hasattr(new_ds, "ROIContourSequence"):
        for roi_contour in new_ds.ROIContourSequence:
            if not hasattr(roi_contour, "ContourSequence"):
                continue
            for contour in roi_contour.ContourSequence:
                if not hasattr(contour, "ContourData"):
                    continue
                contour.ContourData = _apply_T_to_flat_coords(contour.ContourData, T)
                if hasattr(contour, "ContourImageSequence"):
                    _rewrite_referenced_sops(contour.ContourImageSequence, sop_map)

    # 2. Top-Level ReferencedFrameOfReferenceSequence -> neue CT-Serie + SOPs
    if hasattr(new_ds, "ReferencedFrameOfReferenceSequence"):
        for ref in new_ds.ReferencedFrameOfReferenceSequence:
            if new_for_uid is not None and hasattr(ref, "FrameOfReferenceUID"):
                ref.FrameOfReferenceUID = new_for_uid
            if not hasattr(ref, "RTReferencedStudySequence"):
                continue
            for study in ref.RTReferencedStudySequence:
                if not hasattr(study, "RTReferencedSeriesSequence"):
                    continue
                for series in study.RTReferencedSeriesSequence:
                    series.SeriesInstanceUID = new_ct_series_uid
                    if hasattr(series, "ContourImageSequence"):
                        _rewrite_referenced_sops(series.ContourImageSequence, sop_map)

    # 3. RTSTRUCT-Top-Level FoR + ROI-FoR-Refs
    if new_for_uid is not None:
        if hasattr(new_ds, "FrameOfReferenceUID"):
            new_ds.FrameOfReferenceUID = new_for_uid
        if hasattr(new_ds, "StructureSetROISequence"):
            for roi in new_ds.StructureSetROISequence:
                if hasattr(roi, "ReferencedFrameOfReferenceUID"):
                    roi.ReferencedFrameOfReferenceUID = new_for_uid

    # 4. Neue UIDs fuer das RS selbst (SOPInstanceUID auch im File-Meta-Header)
    mod.set_sop_instance_uid(new_ds, rs_sop_uid or generate_uid())
    new_ds.SeriesInstanceUID = rs_series_uid or generate_uid()

    now = datetime.now()
    new_ds.InstanceCreationDate = now.strftime("%Y%m%d")
    new_ds.InstanceCreationTime = now.strftime("%H%M%S.%f")[:13]

    # 4b. Aria-/TPS-sichtbare Metadaten
    if label_suffix:
        orig_label = str(getattr(new_ds, "StructureSetLabel", ""))
        new_ds.StructureSetLabel = _label_with_suffix(orig_label, label_suffix, 16)
        orig_name = str(getattr(new_ds, "StructureSetName", ""))
        if orig_name:
            new_ds.StructureSetName = _truncate(orig_name + label_suffix, 64)
        orig_series_desc = str(getattr(new_ds, "SeriesDescription", ""))
        new_ds.SeriesDescription = _truncate(orig_series_desc + label_suffix, 64)

    if description is not None:
        new_ds.StructureSetDescription = _truncate(description, 64)

    if series_number_offset:
        try:
            sn = int(getattr(new_ds, "SeriesNumber", 0) or 0)
            new_ds.SeriesNumber = sn + series_number_offset
        except (TypeError, ValueError):
            new_ds.SeriesNumber = series_number_offset

    # 5. Drehpunkt-Marker einfuegen
    if drehpunkt_position is not None:
        marker_for = new_for_uid if new_for_uid is not None \
            else str(getattr(new_ds, "FrameOfReferenceUID", ""))
        add_drehpunkt_marker(new_ds, drehpunkt_position, for_uid=marker_for)

    return new_ds


def align_contour_images(orig_ds: pydicom.Dataset, new_ds: pydicom.Dataset, ct_headers: list,
                         sop_map: dict, tol: float = PLANE_TOL_MM) -> dict:
    """
    Nur ``resample``: das CT-Raster bleibt, die Konturen wandern.  Jede Kontur
    verweist danach (erster Eintrag der ``ContourImageSequence``) auf die
    Schicht an ihrer neuen Lage, die naechste Schichtebene entlang der
    Schichtnormalen, statt auf die Schicht mit dem alten Index.

    Zaehlt die Konturen (ohne POINT), die vorher in einer Schichtebene lagen
    und jetzt nicht mehr: gekippt (Rotation um X/Y) oder zwischen zwei
    Schichten (z-Verschiebung kein Vielfaches des Schichtabstands).  Konturen
    jenseits der ersten oder letzten Schicht zaehlen nicht; die meldet das
    Clipping.
    """
    iop = np.asarray(ct_headers[0].ImageOrientationPatient, dtype=np.float64)
    normal = np.cross(iop[:3], iop[3:])
    normal /= np.linalg.norm(normal)
    planes = np.array([float(normal @ np.asarray(h.ImagePositionPatient, dtype=np.float64))
                       for h in ct_headers])
    sops = [str(h.SOPInstanceUID) for h in ct_headers]

    def nearest(pos: np.ndarray) -> tuple:
        k = int(np.argmin(np.abs(planes - pos.mean())))
        return k, float(np.max(np.abs(pos - planes[k])))

    stats = {"n_contours": 0, "n_off_plane": 0, "n_tilted": 0, "max_offset_mm": 0.0}
    # das neue RS hat am Ende zusaetzlich den Drehpunkt (ohne Vorgaenger)
    for rc_old, rc_new in zip(orig_ds.get("ROIContourSequence", []), new_ds.get("ROIContourSequence", []),
                              strict=False):
        for c_old, c_new in zip(rc_old.get("ContourSequence", []), rc_new.get("ContourSequence", []),
                                strict=True):
            if "ContourData" not in c_new:
                continue
            pos = np.asarray(c_new.ContourData, dtype=np.float64).reshape(-1, 3) @ normal
            k, offset = nearest(pos)
            refs = c_new.get("ContourImageSequence")
            if refs:
                refs[0].ReferencedSOPInstanceUID = sop_map[sops[k]]
            if str(c_new.get("ContourGeometricType", "")).upper() == "POINT":
                continue
            stats["n_contours"] += 1
            before = np.asarray(c_old.ContourData, dtype=np.float64).reshape(-1, 3) @ normal
            inside = planes.min() - tol <= pos.mean() <= planes.max() + tol
            if offset <= tol or not inside or nearest(before)[1] > tol:
                continue
            stats["n_off_plane"] += 1
            stats["n_tilted"] += int(float(pos.max() - pos.min()) > tol)
            stats["max_offset_mm"] = max(stats["max_offset_mm"], offset)
    return stats


# ---------------------------------------------------------------------------
# Marker / Rotationszentrum (Stage 3)
# ---------------------------------------------------------------------------

def print_marker_table(markers: list[tuple[str, np.ndarray]]) -> None:
    """Druckt eine kompakte Tabelle aller POINT-Marker."""
    print("\nVerfuegbare POINT-Marker im RTSTRUCT:")
    if not markers:
        print("  (keine)")
        return
    print(f"  {'Idx':<4} {'Name':<28} {'X (L)':>10} {'Y (P)':>10} {'Z (S)':>10}")
    for i, (name, pos) in enumerate(markers, 1):
        print(f"  [{i:<2}] {name:<28} {pos[0]:>10.2f} {pos[1]:>10.2f} {pos[2]:>10.2f}")


def parse_center_spec(
    spec: str,
    rs_ds: pydicom.Dataset,
    volume_center: np.ndarray,
) -> np.ndarray:
    """
    Loest einen Center-Spec-String in eine 3D-Position auf.

    Erlaubte Formen:
      - ``"volume"``               -> Volumenzentrum
      - ``"marker:NAME"``          -> POINT-Marker mit Namen NAME (case-insensitive)
      - ``"x,y,z"``                -> drei kommagetrennte Floats (LPS, mm)
    """
    raw = spec.strip()
    if not raw:
        raise ValueError("Leerer --center Wert")

    if raw.lower() == "volume":
        return volume_center.copy()

    if raw.lower().startswith("marker:"):
        target = raw[len("marker:"):].strip()
        markers = find_point_markers(rs_ds)
        for name, pos in markers:
            if name.lower() == target.lower():
                return pos.copy()
        avail = ", ".join(n for n, _ in markers) or "(keine)"
        raise ValueError(
            f"Marker '{target}' nicht im RTSTRUCT gefunden. Verfuegbar: {avail}"
        )

    parts = [p.strip() for p in raw.split(",")]
    if len(parts) != 3:
        raise ValueError(
            f"Ungueltiger --center Wert {raw!r}. "
            "Erlaubt: 'volume', 'marker:NAME', oder 'x,y,z'."
        )
    try:
        return np.array([float(p) for p in parts], dtype=np.float64)
    except ValueError as e:
        raise ValueError(f"Ungueltige Koordinaten in --center {raw!r}: {e}")


def interactive_center_prompt(
    rs_ds: pydicom.Dataset,
    volume_center: np.ndarray,
) -> np.ndarray:
    """
    Interaktiver Prompt: Marker-Liste anzeigen, Auswahl per Index entgegennehmen,
    plus Optionen ``v`` (Volumenzentrum, Default) und ``m`` (manuelle Eingabe).
    """
    markers = find_point_markers(rs_ds)
    print()
    print("Bitte Rotationszentrum waehlen:")
    for i, (name, pos) in enumerate(markers, 1):
        print(f"  [{i:<2}] Marker {name:<24} ({pos[0]:7.2f}, {pos[1]:7.2f}, {pos[2]:7.2f}) mm")
    print(f"  [v ] Volumenzentrum                 "
          f"({volume_center[0]:7.2f}, {volume_center[1]:7.2f}, {volume_center[2]:7.2f}) mm  [Default]")
    print(f"  [m ] Manueller Punkt (Eingabe x,y,z)")

    while True:
        try:
            raw = input("Auswahl [v]: ").strip()
        except EOFError:
            print("\n  Keine Eingabe -> Volumenzentrum.")
            return volume_center.copy()

        if raw == "" or raw.lower() == "v":
            return volume_center.copy()

        if raw.lower() == "m":
            coords = input("  Koordinaten x,y,z [mm, LPS]: ").strip()
            try:
                vals = [float(p) for p in coords.split(",")]
                if len(vals) != 3:
                    raise ValueError("Genau 3 Werte erwartet (x,y,z)")
                return np.array(vals, dtype=np.float64)
            except ValueError as e:
                print(f"  Ungueltig: {e}.  Bitte erneut.")
                continue

        if raw.isdigit() and markers:
            idx = int(raw)
            if 1 <= idx <= len(markers):
                return markers[idx - 1][1].copy()
            print(f"  Index ausserhalb [1..{len(markers)}].  Bitte erneut.")
            continue

        print("  Bitte eine Zahl 1..N, 'v' oder 'm' eingeben.")


# ---------------------------------------------------------------------------
# Hauptablauf
# ---------------------------------------------------------------------------

@dataclass
class CasePreflight:
    """Ergebnis von ``preflight_case``: geprueft, nur Header gelesen, nichts geschrieben."""
    case_dir: Path
    case_id: str
    ct_dir: Path
    rs_path: Path
    siblings: list                 # RP*/RD*-Dateien daneben (werden nicht transformiert)
    ct_headers: list               # CT-Header ohne Pixel, nach z sortiert
    rs_ds: pydicom.Dataset
    ct_for_uid: str
    geom: dict                     # modifier.extract_geometry der Header
    volume_center: np.ndarray
    markers: list                  # find_point_markers(rs_ds)
    issues: list = field(default_factory=list)


def _printer(quiet: bool):
    return (lambda *a, **k: None) if quiet else print


def preflight_case(case_dir: "str | None", rs_override: "str | None" = None, label: str = "_RB",
                   *, ct_files: "list | None" = None, case_id: "str | None" = None,
                   siblings: "list | None" = None, quiet: bool = False) -> CasePreflight:
    """
    Stufe 1: Discovery, CT-Header (ohne Pixel), Geometrie- und
    Orientierungspruefung, RTSTRUCT, FoR-Konsistenz, ``--label``.  Druckt den
    Kopf des Laufs (``quiet``: nichts).  Eingabefehler ->
    ``ValueError``/``FileNotFoundError``.

    ``ct_files`` (nur mit ``rs_override``) ersetzt die Ordner-Konvention
    ``<case>/CT/`` (z.B. flacher Export); ``case_id`` ersetzt dann den
    Ordnernamen, ``siblings`` die gefundenen RP-/RD-Dateien.
    """
    ctx = _runtime.current()
    ctx.stage("preflight", "Fall pruefen")
    say = _printer(quiet)
    if ct_files is None:
        ct_dir, rs_path, siblings = discover_case(case_dir, rs_override=rs_override,
                                                  return_siblings=True)
        ct_source = ct_dir
    else:
        if rs_override is None:
            raise ValueError("Mit einer CT-Dateiliste muss das RTSTRUCT angegeben werden (--rs).")
        rs_path = Path(rs_override)
        if not rs_path.is_file():
            raise FileNotFoundError(f"RS-Datei nicht gefunden: {rs_override!r}")
        ct_source = [str(f) for f in ct_files]
        parents = {Path(f).parent for f in ct_source}
        ct_dir = parents.pop() if len(parents) == 1 else rs_path.parent
        siblings = [Path(p) for p in (siblings or [])]
    base = Path(case_dir) if case_dir is not None else rs_path.parent
    say(f"\nLade Case '{case_dir if case_dir is not None else base}' …")
    issues = []
    if siblings:
        say(_siblings_note(siblings))
        issues.append(Issue(
            "info", "CASE.SIBLINGS_NOT_TRANSFORMED",
            f"{len(siblings)} Plan- und Dosisdatei(en) im Fallordner werden nicht mit-transformiert.",
            hint_de="RTPLAN/RTDOSE passen nach der Transformation nicht mehr zum CT.",
            detail="\n".join(Path(p).name for p in siblings)))
    say(f"  CT-Ordner   : {ct_dir}")
    say(f"  RTSTRUCT    : {rs_path}")

    headers = load_ct_headers(ct_source)
    say(f"  {len(headers)} CT-Slices geladen")
    validate_ct_geometry(headers)

    rs_ds = pydicom.dcmread(str(rs_path))
    if getattr(rs_ds, "Modality", None) != "RTSTRUCT":
        raise ValueError(
            f"Datei {rs_path!r} ist keine RTSTRUCT (Modalitaet: "
            f"{getattr(rs_ds, 'Modality', '?')})."
        )
    ct_for_uid = get_ct_frame_of_reference(headers)
    validate_for_consistency(ct_for_uid, rs_ds)
    say(f"  FrameOfReferenceUID OK ({ct_for_uid[:24]}…)")

    case_id = case_id or base.resolve().name
    validate_label(label, case_id)
    geom = mod.extract_geometry(headers)
    vol_c = mod.volume_center(geom)
    nz, ny, nx = geom["shape"]
    say(f"  Volumengroesse: {nz} x {ny} x {nx}  Voxel")
    say(f"  Volumen-Mitte : ({vol_c[0]:.1f}, {vol_c[1]:.1f}, {vol_c[2]:.1f}) mm")
    return CasePreflight(
        case_dir=base, case_id=case_id, ct_dir=ct_dir, rs_path=rs_path,
        siblings=siblings, ct_headers=headers, rs_ds=rs_ds, ct_for_uid=ct_for_uid, geom=geom,
        volume_center=vol_c, markers=find_point_markers(rs_ds), issues=issues,
    )


def resolve_center(spec: "str | None", rs_ds: pydicom.Dataset, volume_center_lps: np.ndarray,
                   interactive: bool = True) -> "tuple[np.ndarray | None, str]":
    """
    Rotationszentrum aus ``spec`` ('volume', 'marker:NAME', 'x,y,z').  Ohne
    ``spec`` fragt ein Prompt nach, aber nur wenn ``interactive`` und stdin ein
    Terminal ist (ohne Konsole, z.B. pythonw oder GUI-EXE, ist stdin None);
    sonst Volumenmitte.  Rueckgabe: (Position oder None fuer die Volumenmitte,
    Label fuer Ausgabe und Beschreibung).
    """
    if spec is not None:
        pos = parse_center_spec(spec, rs_ds, volume_center_lps)
        key = spec.lower().strip()
        if key == "volume":
            return None, "Volumenmitte"
        if key.startswith("marker:"):
            return pos, f"Marker '{spec.split(':', 1)[1].strip()}'"
        return pos, "manuell"
    if not interactive or sys.stdin is None or not sys.stdin.isatty():
        return None, "Volumenmitte"
    pos = interactive_center_prompt(rs_ds, volume_center_lps)
    if np.allclose(pos, volume_center_lps, atol=1e-9):
        return None, "Volumenmitte"
    return pos, "interaktiv"


@dataclass
class TransformPlan:
    """
    Ergebnis von ``plan_transform``: Matrix, Pfade, vorab vergebene UIDs, das
    fertig transformierte RTSTRUCT und der Clipping-Befund.  Noch ist nichts
    geschrieben; Fehler im RTSTRUCT (z.B. Verweise auf fremde CT-Schichten)
    sind hier bereits aufgefallen.
    """
    params: dict                   # tx, ty, tz, rx, ry, rz
    method: str
    order: int
    label: str
    T: np.ndarray
    center: np.ndarray
    center_label: str
    drehpunkt_pos: np.ndarray
    case_out: Path
    ct_out: Path
    rs_out: Path
    ct_series_uid: str
    sop_map: dict                  # alte CT-SOP -> neue CT-SOP
    new_for_uid: "str | None"
    for_strategy: str              # keep | new
    series_number_offset: int
    new_rs: pydicom.Dataset
    clipping: list                 # [(roi_name, n_outside, n_total, frac)]
    issues: list = field(default_factory=list)


def plan_transform(pre: CasePreflight, tx: float, ty: float, tz: float,
                   rx: float, ry: float, rz: float, *, output_dir: "str | None" = None,
                   method: str = "resample", order: int = 1, label: str = "_RB",
                   center: "np.ndarray | None" = None, center_label: str = "Volumenmitte",
                   new_frame_of_reference: bool = False,
                   series_number_offset: int = 1000, out_dir: "str | None" = None,
                   quiet: bool = False) -> TransformPlan:
    """
    Stufe 2: Transformationsmatrix, Ausgabepfade und alle neuen UIDs; das
    RTSTRUCT wird schon hier transformiert und auf Clipping geprueft, damit
    beim Schreiben kein halber Ausgabeordner entstehen kann.  Druckt den
    Transformationsblock (``quiet``: nichts).  Ergebnisordner ist ``out_dir``
    oder ``<output_dir>/<case_id><label>``.
    """
    ctx = _runtime.current()
    ctx.check_cancel()
    ctx.stage("plan", "Transformation planen")
    if method not in ("resample", "metadata"):
        raise ValueError(f"Unbekannte Methode: {method!r}")
    vol_c = pre.volume_center
    if center is None:
        center = vol_c
        resolved_label = center_label if center_label != "Volumenmitte" else "Volumenmitte"
    else:
        center = np.asarray(center, dtype=np.float64).reshape(3)
        resolved_label = center_label

    say = _printer(quiet)
    say(f"\nTransformation:")
    say(f"  Translation : tx={tx} mm, ty={ty} mm, tz={tz} mm")
    say(f"  Rotation    : rx={rx} deg, ry={ry} deg, rz={rz} deg  [intrinsisch XYZ]")
    say(f"  Methode     : {method}")
    say(f"  Zentrum     : {resolved_label}  "
        f"({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f}) mm")
    T = mod.build_rigid_transform(rx, ry, rz, tx, ty, tz, center)

    if out_dir is None and output_dir is None:
        raise ValueError("Ausgabeordner fehlt (output_dir oder out_dir).")
    case_out = Path(out_dir) if out_dir is not None else Path(output_dir) / f"{pre.case_id}{label}"
    new_for_uid = str(generate_uid()) if new_frame_of_reference else None
    for_strategy = "new" if new_frame_of_reference else "keep"
    ct_series_uid = str(generate_uid())
    sop_map = {str(s.SOPInstanceUID): str(generate_uid())
               for s in pre.ct_headers if getattr(s, "SOPInstanceUID", None)}
    # Drehpunkt im transformierten System: die Rotation laesst das Zentrum
    # invariant, also T(centre) = centre + (tx, ty, tz).
    drehpunkt_pos = center + np.array([tx, ty, tz])
    description = build_transform_description(
        tx, ty, tz, rx, ry, rz,
        center_label=resolved_label, method=method, for_strategy=for_strategy,
    )
    new_rs = transform_rtstruct(
        rs_ds=pre.rs_ds, T=T, sop_map=sop_map, new_ct_series_uid=ct_series_uid,
        new_for_uid=new_for_uid, drehpunkt_position=drehpunkt_pos, label_suffix=label,
        description=description, series_number_offset=series_number_offset,
    )
    # resample: das Raster bleibt, also zeigen die Bildverweise auf die Schicht an der neuen Lage
    align = align_contour_images(pre.rs_ds, new_rs, pre.ct_headers, sop_map) if method == "resample" else None
    clipping = check_contour_clipping(new_rs, pre.geom, method)

    issues = []
    if align and align["n_off_plane"]:
        n_off, n_all, n_tilt = align["n_off_plane"], align["n_contours"], align["n_tilted"]
        hints = ["Ein TPS kann solche Konturen beim Import verwerfen oder auf die naechste Schicht legen."]
        if n_tilt:
            hints.append("Rotationen um die Links-Rechts- und die anterior-posteriore Achse bildet nur "
                         "die Methode 'metadata' exakt ab (schraege Schichten).")
        if n_off > n_tilt:
            hints.append(f"Verschiebung Z als Vielfaches des Schichtabstands ({pre.geom['dz']:g} mm) waehlen.")
        issues.append(Issue(
            "warning", "CASE.CONTOURS_OFF_PLANE",
            f"{n_off} von {n_all} Konturen liegen nach der Bewegung nicht mehr in einer CT-Schichtebene"
            + (f" ({n_tilt} gekippt)." if n_tilt else "."),
            hint_de=" ".join(hints), field="method" if n_tilt else "tz"))
        say(f"  ! Konturen  : {n_off} von {n_all} nicht mehr in einer CT-Schichtebene"
            + (f" ({n_tilt} gekippt)" if n_tilt else "")
            + f", bis {align['max_offset_mm']:.2f} mm daneben; ein TPS kann sie verwerfen.")
    if for_strategy == "keep":                     # gewollter Standard, daher nur Info
        issues.append(Issue(
            "info", "CASE.FOR_KEPT",
            "Original und transformierter Datensatz tragen dieselbe FrameOfReferenceUID.",
            hint_de="Das TPS legt vorhandene Plaene/Dosen des Originals ungeprueft auf das "
                    "transformierte CT; so laesst sich der Originalplan bei einem Lagerungsfehler "
                    "auswerten. Fuer getrennte Planung eine neue FoR vergeben.",
            field="new_frame_of_reference"))
    if clipping:
        issues.append(Issue(
            "warning", "CASE.CONTOUR_CLIPPING",
            f"{len(clipping)} Struktur(en) ragen nach der Transformation aus dem CT-Volumen.",
            hint_de="Dort zeigt das CT Luft statt Anatomie; die Methode 'metadata' vermeidet das.",
            field="method",
            detail="\n".join(f"{n}: {o} von {t} Punkten ({f:.1%})" for n, o, t, f in clipping)))
    return TransformPlan(
        params={"tx": tx, "ty": ty, "tz": tz, "rx": rx, "ry": ry, "rz": rz},
        method=method, order=order, label=label, T=T, center=center, center_label=resolved_label,
        drehpunkt_pos=drehpunkt_pos, case_out=case_out, ct_out=case_out / "CT",
        rs_out=case_out / f"RS{label}.dcm", ct_series_uid=ct_series_uid, sop_map=sop_map,
        new_for_uid=new_for_uid, for_strategy=for_strategy,
        series_number_offset=series_number_offset, new_rs=new_rs, clipping=clipping, issues=issues,
    )


def print_dry_run(pre: CasePreflight, plan: TransformPlan) -> dict:
    """``--dry-run``: Matrix und geplante Pfade ausgeben; nichts wird geschrieben."""
    print("\nDRY-RUN  ----  es werden KEINE Dateien geschrieben.")
    print("\nT-Matrix (Patient -> Patient):")
    print(np.array2string(plan.T, precision=4, suppress_small=True))
    print(f"\nGeplante Output-Pfade:")
    print(f"  CT-Verzeichnis : {plan.ct_out}")
    print(f"  RTSTRUCT       : {plan.rs_out}")
    return {
        "case_id":    pre.case_id,
        "dry_run":    True,
        "T":          plan.T.tolist(),
        "rotation_center": plan.center.tolist(),
        "rotation_center_label": plan.center_label,
        "planned_ct_output_dir": str(plan.ct_out),
        "planned_rs_output_path": str(plan.rs_out),
        "clipping": _clipping_dicts(plan.clipping),
        "issues": [i.to_dict() for i in pre.issues + plan.issues],
    }


def _clipping_dicts(clipping: list) -> list:
    return [{"roi": n, "n_outside": o, "n_total": t, "fraction": f} for n, o, t, f in clipping]


def execute_transform(pre: CasePreflight, plan: TransformPlan, *, verify: bool = False,
                      no_viz: bool = False, viz_ct_surface: bool = False) -> dict:
    """
    Stufe 3: die geprueften CT-Dateien mit Pixeln laden, transformieren und in
    einem Durchgang schreiben (SeriesNumber-Offset und vorab vergebene UIDs),
    danach das fertige RTSTRUCT; optional ``--verify`` und die
    Vorher/Nachher-Ansicht.
    """
    ctx = _runtime.current()
    ctx.check_cancel()
    ctx.stage("ct_transform", "CT transformieren")
    if plan.new_for_uid is not None:
        print(f"  FoR-Strategie: NEU ({plan.new_for_uid[:24]}…)")
    else:
        print(
            "  FoR-Strategie: KEEP (alte FoR wird beibehalten).\n"
            "  ! WARNUNG ! ----------------------------------------------------------\n"
            "    Beide Datensaetze (Original + transformiert) tragen DIESELBE\n"
            "    FrameOfReferenceUID, obwohl sie sich physikalisch unterscheiden.\n"
            "    Das TPS verlinkt sie automatisch ohne Geometrievalidierung.\n"
            "    Folgen:\n"
            "      - Vorhandene RTPLAN/RTDOSE des Originals werden auf das\n"
            "        transformierte CT geworfen, obwohl sie geometrisch nicht\n"
            "        mehr dazu passen.  Dosis-Overlays sind dann irrefuehrend.\n"
            "      - Mischen von Konturen aus Original-RS und transformiertem RS\n"
            "        in einem Plan ergibt klinisch falsche DVHs.\n"
            "    Verwende --new-frame-of-reference fuer eine saubere Trennung,\n"
            "    falls die Datensaetze unabhaengig voneinander geplant werden.\n"
            "  ----------------------------------------------------------------------"
        )

    # genau die geprueften Dateien, jetzt mit Pixeln
    slices = mod.load_ct_series_files([h.filename for h in pre.ct_headers])
    if [str(getattr(s, "SOPInstanceUID", "")) for s in slices] != \
            [str(getattr(h, "SOPInstanceUID", "")) for h in pre.ct_headers]:
        raise ValueError(f"CT-Ordner {str(pre.ct_dir)!r} hat sich seit der Pruefung geaendert; "
                         "bitte den Lauf neu starten.")
    need_hu = plan.method == "resample" or (not no_viz and viz_ct_surface)
    volume_hu = mod.slices_to_hu(slices) if need_hu else None
    plan.case_out.mkdir(parents=True, exist_ok=True)
    save_kw = dict(series_description_suffix=plan.label, frame_of_reference_uid=plan.new_for_uid,
                   series_uid=plan.ct_series_uid, sop_map=plan.sop_map,
                   series_number_offset=plan.series_number_offset)
    if plan.method == "metadata":
        print("\nAktualisiere DICOM-Metadaten (HU-Werte exakt erhalten) …")
        out_slices = mod.apply_metadata_transform(slices, plan.T)
        save_info = mod.save_ct_series(out_slices, str(plan.ct_out), **save_kw)
    else:
        print(f"\nNeuabtastung (Interpolationsordnung {plan.order}) …")
        new_volume = mod.resample_volume(volume_hu, pre.geom["affine"], plan.T, order=plan.order)
        save_info = mod.save_ct_series(slices, str(plan.ct_out), new_volume_hu=new_volume, **save_kw)

    ctx.check_cancel()
    ctx.stage("rs_write", "RTSTRUCT schreiben")
    print("\nTransformiere RTSTRUCT …")
    new_rs = plan.new_rs
    new_rs.save_as(str(plan.rs_out))
    dp = plan.drehpunkt_pos
    print(f"  RTSTRUCT geschrieben -> {plan.rs_out}")
    print(f"  Drehpunkt-Marker eingefuegt bei ({dp[0]:.2f}, {dp[1]:.2f}, {dp[2]:.2f}) mm")

    # Clipping (nur resample-Mode; in plan_transform berechnet)
    if plan.clipping:
        print("\n  ! CLIPPING-WARNUNG ! Konturpunkte ausserhalb des Output-CT-Grids:")
        print(f"  {'ROI':<28} {'aussen':>8} {'gesamt':>8} {'Anteil':>8}")
        for roi_name, n_out, n_total, frac in sorted(plan.clipping, key=lambda x: -x[3]):
            print(f"    {roi_name:<26} {n_out:>8} {n_total:>8} {frac:>7.1%}")
        print(
            "    Diese Strukturen ragen aus dem aufgenommenen Bildvolumen heraus.\n"
            "    Im TPS sind die Konturen sichtbar, aber das CT zeigt dort -1000 HU\n"
            "    (Luft) statt Anatomie.  DVH-Auswertungen werden 'kuenstlich besser'\n"
            "    aussehen, weil Volumenanteile schlicht fehlen.\n"
            "    Tipp: --method metadata vermeidet das (oblique Slices, exakte\n"
            "    Geometrie), wird aber von manchen aelteren TPS abgelehnt."
        )

    issues = list(pre.issues) + list(plan.issues)
    verify_report = None
    if verify:
        ctx.stage("verify", "Centroide pruefen")
        verify_report = _verify_rs_centroids(pre.rs_ds, str(plan.rs_out), plan.T)
        if not verify_report["passed"]:
            issues.append(Issue(
                "warning", "CASE.VERIFY_FAILED",
                f"Centroid-Pruefung: Abweichung {verify_report['max_err_mm']:.3e} mm ueber der "
                f"Schwelle {verify_report['threshold_mm']:g} mm.",
                hint_de="Transformierte Konturen nicht verwenden und den Fall melden."))

    # Vorher/Nachher-Visualisierung; Fehler duerfen den geschriebenen Transform
    # nicht entwerten -> defensiv abgefangen.
    viz_report = None
    if not no_viz:
        ctx.check_cancel()
        ctx.stage("viz", "Vorher/Nachher-Ansicht erstellen")
        try:
            from . import visualizer as viz
            viz_report = viz.run_case_visualization(
                orig_ds=pre.rs_ds,
                new_ds=new_rs,
                center=plan.center,
                drehpunkt_pos=plan.drehpunkt_pos,
                translation=(plan.params["tx"], plan.params["ty"], plan.params["tz"]),
                T=plan.T,
                output_dir=plan.case_out,
                markers=pre.markers,
                geom=pre.geom,
                volume_hu=volume_hu,
                ct_surface=viz_ct_surface,
            )
        except _runtime.JobCancelled:
            raise
        except Exception as e:
            print(f"\n  Hinweis: Visualisierung fehlgeschlagen ({e}). "
                  "Transform-Dateien sind dennoch gueltig geschrieben.")
            issues.append(Issue("warning", "CASE.VIZ_FAILED",
                                f"Visualisierung fehlgeschlagen ({type(e).__name__}: {e}).",
                                hint_de="Die Transform-Dateien sind trotzdem gueltig."))

    print("\nFertig.")
    return {
        "case_id":         pre.case_id,
        "output_dir":      str(plan.case_out),
        "ct_output_dir":   str(plan.ct_out),
        "rs_output_path":  str(plan.rs_out),
        "ct_series_uid":   save_info["series_uid"],
        "ct_for_uid":      save_info["frame_of_reference_uid_used"],
        "sop_map":         save_info["sop_map"],
        "rs_sop_uid":      str(new_rs.SOPInstanceUID),
        "rs_series_uid":   str(new_rs.SeriesInstanceUID),
        "rotation_center": plan.center.tolist(),
        "rotation_center_label": plan.center_label,
        "drehpunkt_pos":   plan.drehpunkt_pos.tolist(),
        "verify":          verify_report,
        "method":          plan.method,
        "for_strategy":    plan.for_strategy,
        "clipping":        _clipping_dicts(plan.clipping),
        "viz":             viz_report,
        "issues":          [i.to_dict() for i in issues],
    }


def run_case_transform(
    case_dir: str,
    output_dir: str,
    tx: float, ty: float, tz: float,
    rx: float, ry: float, rz: float,
    method: str = "resample",
    order: int = 1,
    label: str = "_RB",
    rs_override: "str | None" = None,
    center: "np.ndarray | None" = None,
    center_label: str = "Volumenmitte",
    new_frame_of_reference: bool = False,
    series_number_offset: int = 1000,
    dry_run: bool = False,
    verify: bool = False,
    no_viz: bool = False,
    viz_ct_surface: bool = False,
) -> dict:
    """
    Kompletter Lauf in drei Stufen: ``preflight_case`` -> ``plan_transform`` ->
    ``execute_transform`` (bzw. ``print_dry_run``).  Gibt ein Dict mit den
    neuen UIDs, Pfaden, Clipping-Befund und ``issues`` zurueck.
    """
    pre = preflight_case(case_dir, rs_override=rs_override, label=label)
    plan = plan_transform(pre, tx, ty, tz, rx, ry, rz, output_dir=output_dir, method=method,
                          order=order, label=label, center=center, center_label=center_label,
                          new_frame_of_reference=new_frame_of_reference,
                          series_number_offset=series_number_offset)
    if dry_run:
        return print_dry_run(pre, plan)
    return execute_transform(pre, plan, verify=verify, no_viz=no_viz, viz_ct_surface=viz_ct_surface)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Rigid Body Transform fuer CT + RTSTRUCT im Verbund. "
            "Stage 1: CT wird transformiert, RTSTRUCT wird unveraendert kopiert."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Beispiele:\n"
            "  python -m dicom_file_modifier.case_modifier data/<case-id>\n"
            "  python -m dicom_file_modifier.case_modifier data/<case-id> "
            "--tx 10 --rz 15\n"
        ),
    )
    p.add_argument("case_dir", help="Case-Ordner mit CT/-Unterordner und RS*.dcm")
    p.add_argument("--output", "-o", default="output",
                   metavar="DIR",
                   help="Basis-Ausgabeverzeichnis (Standard: output)")
    p.add_argument("--rs", dest="rs_override", default=None, metavar="PATH",
                   help="Explizite RTSTRUCT-Datei (falls mehrere RS*.dcm vorliegen)")
    p.add_argument("--label", default="_RB", metavar="TEXT",
                   help="Suffix fuer Output-Ordner und Dateinamen (Standard: _RB)")

    grp_t = p.add_argument_group("Translation [mm]")
    grp_t.add_argument("--tx", type=float, default=0.0, metavar="mm")
    grp_t.add_argument("--ty", type=float, default=0.0, metavar="mm")
    grp_t.add_argument("--tz", type=float, default=0.0, metavar="mm")

    grp_r = p.add_argument_group("Rotation [deg]  -  intrinsisch XYZ um Volumenmitte")
    grp_r.add_argument("--rx", type=float, default=0.0, metavar="deg")
    grp_r.add_argument("--ry", type=float, default=0.0, metavar="deg")
    grp_r.add_argument("--rz", type=float, default=0.0, metavar="deg")

    grp_m = p.add_argument_group("Methode und Qualitaet")
    grp_m.add_argument("--method", choices=["resample", "metadata"], default="resample")
    grp_m.add_argument("--order", type=int, choices=[0, 1, 3], default=1, metavar="N")

    grp_c = p.add_argument_group("Rotationszentrum")
    grp_c.add_argument("--center", default=None, metavar="SPEC",
                       help=("Rotationszentrum.  Erlaubte Formen: "
                             "'volume' (Default), 'marker:NAME', oder 'x,y,z' (LPS, mm).  "
                             "Ohne Angabe: interaktive Auswahl falls Marker vorhanden."))
    grp_c.add_argument("--list-markers", action="store_true",
                       help="POINT-Marker im RTSTRUCT auflisten und beenden.")
    grp_c.add_argument("--non-interactive", action="store_true",
                       help="Kein interaktiver Prompt; bei fehlendem --center "
                            "wird das Volumenzentrum genutzt.")

    grp_f = p.add_argument_group("FrameOfReferenceUID")
    grp_f.add_argument("--new-frame-of-reference", action="store_true",
                       help=("Neue FrameOfReferenceUID fuer transformiertes CT+RS "
                             "vergeben (verhindert versehentliche Ueberlagerung mit "
                             "alten Plaenen/Dosen).  Default: alte FoR beibehalten."))

    grp_viz = p.add_argument_group("Visualisierung")
    grp_viz.add_argument("--no-viz", action="store_true",
                         help="Keine Vorher/Nachher-Plots erzeugen "
                              "(transform_3d.html / transform_overview.png / "
                              "displacement.png).")
    grp_viz.add_argument("--viz-ct-surface", action="store_true",
                         help="Im 3D-HTML zusaetzlich die CT-Koerperoberflaeche "
                              "extrahieren (Marching Cubes, rechenintensiv; "
                              "Default: aus).")

    grp_v = p.add_argument_group("Validierung")
    grp_v.add_argument("--dry-run", action="store_true",
                       help="Nur validieren und Plan ausgeben; keine Dateien schreiben.")
    grp_v.add_argument("--verify", action="store_true",
                       help="Nach dem Schreiben Centroid-Linearitaet pro ROI pruefen.")
    grp_v.add_argument("--self-test", action="store_true",
                       help=("Identitaets-Transform end-to-end laufen lassen und "
                             "asserten, dass alle ContourData-Werte sich um <1e-4 mm "
                             "vom Original unterscheiden.  Exit 0 = pass, 1 = fail."))
    return p


def _resolve_center_from_args(
    args: argparse.Namespace,
    rs_ds: pydicom.Dataset,
    volume_center_lps: np.ndarray,
) -> "tuple[np.ndarray | None, str]":
    """CLI-Huelle um ``resolve_center`` (``--center``, ``--non-interactive``)."""
    return resolve_center(args.center, rs_ds, volume_center_lps,
                          interactive=not args.non_interactive)


def _run_self_test(args: argparse.Namespace) -> int:
    """
    Mehrstufiger Self-Test, der typische Bug-Spots abdeckt:

      1. Identitaets-Transform   ->  ContourData byte-aequivalent zum Original.
      2. Z-Rotation 15 deg       ->  alle Pruef-Metriken muessen identisch zum
                                     Original sein (Z-Achse ist die Schicht-
                                     normale, daher trivialerweise erhalten).
      3. X-Rotation 5 deg        ->  Punktabstaende (nicht Volumina!) muessen
                                     bis auf Float-Round-Trip exakt erhalten
                                     bleiben - das ist die echte Rigid-Body-
                                     Eigenschaft.

    Ergebnis 0 = alle drei PASS, sonst 1.
    """
    import tempfile
    from .analyzer import load_rtstruct, extract_contours, get_structure_names

    _, rs_path, _ = discover_case(args.case_dir, rs_override=args.rs_override, return_siblings=True)
    orig = load_rtstruct(str(rs_path))
    names = get_structure_names(orig)

    print(f"\n--- Self-Test auf '{args.case_dir}' ---")

    def _max_point_diff(orig_ds, new_ds) -> tuple[float, str]:
        max_d = 0.0
        worst = None
        for roi_num, roi_name in names.items():
            co = extract_contours(orig_ds, roi_num)
            cn = extract_contours(new_ds, roi_num)
            if not co or not cn or len(co) != len(cn):
                continue
            for a, b in zip(co, cn):
                if a.shape != b.shape:
                    continue
                d = float(np.max(np.abs(a - b)))
                if d > max_d:
                    max_d = d
                    worst = roi_name
        return max_d, worst

    def _max_pairwise_distance_drift(orig_ds, new_ds, n_samples=200) -> float:
        """
        Rigid transforms preserve distances.  Wir samplen pro ROI bis zu n_samples
        Punktepaare und vergleichen ||p_i - p_j|| zwischen Original und neu.
        """
        rng = np.random.default_rng(0)
        max_drift = 0.0
        for roi_num in names:
            co = extract_contours(orig_ds, roi_num)
            cn = extract_contours(new_ds, roi_num)
            if not co or not cn:
                continue
            po = np.vstack(co)
            pn = np.vstack(cn)
            if po.shape != pn.shape or po.shape[0] < 2:
                continue
            n = po.shape[0]
            k = min(n_samples, n)
            idx = rng.choice(n, size=k, replace=False)
            do = np.linalg.norm(po[idx][:, None] - po[idx][None, :], axis=-1)
            dn = np.linalg.norm(pn[idx][:, None] - pn[idx][None, :], axis=-1)
            max_drift = max(max_drift, float(np.max(np.abs(do - dn))))
        return max_drift

    overall_pass = True
    with tempfile.TemporaryDirectory(prefix="case_modifier_selftest_") as tmp:
        # Test 1: Identity
        info = run_case_transform(
            case_dir=args.case_dir, output_dir=tmp,
            tx=0, ty=0, tz=0, rx=0, ry=0, rz=0,
            method="metadata", label=args.label,
            rs_override=args.rs_override,
            no_viz=True,
        )
        new = load_rtstruct(info["rs_output_path"])
        d1, w1 = _max_point_diff(orig, new)
        ok1 = d1 < 1e-4
        print(f"  [1/3] Identity round-trip:    "
              f"max point dev {d1:.3e} mm  "
              f"({'PASS' if ok1 else 'FAIL'}, worst: {w1})")
        overall_pass &= ok1

    with tempfile.TemporaryDirectory(prefix="case_modifier_selftest_") as tmp:
        # Test 2: Z-Rotation
        info = run_case_transform(
            case_dir=args.case_dir, output_dir=tmp,
            tx=0, ty=0, tz=0, rx=0, ry=0, rz=15,
            method="metadata", label=args.label,
            rs_override=args.rs_override,
            no_viz=True,
        )
        new = load_rtstruct(info["rs_output_path"])
        drift2 = _max_pairwise_distance_drift(orig, new)
        ok2 = drift2 < 1e-3
        print(f"  [2/3] Z-rotation 15 deg:      "
              f"max distance drift {drift2:.3e} mm  "
              f"({'PASS' if ok2 else 'FAIL'})")
        overall_pass &= ok2

    with tempfile.TemporaryDirectory(prefix="case_modifier_selftest_") as tmp:
        # Test 3: X-Rotation (echter Rigid-Body-Test)
        info = run_case_transform(
            case_dir=args.case_dir, output_dir=tmp,
            tx=0, ty=0, tz=0, rx=5, ry=0, rz=0,
            method="metadata", label=args.label,
            rs_override=args.rs_override,
            no_viz=True,
        )
        new = load_rtstruct(info["rs_output_path"])
        drift3 = _max_pairwise_distance_drift(orig, new)
        ok3 = drift3 < 1e-3
        print(f"  [3/3] X-rotation 5 deg:       "
              f"max distance drift {drift3:.3e} mm  "
              f"({'PASS' if ok3 else 'FAIL'})")
        overall_pass &= ok3

    print(f"\n  Gesamt: {'PASS' if overall_pass else 'FAIL'}")
    print(
        "\n  Hinweis: Punkte-zu-Punkte-Distanzen sind die definitive Rigid-Body-\n"
        "  Eigenschaft (immer erhalten).  Die Volumina aus analyzer.py sind das\n"
        "  NICHT - die planar-axiale Shoelace-Formel gibt bei Nicht-Z-Rotationen\n"
        "  Abweichungen im Prozent-Bereich, weil Konturen schraeg zur Z-Achse\n"
        "  liegen.  Das ist ein Mess-Artefakt der Volumenformel, kein Bug der\n"
        "  Transformation."
    )
    return 0 if overall_pass else 1


def main(argv: "list[str] | None" = None) -> int:
    args = _build_parser().parse_args(argv)
    try:
        # --self-test: dedizierter Code-Pfad, kein User-Output
        if args.self_test:
            return _run_self_test(args)

        # --list-markers: nur RS oeffnen und Marker auflisten, dann beenden.
        if args.list_markers:
            _, rs_path, _ = discover_case(args.case_dir, rs_override=args.rs_override,
                                          return_siblings=True)
            rs_ds = pydicom.dcmread(str(rs_path))
            print_marker_table(find_point_markers(rs_ds))
            return 0

        # Stufen wie run_case_transform; das Rotationszentrum wird nach der
        # (Header-)Pruefung gewaehlt, damit der Prompt Fall und Volumenmitte zeigt.
        pre = preflight_case(args.case_dir, rs_override=args.rs_override, label=args.label)
        center, center_label = _resolve_center_from_args(args, pre.rs_ds, pre.volume_center)
        plan = plan_transform(
            pre, args.tx, args.ty, args.tz, args.rx, args.ry, args.rz,
            output_dir=args.output, method=args.method, order=args.order, label=args.label,
            center=center, center_label=center_label,
            new_frame_of_reference=args.new_frame_of_reference,
        )
        if args.dry_run:
            print_dry_run(pre, plan)
        else:
            execute_transform(pre, plan, verify=args.verify, no_viz=args.no_viz,
                              viz_ct_surface=args.viz_ct_surface)
    except (FileNotFoundError, ValueError, KeyError) as e:
        print(f"\nFehler: {e}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
