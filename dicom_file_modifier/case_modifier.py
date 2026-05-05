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
+ extrinsische XYZ-Euler-Rotation) angewendet, die ``modifier.py`` bereits
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
  python -m dicom_file_modifier.case_modifier data/0000000171

  # Reale Transformation (Stage 1 transformiert nur das CT)
  python -m dicom_file_modifier.case_modifier data/0000000171 \
      --tx 10 --ty 0 --tz -5 --rx 0 --ry 0 --rz 15 --output output/run1
"""

from __future__ import annotations

import argparse
import copy
import os
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pydicom
from pydicom.dataset import Dataset
from pydicom.sequence import Sequence
from pydicom.uid import generate_uid

from . import modifier as mod


# ---------------------------------------------------------------------------
# Auto-Discovery
# ---------------------------------------------------------------------------

def discover_case(case_dir: str, rs_override: "str | None" = None) -> tuple[Path, Path]:
    """
    Sucht im ``case_dir`` den ``CT/``-Unterordner und genau eine ``RS*.dcm``-Datei.

    Wirft ``FileNotFoundError`` / ``ValueError`` mit klarer Fehlermeldung,
    falls die Konvention verletzt ist.  Bei mehreren RS-Dateien wird der Aufrufer
    aufgefordert, mit ``--rs <pfad>`` explizit auszuwaehlen.
    Hinweis-Druck (kein Error), wenn parallel RP*/RD*-Dateien vorhanden sind.
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
    if extras:
        names = ", ".join(p.name for p in extras)
        print(
            f"  Hinweis: Zusaetzliche Plan-/Dosis-Dateien gefunden "
            f"({names}). Diese werden NICHT mit-transformiert."
        )

    return ct_dir, rs_path


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


def get_rs_frame_of_references(rs_ds: pydicom.Dataset) -> set[str]:
    """Sammelt alle im RTSTRUCT referenzierten FrameOfReferenceUIDs."""
    uids: set[str] = set()
    if hasattr(rs_ds, "FrameOfReferenceUID"):
        uids.add(str(rs_ds.FrameOfReferenceUID))
    if hasattr(rs_ds, "ReferencedFrameOfReferenceSequence"):
        for ref in rs_ds.ReferencedFrameOfReferenceSequence:
            if hasattr(ref, "FrameOfReferenceUID"):
                uids.add(str(ref.FrameOfReferenceUID))
    return uids


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

def _verify_rs_centroids(orig_ds: pydicom.Dataset,
                         new_rs_path: str,
                         T: np.ndarray) -> dict:
    """
    Liest das frisch geschriebene RTSTRUCT wieder ein, berechnet pro ROI den
    Centroid und vergleicht ihn mit ``T @ centroid_orig``.  Druckt die maximale
    Abweichung pro ROI sowie die Gesamtstatistik.

    Centroids sind unter rigiden Transformationen linear: T(mean(p_i)) =
    mean(T(p_i)).  Eine signifikante Abweichung deutet auf einen Indizierungs-
    oder Reshape-Bug im Transform-Pfad hin.

    Rueckgabe: dict mit ``"max_err_mm"``, ``"checked"``, ``"worst_roi"``.
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
    return {"max_err_mm": max_err, "checked": checked, "worst_roi": worst}


def validate_ct_geometry(slices: list, atol_iop: float = 1e-3,
                        rel_tol_spacing: float = 0.01) -> None:
    """
    Pre-Flight-Check der CT-Geometrie.  Wirft ``ValueError`` wenn:
      - weniger als 2 Slices
      - ImageOrientationPatient variiert ueber Slices (>= ``atol_iop``)
      - PixelSpacing nicht einheitlich
      - Slice-Spacing nicht uniform (relative Abweichung >= ``rel_tol_spacing``)
    """
    if len(slices) < 2:
        raise ValueError(
            f"Mindestens 2 CT-Slices benoetigt fuer Geometrie-Berechnung "
            f"(gefunden: {len(slices)})."
        )

    iop_ref = np.array([float(x) for x in slices[0].ImageOrientationPatient])
    ps_ref  = [float(x) for x in slices[0].PixelSpacing]
    for i, s in enumerate(slices[1:], 1):
        iop = np.array([float(x) for x in s.ImageOrientationPatient])
        if not np.allclose(iop, iop_ref, atol=atol_iop):
            raise ValueError(
                f"ImageOrientationPatient variiert ueber Slices.\n"
                f"  Slice 0: {iop_ref.tolist()}\n  Slice {i}: {iop.tolist()}"
            )
        ps = [float(x) for x in s.PixelSpacing]
        if not np.allclose(ps, ps_ref, atol=1e-4):
            raise ValueError(
                f"PixelSpacing variiert ueber Slices.\n"
                f"  Slice 0: {ps_ref}\n  Slice {i}: {ps}"
            )

    # Slice-Spacing (Distanzen zwischen aufeinanderfolgenden IPPs)
    ipps = np.array([[float(x) for x in s.ImagePositionPatient] for s in slices])
    diffs = np.linalg.norm(np.diff(ipps, axis=0), axis=1)
    mean = float(np.mean(diffs))
    if mean <= 0:
        raise ValueError("Slice-Spacing-Mittelwert ist 0  Slices ueberlagern sich?")
    rel_dev = float(np.max(np.abs(diffs - mean) / mean))
    if rel_dev > rel_tol_spacing:
        raise ValueError(
            f"Slice-Spacing nicht uniform (relative Abweichung {rel_dev:.2%} "
            f"> {rel_tol_spacing:.0%}).\n"
            f"  min={diffs.min():.3f} mm, max={diffs.max():.3f} mm, "
            f"mean={mean:.3f} mm"
        )


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


def _truncate(value: str, max_len: int) -> str:
    """Schneidet einen String auf die DICOM-VR-Laenge ohne Encoding-Tricks."""
    return value[:max_len]


def _label_with_suffix(orig: str, suffix: str, max_len: int) -> str:
    """
    Haengt ``suffix`` an ``orig`` an und kuerzt das Ergebnis sauber auf
    ``max_len``.  Wenn ``orig + suffix`` zu lang ist, wird ``orig`` so weit
    gekuerzt, dass das Suffix vollstaendig erhalten bleibt.
    """
    if len(orig) + len(suffix) <= max_len:
        return orig + suffix
    keep = max(0, max_len - len(suffix))
    return (orig[:keep] + suffix)[:max_len]


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
    - Frische ``SOPInstanceUID`` und ``SeriesInstanceUID`` fuer das RS selbst.
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

    # 4. Neue UIDs fuer das RS selbst
    new_ds.SOPInstanceUID    = generate_uid()
    new_ds.SeriesInstanceUID = generate_uid()

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


# ---------------------------------------------------------------------------
# Marker / Rotationszentrum (Stage 3)
# ---------------------------------------------------------------------------

def find_point_markers(rs_ds: pydicom.Dataset) -> list[tuple[str, np.ndarray]]:
    """
    Liefert ``[(roi_name, position_lps), ...]`` fuer alle ROIs, deren
    ``ROIContourSequence`` mindestens eine POINT-Type-Kontur enthaelt.
    Reihenfolge entspricht der ``StructureSetROISequence``.
    """
    if not hasattr(rs_ds, "StructureSetROISequence"):
        return []
    name_map = {int(r.ROINumber): str(r.ROIName) for r in rs_ds.StructureSetROISequence}

    markers: list[tuple[str, np.ndarray]] = []
    if not hasattr(rs_ds, "ROIContourSequence"):
        return markers
    for rc in rs_ds.ROIContourSequence:
        if not hasattr(rc, "ContourSequence"):
            continue
        for c in rc.ContourSequence:
            if str(getattr(c, "ContourGeometricType", "")) != "POINT":
                continue
            pts = np.array(c.ContourData, dtype=np.float64).reshape(-1, 3)
            if pts.shape[0] == 0:
                continue
            roi_num = int(getattr(rc, "ReferencedROINumber", -1))
            markers.append((name_map.get(roi_num, f"ROI#{roi_num}"), pts[0].copy()))
            break  # erster POINT je ROI reicht
    return markers


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
) -> dict:
    """
    Stage-1-Implementierung: CT transformieren, RS unveraendert kopieren.

    Gibt ein Dict mit den wichtigsten neuen UIDs / Pfaden zurueck (fuer Tests
    und nachgelagerte Stages nutzbar).
    """
    # ── 1. Discovery + Pre-Flight ────────────────────────────────────────────
    print(f"\nLade Case '{case_dir}' …")
    ct_dir, rs_path = discover_case(case_dir, rs_override=rs_override)
    print(f"  CT-Ordner   : {ct_dir}")
    print(f"  RTSTRUCT    : {rs_path}")

    slices = mod.load_ct_series(str(ct_dir))
    print(f"  {len(slices)} CT-Slices geladen")
    validate_ct_geometry(slices)

    rs_ds = pydicom.dcmread(str(rs_path))
    if getattr(rs_ds, "Modality", None) != "RTSTRUCT":
        raise ValueError(
            f"Datei {rs_path!r} ist keine RTSTRUCT (Modalitaet: "
            f"{getattr(rs_ds, 'Modality', '?')})."
        )

    ct_for_uid = get_ct_frame_of_reference(slices)
    validate_for_consistency(ct_for_uid, rs_ds)
    print(f"  FrameOfReferenceUID OK ({ct_for_uid[:24]}…)")

    # ── 2. Output-Layout vorbereiten ─────────────────────────────────────────
    case_id = Path(case_dir).resolve().name
    case_out = Path(output_dir) / f"{case_id}{label}"
    ct_out   = case_out / "CT"
    rs_out   = case_out / f"RS{label}.dcm"
    if not dry_run:
        case_out.mkdir(parents=True, exist_ok=True)

    # ── 3. CT-Geometrie & Transform-Matrix ───────────────────────────────────
    volume_hu = mod.slices_to_hu(slices)
    geom      = mod.extract_geometry(slices)
    vol_c     = mod.volume_center(geom)
    if center is None:
        center = vol_c
        resolved_label = center_label if center_label != "Volumenmitte" else "Volumenmitte"
    else:
        center = np.asarray(center, dtype=np.float64).reshape(3)
        resolved_label = center_label

    nz, ny, nx = geom["shape"]
    print(f"  Volumengroesse: {nz} x {ny} x {nx}  Voxel")
    print(f"  Volumen-Mitte : ({vol_c[0]:.1f}, {vol_c[1]:.1f}, {vol_c[2]:.1f}) mm")

    print(f"\nTransformation:")
    print(f"  Translation : tx={tx} mm, ty={ty} mm, tz={tz} mm")
    print(f"  Rotation    : rx={rx} deg, ry={ry} deg, rz={rz} deg  [extrinsisch XYZ]")
    print(f"  Methode     : {method}")
    print(f"  Zentrum     : {resolved_label}  "
          f"({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f}) mm")

    T = mod.build_rigid_transform(rx, ry, rz, tx, ty, tz, center)

    # ── Dry-Run: nur Plan ausgeben, nichts schreiben ─────────────────────────
    if dry_run:
        print("\nDRY-RUN  ----  es werden KEINE Dateien geschrieben.")
        print("\nT-Matrix (Patient -> Patient):")
        print(np.array2string(T, precision=4, suppress_small=True))
        print(f"\nGeplante Output-Pfade:")
        case_id_dr = Path(case_dir).resolve().name
        case_out_dr = Path(output_dir) / f"{case_id_dr}{label}"
        print(f"  CT-Verzeichnis : {case_out_dr / 'CT'}")
        print(f"  RTSTRUCT       : {case_out_dr / f'RS{label}.dcm'}")
        return {
            "case_id":    Path(case_dir).resolve().name,
            "dry_run":    True,
            "T":          T.tolist(),
            "rotation_center": center.tolist(),
            "rotation_center_label": resolved_label,
        }

    # ── 3b. FoR-Strategie ────────────────────────────────────────────────────
    if new_frame_of_reference:
        new_for_uid = str(generate_uid())
        for_strategy = "new"
        print(f"  FoR-Strategie: NEU ({new_for_uid[:24]}…)")
    else:
        new_for_uid = None  # save_ct_series und transform_rtstruct lassen FoR unveraendert
        for_strategy = "keep"
        print(
            "  FoR-Strategie: KEEP (alte FoR wird beibehalten).\n"
            "    Hinweis: Vorhandene RTPLAN/RTDOSE mit derselben FoR werden im\n"
            "    TPS automatisch auch auf das transformierte CT ueberlagert."
        )

    # ── 4. CT transformieren + speichern ─────────────────────────────────────
    if method == "metadata":
        print("\nAktualisiere DICOM-Metadaten (HU-Werte exakt erhalten) …")
        out_slices = mod.apply_metadata_transform(slices, T)
        save_info  = mod.save_ct_series(
            out_slices, str(ct_out),
            series_description_suffix=label,
            frame_of_reference_uid=new_for_uid,
        )
    elif method == "resample":
        print(f"\nNeuabtastung (Interpolationsordnung {order}) …")
        new_volume = mod.resample_volume(volume_hu, geom["affine"], T, order=order)
        save_info  = mod.save_ct_series(
            slices, str(ct_out), new_volume_hu=new_volume,
            series_description_suffix=label,
            frame_of_reference_uid=new_for_uid,
        )
    else:
        raise ValueError(f"Unbekannte Methode: {method!r}")

    # SeriesNumber-Offset auf alle CT-Slices anwenden (unabhaengig von Methode)
    if series_number_offset:
        for fname in sorted(os.listdir(str(ct_out))):
            if not fname.endswith(".dcm"):
                continue
            fpath = os.path.join(str(ct_out), fname)
            ds_ct = pydicom.dcmread(fpath)
            try:
                sn = int(getattr(ds_ct, "SeriesNumber", 0) or 0)
                ds_ct.SeriesNumber = sn + series_number_offset
                ds_ct.save_as(fpath)
            except (TypeError, ValueError):
                pass

    # ── 5. RTSTRUCT transformieren (Stage 2) ─────────────────────────────────
    # Drehpunkt-Position im transformierten Koordinatensystem.  Da die Rotation
    # das gewaehlte Zentrum invariant laesst, gilt: T(centre) = centre + (tx,ty,tz).
    drehpunkt_pos = center + np.array([tx, ty, tz])

    print("\nTransformiere RTSTRUCT …")
    description = build_transform_description(
        tx, ty, tz, rx, ry, rz,
        center_label=resolved_label, method=method, for_strategy=for_strategy,
    )
    new_rs = transform_rtstruct(
        rs_ds=rs_ds,
        T=T,
        sop_map=save_info["sop_map"],
        new_ct_series_uid=save_info["series_uid"],
        new_for_uid=new_for_uid,
        drehpunkt_position=drehpunkt_pos,
        label_suffix=label,
        description=description,
        series_number_offset=series_number_offset,
    )
    new_rs.save_as(str(rs_out))
    print(f"  RTSTRUCT geschrieben -> {rs_out}")
    print(f"  Drehpunkt-Marker eingefuegt bei "
          f"({drehpunkt_pos[0]:.2f}, {drehpunkt_pos[1]:.2f}, {drehpunkt_pos[2]:.2f}) mm")

    # ── 6. Optional: Verifikation der Centroid-Linearitaet ───────────────────
    verify_report: dict | None = None
    if verify:
        verify_report = _verify_rs_centroids(rs_ds, str(rs_out), T)

    print("\nFertig.")
    return {
        "case_id":         case_id,
        "output_dir":      str(case_out),
        "ct_output_dir":   str(ct_out),
        "rs_output_path":  str(rs_out),
        "ct_series_uid":   save_info["series_uid"],
        "ct_for_uid":      save_info["frame_of_reference_uid_used"],
        "sop_map":         save_info["sop_map"],
        "rs_sop_uid":      str(new_rs.SOPInstanceUID),
        "rs_series_uid":   str(new_rs.SeriesInstanceUID),
        "rotation_center": center.tolist(),
        "rotation_center_label": resolved_label,
        "drehpunkt_pos":   drehpunkt_pos.tolist(),
        "verify":          verify_report,
    }


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
            "  python -m dicom_file_modifier.case_modifier data/0000000171\n"
            "  python -m dicom_file_modifier.case_modifier data/0000000171 "
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

    grp_r = p.add_argument_group("Rotation [deg]  -  extrinsisch XYZ um Volumenmitte")
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
    """
    Bestimmt aus den CLI-Argumenten + dem RS das gewuenschte Rotationszentrum.

    Rueckgabe: (Position oder None fuer Volumenzentrum, Label fuer Logging).
    """
    if args.center is not None:
        pos = parse_center_spec(args.center, rs_ds, volume_center_lps)
        if args.center.lower().strip() == "volume":
            return None, "Volumenmitte"
        if args.center.lower().strip().startswith("marker:"):
            return pos, f"Marker '{args.center.split(':', 1)[1].strip()}'"
        return pos, "manuell"

    if args.non_interactive or not sys.stdin.isatty():
        return None, "Volumenmitte"

    pos = interactive_center_prompt(rs_ds, volume_center_lps)
    if np.allclose(pos, volume_center_lps, atol=1e-9):
        return None, "Volumenmitte"
    return pos, "interaktiv"


def _run_self_test(args: argparse.Namespace) -> int:
    """
    Identitaets-Transform end-to-end + Pruefung, dass alle ContourData-Werte
    bis auf Float-Round-Trip-Rauschen mit dem Original uebereinstimmen.
    """
    import tempfile

    print(f"\n--- Self-Test (Identitaets-Transform) auf '{args.case_dir}' ---")
    with tempfile.TemporaryDirectory(prefix="case_modifier_selftest_") as tmp:
        info = run_case_transform(
            case_dir=args.case_dir,
            output_dir=tmp,
            tx=0, ty=0, tz=0, rx=0, ry=0, rz=0,
            method="metadata",
            label=args.label,
            rs_override=args.rs_override,
            new_frame_of_reference=False,
            verify=True,
        )
        # Konturen Original vs. neu vergleichen
        from .analyzer import (
            load_rtstruct, extract_contours, get_structure_names,
        )
        _, rs_path = discover_case(args.case_dir, rs_override=args.rs_override)
        orig = load_rtstruct(str(rs_path))
        new  = load_rtstruct(info["rs_output_path"])
        names_o = get_structure_names(orig)

        max_dev = 0.0
        worst   = None
        for roi_num, roi_name in names_o.items():
            c_o = extract_contours(orig, roi_num)
            c_n = extract_contours(new, roi_num)
            if not c_o or not c_n or len(c_o) != len(c_n):
                continue
            for a, b in zip(c_o, c_n):
                if a.shape != b.shape:
                    continue
                d = float(np.max(np.abs(a - b)))
                if d > max_dev:
                    max_dev = d
                    worst = roi_name
        print(f"  Max Konturpunkt-Abweichung: {max_dev:.3e} mm  (worst: {worst})")
        passed = max_dev < 1e-4
        print(f"  Ergebnis: {'PASS' if passed else 'FAIL'}")
    return 0 if passed else 1


def main(argv: "list[str] | None" = None) -> int:
    args = _build_parser().parse_args(argv)
    try:
        # --self-test: dedizierter Code-Pfad, kein User-Output
        if args.self_test:
            return _run_self_test(args)

        # --list-markers: nur RS oeffnen und Marker auflisten, dann beenden.
        if args.list_markers:
            _, rs_path = discover_case(args.case_dir, rs_override=args.rs_override)
            rs_ds = pydicom.dcmread(str(rs_path))
            print_marker_table(find_point_markers(rs_ds))
            return 0

        # Center-Resolution braucht das RS und das Volumenzentrum.  Damit
        # der interaktive Prompt sinnvolle Defaults zeigen kann, laden wir
        # CT (lite, ohne Pixel) + RS bereits hier.
        ct_dir, rs_path = discover_case(args.case_dir, rs_override=args.rs_override)
        rs_ds = pydicom.dcmread(str(rs_path))

        # CT-Metadaten ohne Pixel laden -> Volumenzentrum berechnen.
        import glob
        ct_files = sorted(glob.glob(os.path.join(str(ct_dir), "*.dcm")))
        meta_slices = []
        for f in ct_files:
            try:
                ds = pydicom.dcmread(f, stop_before_pixels=True)
                if hasattr(ds, "ImagePositionPatient"):
                    meta_slices.append(ds)
            except Exception:
                pass
        if not meta_slices:
            raise ValueError(f"Keine gueltigen CT-Slices in {ct_dir}")
        meta_slices.sort(key=lambda s: float(s.ImagePositionPatient[2]))
        geom_lite = mod.extract_geometry(meta_slices)
        vol_c = mod.volume_center(geom_lite)

        center, center_label = _resolve_center_from_args(args, rs_ds, vol_c)

        run_case_transform(
            case_dir=args.case_dir,
            output_dir=args.output,
            tx=args.tx, ty=args.ty, tz=args.tz,
            rx=args.rx, ry=args.ry, rz=args.rz,
            method=args.method, order=args.order,
            label=args.label,
            rs_override=args.rs_override,
            center=center,
            center_label=center_label,
            new_frame_of_reference=args.new_frame_of_reference,
            dry_run=args.dry_run,
            verify=args.verify,
        )
    except (FileNotFoundError, ValueError, KeyError) as e:
        print(f"\nFehler: {e}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
