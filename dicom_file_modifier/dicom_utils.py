"""
dicom_utils.py - Kleine DICOM-Helfer ohne Abhaengigkeit von den CLI-Modulen.

Mehrere Module brauchen diese Funktionen (``case_modifier``, ``modifier``,
``dose``, ``rtstruct_writer``, ``analyzer``).  Frueher lagen sie in
``case_modifier`` bzw. ``modifier`` (``find_point_markers`` bis P0.7 in
``case_modifier``); dort werden sie weiter re-exportiert, damit bestehende
Importe gueltig bleiben (Plan P0.3).  Das Modul importiert nichts aus dem Paket.
"""

from __future__ import annotations

import numpy as np
import pydicom


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


def set_sop_instance_uid(ds: pydicom.Dataset, uid: str) -> None:
    """
    Setzt die SOPInstanceUID und zieht die MediaStorageSOPInstanceUID im
    File-Meta-Header (DICOM Part 10) mit.  Beide müssen identisch sein, und
    ``save_as`` gleicht sie nicht ab.
    """
    ds.SOPInstanceUID = uid
    # Ohne Header (fehlt oder leer) keinen unvollständigen Header anlegen
    if getattr(ds, "file_meta", None):
        ds.file_meta.MediaStorageSOPInstanceUID = uid


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
