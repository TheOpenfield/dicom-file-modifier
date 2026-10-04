"""
dicom_utils.py - Kleine DICOM-Helfer ohne Abhaengigkeit von den CLI-Modulen.

Mehrere Module brauchen diese Funktionen (``case_modifier``, ``modifier``,
``dose``, ``rtstruct_writer``, ``analyzer``, ``demo``).  Das Modul importiert
nichts aus dem Paket.
"""

from __future__ import annotations

import numpy as np
import pydicom
from pydicom.uid import ExplicitVRLittleEndian


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


CONTOUR_DATA_TAG = 0x30060050


def contour_points(contour: pydicom.Dataset) -> np.ndarray:
    """
    ``ContourData`` einer Kontur als ``(N, 3)``-float64-Array, bitgleich zu
    ``np.asarray(contour.ContourData, dtype=np.float64)``.  Solange das Element
    noch roh ist, werden die Bytes direkt gelesen: pydicom legt sonst je Wert ein
    ``DSfloat``-Objekt an (langsam, viel Speicher) und behaelt es im Dataset.
    Alles Unerwartete (keine ASCII-Bytes, ungueltige Werte, Anzahl kein
    Vielfaches von 3) nimmt den pydicom-Weg; ein fehlendes Element wirft wie dort
    ``AttributeError``.
    """
    elem = contour.get_item(CONTOUR_DATA_TAG)
    raw = elem.value if elem is not None else None
    if isinstance(raw, (bytes, bytearray)) and raw:
        try:
            vals = np.array(raw.decode("ascii").split("\\"), dtype=np.float64)
        except (UnicodeDecodeError, ValueError):
            vals = None
        if vals is not None and vals.size % 3 == 0:
            return vals.reshape(-1, 3)
    return np.asarray(contour.ContourData, dtype=np.float64).reshape(-1, 3)


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
            pts = contour_points(c)
            if pts.shape[0] == 0:
                continue
            roi_num = int(getattr(rc, "ReferencedROINumber", -1))
            markers.append((name_map.get(roi_num, f"ROI#{roi_num}"), pts[0].copy()))
            break  # erster POINT je ROI reicht
    return markers


# Tags, deren Werte vom alten Pixelinhalt bzw. von der alten PixelRepresentation
# (VR US/SS) abhaengen und nach dem Neuschreiben der Pixel nicht mehr stimmen.
_STALE_PIXEL_TAGS = ("SmallestImagePixelValue", "LargestImagePixelValue",
                     "PixelPaddingValue", "PixelPaddingRangeLimit")


def replace_pixel_data(ds: pydicom.Dataset, stored: np.ndarray) -> None:
    """
    Ersetzt die Pixeldaten von ``ds`` durch das 2D-Array ``stored`` als
    vorzeichenbehaftetes int16 (unkomprimiert, VR OW).

    Eine komprimierte oder Big-Endian-Transfer-Syntax der Quelle wird auf
    Explicit VR Little Endian umgestellt (sonst entstuende beim Speichern ein
    Absturz bzw. eine inkonsistente Datei); ``Smallest/LargestImagePixelValue``
    und ``PixelPadding*`` werden entfernt, weil sie nach dem Resampling nicht
    mehr gelten.  Rows/Columns/Bits*/PixelRepresentation werden gesetzt.
    """
    arr = np.ascontiguousarray(stored, dtype=np.int16)
    fm = getattr(ds, "file_meta", None)
    ts = fm.get("TransferSyntaxUID") if fm is not None else None
    if ts is not None and (ts.is_compressed or not ts.is_little_endian):
        fm.TransferSyntaxUID = ExplicitVRLittleEndian
    for kw in _STALE_PIXEL_TAGS:
        if kw in ds:
            del ds[kw]
    ds.PixelData = arr.tobytes()
    elem = ds["PixelData"]
    elem.VR = "OW"
    elem.is_undefined_length = False
    ds.Rows, ds.Columns = arr.shape
    ds.BitsAllocated = 16
    ds.BitsStored = 16
    ds.HighBit = 15
    ds.PixelRepresentation = 1
