"""
_compat.py - Kleine Kompatibilitaets-Shims fuer pydicom 2.x/3.x.

Nur Funktionen, die auf beiden Versionen identische Dateien erzeugen sollen.
Genutzt von ``demo``, ``rtstruct_writer`` (``dcmwrite_file_format``) und
``modifier`` (``replace_pixel_data``).  Das Paket verlangt seit P0.2
pydicom >= 3; die Zweige fuer pydicom 2 halten die alte Golden-Umgebung
(``tools/legacy/requirements-py38.txt``) lauffaehig.
"""

from __future__ import annotations

import numpy as np
import pydicom
from pydicom.uid import ExplicitVRLittleEndian, ImplicitVRLittleEndian

PYDICOM_MAJOR = int(str(pydicom.__version__).split(".")[0])

# Tags, deren Werte vom alten Pixelinhalt bzw. von der alten PixelRepresentation
# (VR US/SS) abhaengen und nach dem Neuschreiben der Pixel nicht mehr stimmen.
_STALE_PIXEL_TAGS = ("SmallestImagePixelValue", "LargestImagePixelValue",
                     "PixelPaddingValue", "PixelPaddingRangeLimit")


def replace_pixel_data(ds: pydicom.Dataset, stored: np.ndarray) -> None:
    """
    Ersetzt die Pixeldaten von ``ds`` durch das 2D-Array ``stored`` als
    vorzeichenbehaftetes int16 (unkomprimiert, VR OW).

    Eine komprimierte oder Big-Endian-Transfer-Syntax der Quelle wird auf
    Explicit VR Little Endian umgestellt (vorher entstand dort ein Absturz beim
    Speichern bzw. eine inkonsistente Datei); ``Smallest/LargestImagePixelValue``
    und ``PixelPadding*`` werden entfernt, weil sie nach dem Resampling nicht
    mehr gelten.  Rows/Columns/Bits*/PixelRepresentation werden gesetzt.
    """
    arr = np.ascontiguousarray(stored, dtype=np.int16)
    fm = getattr(ds, "file_meta", None)
    ts = fm.get("TransferSyntaxUID") if fm is not None else None
    if ts is not None and (ts.is_compressed or not ts.is_little_endian):
        fm.TransferSyntaxUID = ExplicitVRLittleEndian
        if PYDICOM_MAJOR < 3:
            ds.is_implicit_VR = False
            ds.is_little_endian = True
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


def dcmwrite_file_format(path, ds: pydicom.Dataset) -> None:
    """
    Schreibt ``ds`` im DICOM-File-Format (Praeambel + File-Meta), kodiert
    gemaess ``ds.file_meta.TransferSyntaxUID`` (nur Little Endian).

    pydicom 2.x braucht die Kodierflags am Dataset und ``write_like_original=False``;
    pydicom 3.x leitet die Kodierung aus der Transfer Syntax ab und kennt dafuer
    ``enforce_file_format=True`` (die Flag-Setter sind dort deprecated).
    """
    ts = ds.file_meta.TransferSyntaxUID
    if PYDICOM_MAJOR < 3:
        ds.is_little_endian = True
        ds.is_implicit_VR = (ts == ImplicitVRLittleEndian)
        pydicom.dcmwrite(str(path), ds, write_like_original=False)
    else:
        pydicom.dcmwrite(str(path), ds, enforce_file_format=True)
