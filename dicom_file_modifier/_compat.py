"""
_compat.py - Kleine Kompatibilitaets-Shims fuer pydicom 2.x/3.x.

Nur Funktionen, die auf beiden Versionen identische Dateien erzeugen sollen.
Der Rest des Pakets ruft pydicom weiterhin direkt auf; die Umstellung der
bestehenden Schreibpfade folgt mit dem Umgebungs-Upgrade (siehe Plan P0.2).
"""

from __future__ import annotations

import pydicom
from pydicom.uid import ImplicitVRLittleEndian

PYDICOM_MAJOR = int(str(pydicom.__version__).split(".")[0])


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
