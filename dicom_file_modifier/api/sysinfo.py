"""Systemangaben fuer Pruefungen vor dem Start (ohne Zusatzpakete)."""

from __future__ import annotations

import ctypes
import os
import sys
from typing import Optional


def available_memory_bytes() -> Optional[int]:
    """Verfuegbarer physischer Arbeitsspeicher; ``None``, wenn nicht ermittelbar."""
    if sys.platform == "win32":
        class _MemoryStatusEx(ctypes.Structure):
            _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                        ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                        ("ullTotalPageFile", ctypes.c_ulonglong),
                        ("ullAvailPageFile", ctypes.c_ulonglong),
                        ("ullTotalVirtual", ctypes.c_ulonglong),
                        ("ullAvailVirtual", ctypes.c_ulonglong),
                        ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]

        stat = _MemoryStatusEx()
        stat.dwLength = ctypes.sizeof(_MemoryStatusEx)
        if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(stat)):
            return int(stat.ullAvailPhys)
        return None
    try:
        return int(os.sysconf("SC_AVPHYS_PAGES") * os.sysconf("SC_PAGE_SIZE"))
    except (AttributeError, OSError, ValueError):
        return None


def format_bytes(n: Optional[float]) -> str:
    """``1.23e9`` -> ``'1.2 GB'`` (Dezimalpunkt wie alle Zahlen der Texte und Berichte)."""
    if n is None:
        return "?"
    for unit, size in (("GB", 1e9), ("MB", 1e6), ("kB", 1e3)):
        if n >= size:
            return f"{n / size:.1f} {unit}"
    return f"{int(n)} B"
