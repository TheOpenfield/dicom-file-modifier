"""
Befunde fuer Oberflaechen: ``Issue``/``UserInputError`` aus dem Kern und die
Uebersetzung bekannter Ausnahmen in deutsche Klartexte mit "Was tun?".

``exit_code_for`` ordnet eine Ausnahme dem Worker-Exit-Code zu: 2 = Eingabe
(vom Nutzer behebbar), 3 = abgebrochen, 1 = interner Fehler.
"""

from __future__ import annotations

import traceback

from .. import _runtime
from ..issues import LEVELS, Issue, UserInputError

__all__ = ["LEVELS", "Issue", "UserInputError", "issue_from_exception", "exit_code_for",
           "is_input_error", "has_errors", "worst_level"]

_DATASET_ATTR = "Dataset' object has no attribute"      # pydicom: fehlendes Pflicht-Tag


def _invalid_dicom_error():
    from pydicom.errors import InvalidDicomError
    return InvalidDicomError


def is_input_error(exc: BaseException) -> bool:
    """Vom Nutzer behebbar: ungueltige Eingaben, fehlende/gesperrte Dateien, kaputte DICOMs."""
    if isinstance(exc, (ValueError, KeyError, OSError, _invalid_dicom_error())):
        return True
    return isinstance(exc, AttributeError) and _DATASET_ATTR in str(exc)


def exit_code_for(exc: BaseException) -> int:
    if isinstance(exc, _runtime.JobCancelled):
        return 3
    return 2 if is_input_error(exc) else 1


def _detail(exc: BaseException) -> str:
    return f"{type(exc).__name__}: {exc}"


def issue_from_exception(exc: BaseException, *, field: str = "") -> Issue:
    """Ausnahme -> ``Issue`` (Klartext, Hinweis, Details); ``UserInputError`` behaelt sein Issue."""
    if isinstance(exc, UserInputError):
        return exc.issue
    if isinstance(exc, _runtime.JobCancelled):
        return Issue("info", "JOB.CANCELLED", "Abgebrochen.",
                     hint_de="Es wurden keine Ergebnisse uebernommen.")
    if isinstance(exc, _invalid_dicom_error()):
        return Issue("error", "DICOM.INVALID", f"Keine gueltige DICOM-Datei ({exc}).",
                     hint_de="Datei pruefen oder aus dem Planungssystem neu exportieren.",
                     field=field, detail=_detail(exc))
    if isinstance(exc, FileNotFoundError):
        msg = (f"Datei oder Ordner nicht gefunden: {exc.filename}" if exc.filename
               else str(exc))
        return Issue("error", "FILE.NOT_FOUND", msg,
                     hint_de="Pfad pruefen; ist ein Netzlaufwerk getrennt oder die Datei verschoben?",
                     field=field, detail=_detail(exc))
    if isinstance(exc, PermissionError):
        return Issue("error", "FILE.PERMISSION",
                     f"Kein Zugriff auf {exc.filename or 'eine Datei'}.",
                     hint_de="Ist die Datei in einem anderen Programm geoeffnet (z.B. die Sammel-CSV "
                             "in Excel) oder schreibgeschuetzt? Schliessen und erneut starten.",
                     field=field, detail=_detail(exc))
    if isinstance(exc, OSError):
        return Issue("error", "FILE.IO", f"Dateifehler: {exc}",
                     hint_de="Freien Speicherplatz und die Pfadlaenge pruefen (Windows: 260 Zeichen).",
                     field=field, detail=_detail(exc))
    if isinstance(exc, MemoryError):
        return Issue("error", "SYSTEM.MEMORY", "Nicht genug Arbeitsspeicher.",
                     hint_de="Andere Programme schliessen oder ein groeberes Raster waehlen.",
                     field=field, detail=_detail(exc))
    if isinstance(exc, AttributeError) and _DATASET_ATTR in str(exc):
        tag = str(exc).rsplit("attribute", 1)[-1].strip(" '\"")
        return Issue("error", "DICOM.MISSING_TAG", f"Pflichtangabe {tag} fehlt in einer DICOM-Datei.",
                     hint_de="Datei aus dem Planungssystem neu exportieren.",
                     field=field, detail=_detail(exc))
    if isinstance(exc, KeyError):
        return Issue("error", "INPUT.MISSING_KEY", f"Erwarteter Eintrag fehlt: {exc}",
                     field=field, detail=_detail(exc))
    if isinstance(exc, ValueError):
        return Issue("error", "INPUT.INVALID", str(exc), field=field, detail=_detail(exc))
    return Issue("error", "INTERNAL", f"Interner Fehler ({type(exc).__name__}: {exc}).",
                 hint_de="Bitte das Diagnosepaket erstellen und melden.", field=field,
                 detail="".join(traceback.format_exception(type(exc), exc, exc.__traceback__)))


def has_errors(issues) -> bool:
    return any(i.level == "error" for i in issues)


def worst_level(issues) -> str:
    """``error`` > ``warning`` > ``info``; ohne Issues ``""``."""
    levels = {i.level for i in issues}
    return next((lv for lv in LEVELS if lv in levels), "")
