"""Einstieg der Desktop-App: ``dfm gui [DATENSATZ]`` bzw. ``dicom-rt-toolkit [DATENSATZ]``."""

from __future__ import annotations

import argparse
import sys
from typing import Optional


def main(argv: Optional[list] = None) -> int:
    parser = argparse.ArgumentParser(prog="dfm gui", description="Desktop-App starten.")
    parser.add_argument("case", nargs="?", metavar="DATENSATZ", help="Datensatz-Ordner direkt öffnen")
    parser.add_argument("--smoke-test", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if sys.platform == "win32":
        _set_app_user_model_id()

    from PySide6.QtWidgets import QApplication

    from ..api.jobs import cleanup_staging
    from .config import APP_NAME, AppConfig
    from .icons import app_icon
    from .window import MainWindow

    app = QApplication.instance() or QApplication(sys.argv[:1])
    app.setApplicationName(APP_NAME)
    app.setWindowIcon(app_icon())                    # Dialoge, Taskleiste, Alt-Tab
    if args.smoke_test:
        return _smoke_test(app, args.case)
    config = AppConfig.load()
    cleanup_staging(config.results_root)             # Reste eines abgestuerzten Laufs
    win = MainWindow(config)
    win.show()
    if args.case:
        win.open_case(args.case)
    return app.exec()


def _set_app_user_model_id() -> None:
    """Beim Start aus Python zeigt die Taskleiste sonst das Python-Symbol (die EXE braucht das nicht)."""
    import ctypes
    try:
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("DicomRtToolkit.App")
    except (AttributeError, OSError):
        pass


def _smoke_test(app, case) -> int:
    """
    Pruefung einer Installation (auch der gefrorenen App): Datensatz oeffnen, alle
    Seiten pruefen lassen und die Strukturanalyse ueber den Worker rechnen,
    Ergebnisse in einem Temp-Ordner.
    Exit 0 ok, 1 Lauf fehlgeschlagen, 2 kein Datensatz, 3 Pruefung nicht bestanden, 4 Zeitueberschreitung,
    5 Symbole fehlen (Assets oder Qt6Svg nicht im Bundle).
    """
    import shutil
    import tempfile
    import time
    from pathlib import Path

    from .config import AppConfig
    from .window import MainWindow

    if not case:
        return 2
    tmp = Path(tempfile.mkdtemp(prefix="dfm_smoke_"))
    win = MainWindow(AppConfig(results_root=tmp / "results", jobs_dir=tmp / "jobs"))
    page = win.structures

    def wait(pred, timeout: float) -> bool:
        end = time.monotonic() + timeout
        while not pred():
            if time.monotonic() > end:
                return False
            app.processEvents()
            time.sleep(0.02)
        return True

    try:
        win.show()
        if win.windowIcon().isNull() or any(win.nav.item(i).icon().pixmap(22, 22).isNull()
                                            for i in range(win.nav.count())):
            return 5
        win.open_case(case)
        if not wait(lambda: all(p.start_button.isEnabled() for p in win.pages), 120):
            return 3
        page.start_button.click()
        if not wait(lambda: page.last_result is not None, 600):
            return 4
        return 0 if page.last_result["status"] in ("ok", "ok_warnings") else 1
    finally:
        win.close()
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
