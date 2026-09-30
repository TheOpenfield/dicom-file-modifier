"""Einstieg der Desktop-App: ``dfm gui [FALL]`` bzw. ``dicom-rt-toolkit [FALL]``."""

from __future__ import annotations

import argparse
import sys
from typing import Optional


def main(argv: Optional[list] = None) -> int:
    parser = argparse.ArgumentParser(prog="dfm gui", description="Desktop-App starten.")
    parser.add_argument("case", nargs="?", help="Fallordner direkt öffnen")
    args = parser.parse_args(argv)

    from PySide6.QtWidgets import QApplication

    from ..api.jobs import cleanup_staging
    from .config import APP_NAME, AppConfig
    from .window import MainWindow

    app = QApplication.instance() or QApplication(sys.argv[:1])
    app.setApplicationName(APP_NAME)
    config = AppConfig.load()
    cleanup_staging(config.results_root)             # Reste eines abgestuerzten Laufs
    win = MainWindow(config)
    win.show()
    if args.case:
        win.open_case(args.case)
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
