"""App-Einstellungen: Ergebnis-Stammordner und Job-Ordner (QSettings, pro Benutzer)."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

from PySide6.QtCore import QSettings

APP_NAME = "DICOM-RT-Toolkit"


def default_results_root() -> Path:
    return Path.home() / APP_NAME / "Ergebnisse"


def default_jobs_dir() -> Path:
    base = os.environ.get("LOCALAPPDATA") or str(Path.home() / "AppData" / "Local")
    return Path(base) / APP_NAME / "jobs"


@dataclass
class AppConfig:
    results_root: Path
    jobs_dir: Path

    @classmethod
    def load(cls) -> "AppConfig":
        s = QSettings(APP_NAME, APP_NAME)
        root = s.value("results_root", "") or str(default_results_root())
        return cls(results_root=Path(root), jobs_dir=default_jobs_dir())

    def save(self) -> None:
        QSettings(APP_NAME, APP_NAME).setValue("results_root", str(self.results_root))
