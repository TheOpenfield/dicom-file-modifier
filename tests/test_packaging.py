"""Lizenzhinweise Dritter (packaging/third_party_notices.py): der Generator laeuft ohne Build und ohne PyInstaller."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "packaging"))

from third_party_notices import notices_text  # noqa: E402


def test_notices_list_the_bundled_packages_and_carry_the_lgpl_texts():
    text = notices_text({"numpy", "shapely", "pydicom", "PySide6", "libcrypto-3-x64"}, app_version="9.9.9")
    assert text.startswith("DICOM-RT-Toolkit 9.9.9")
    packages = text.split("2. Native components")[0]
    assert "numpy" in packages and "shapely" in packages and "pydicom" in packages
    assert "pytest" not in packages and "dicom-file-modifier" not in packages.split("1. Python packages")[1]
    for needle in ("GEOS", "Qt 6", "PySide6", "GNU LESSER GENERAL PUBLIC LICENSE", "Version 3, 29 June 2007",
                   "GNU GENERAL PUBLIC LICENSE", "Version 2.1, February 1999", "Apache License",
                   "LICENSE_GEOS"):
        assert needle in text, needle
    assert "LicenseRef-Qt-Commercial" not in text
