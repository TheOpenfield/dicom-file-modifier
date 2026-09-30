# PyInstaller-Spec (Plan: Packaging): ein onedir-Ordner mit zwei EXEs.
#   DICOM-RT-Toolkit.exe   Desktop-App (fensterlos)
#   dfm.exe                CLI und Worker der App (Konsole, ohne Qt)
# Bauen und pruefen (README: Desktop app):
#   uv sync --extra gui --group build
#   .venv/Scripts/python -m PyInstaller packaging/dfm.spec --noconfirm
#   .venv/Scripts/python packaging/check_bundle.py
# (ein einfaches "uv sync" entfernt Extra und Gruppe wieder)
# Kein onefile (entpackt bei jedem Start nach %TEMP%), kein UPX (Virenscanner).
# Neben den EXEs liegen LICENSE.txt und THIRD-PARTY-NOTICES.txt (third_party_notices.py).

import shutil
import sys
from pathlib import Path

from PyInstaller.utils.hooks import collect_submodules, copy_metadata

ROOT = Path(SPECPATH).parent
PKG = "dicom_file_modifier"
LIBS = ("numpy", "scipy", "pydicom", "shapely", "matplotlib", "scikit-image", "plotly")   # cli.versions()
EXCLUDES = ["tkinter", "IPython", "pandas", "pytest", "PyInstaller",
            "PySide6.QtWebEngineCore", "PySide6.QtWebEngineWidgets", "PySide6.QtQml",
            "PySide6.QtQuick", "PySide6.Qt3DCore", "PySide6.QtCharts", "PySide6.QtMultimedia",
            "PySide6.QtPdf", "pyqtgraph", "OpenGL"]
QT = ["PySide6", "shiboken6"]
metadata = [d for lib in LIBS for d in copy_metadata(lib)]
ASSETS = [(str(ROOT / PKG / "gui" / "assets"), f"{PKG}/gui/assets")]   # App-Icon und SVG-Symbole
ICON = str(ROOT / PKG / "gui" / "assets" / "app.ico")

# cli.py und api.workflow() importieren die Werkzeuge per importlib
cli = Analysis([str(ROOT / "packaging" / "entry_dfm.py")], pathex=[str(ROOT)],
               hiddenimports=[m for m in collect_submodules(PKG) if not m.startswith(f"{PKG}.gui")]
               + ["mpl_toolkits.mplot3d"],
               datas=metadata, excludes=EXCLUDES + QT)
gui = Analysis([str(ROOT / "packaging" / "entry_gui.py")], pathex=[str(ROOT)],
               hiddenimports=collect_submodules(f"{PKG}.gui") + collect_submodules(f"{PKG}.api"),
               datas=metadata + ASSETS, excludes=EXCLUDES)

# Keine DICOM-Dateien im Bundle (check_bundle.py prueft das).  pydicom bringt Testdateien
# mit; von pydicom/data/ bleiben nur die JSON-Dateien: pydicom.examples sucht seine
# Beispiele beim Import und liest dabei urls.json (ohne Download, die Namen stehen nicht darin).
def _keep(dest: str) -> bool:
    d = dest.replace("\\", "/")
    if d.startswith("pydicom/data/"):
        return d.endswith(".json") and d.count("/") == 2
    return not d.lower().endswith(".dcm")


for a in (cli, gui):
    a.datas = [d for d in a.datas if _keep(d[0])]

cli_exe = EXE(PYZ(cli.pure), cli.scripts, [], exclude_binaries=True, name="dfm",
              console=True, upx=False, icon=ICON)
gui_exe = EXE(PYZ(gui.pure), gui.scripts, [], exclude_binaries=True, name="DICOM-RT-Toolkit",
              console=False, upx=False, icon=ICON)
COLLECT(gui_exe, gui.binaries, gui.datas, cli_exe, cli.binaries, cli.datas,
        name="DICOM-RT-Toolkit", upx=False)

# Eigene Lizenz und Lizenzhinweise Dritter neben die EXEs (check_bundle.py prueft beide)
sys.path.insert(0, str(ROOT / "packaging"))
from third_party_notices import bundled_names, write_notices
from dicom_file_modifier import __version__

BUNDLE = Path(DISTPATH) / "DICOM-RT-Toolkit"
shutil.copyfile(ROOT / "LICENSE", BUNDLE / "LICENSE.txt")
write_notices(bundled_names(cli, gui), BUNDLE / "THIRD-PARTY-NOTICES.txt", __version__)
