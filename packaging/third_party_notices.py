"""
Lizenzhinweise Dritter fuer das Bundle: THIRD-PARTY-NOTICES.txt neben den EXEs.

Der Text entsteht beim Build aus den Metadaten der tatsaechlich gebuendelten
Distributionen (Name, Version, Lizenz, Quelle und die Lizenztexte aus dem
dist-info), dazu die nativen Komponenten (Python, Qt 6, GEOS, OpenBLAS, OpenSSL,
libffi, MSVC-Laufzeit, PyInstaller-Bootloader) und die vollen Texte LGPL 3,
GPL 3, LGPL 2.1 und Apache 2.0 aus packaging/licenses/: Qt 6 und
PySide6/shiboken6 werden unter der LGPL 3 genutzt, GEOS unter der LGPL 2.1.

  python packaging/third_party_notices.py [--toc-dir build/dfm] [-o DATEI]

dfm.spec ruft write_notices() nach COLLECT mit den Namen der Analyse auf; die
Kommandozeile liest die Namen aus den TOC-Dateien eines fertigen Builds.
"""

from __future__ import annotations

import argparse
import datetime as dt
import platform
import re
import sys
from importlib.metadata import Distribution, distribution, packages_distributions
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
LICENSES_DIR = HERE / "licenses"
APP_NAME = "DICOM-RT-Toolkit"
OWN_DISTS = {"dicom-file-modifier"}
SKIP_FILES = ("LicenseRef-Qt-Commercial",)          # Qt-Kommerzlizenz: nicht die genutzte Lizenz
LGPL_CHOICE = {"pyside6_essentials": "LGPL-3.0-only", "pyside6": "LGPL-3.0-only",
               "shiboken6": "LGPL-3.0-only"}         # gewaehlte Lizenz bei "LGPL OR GPL"

# Native Komponenten ausserhalb der Python-Distributionen: (Name, Dateien im Bundle, Lizenz, Quelle)
NATIVE = [
    ("Python", "_internal/python3*.dll, base_library.zip, *.pyd",
     "PSF-2.0 (Python Software Foundation License; text in section 3.2)", "https://www.python.org/"),
    ("Qt 6 Essentials", "_internal/PySide6/Qt6*.dll, plugins/, translations/",
     "LGPL-3.0-only (GNU LGPL v3 and GPL v3 in section 3.3)",
     "https://code.qt.io/ ; third-party code inside Qt: https://doc.qt.io/qt-6/licenses-used-in-qt.html"),
    ("GEOS", "_internal/Shapely.libs/geos*.dll",
     "LGPL-2.1-or-later (GNU LGPL v2.1 in section 3.3)", "https://libgeos.org/"),
    ("OpenBLAS", "_internal/numpy.libs/libscipy_openblas*.dll, _internal/scipy.libs/libscipy_openblas*.dll",
     "BSD-3-Clause (see the numpy and scipy licence texts)", "https://www.openblas.net/"),
    ("OpenSSL", "_internal/libcrypto-3-x64.dll, _internal/libssl-3-x64.dll",
     "Apache-2.0 (text in section 3.3)", "https://www.openssl.org/"),
    ("libffi", "_internal/libffi-8.dll", "MIT", "https://sourceware.org/libffi/"),
    ("Microsoft Visual C++ Runtime",
     "_internal/MSVCP140*.dll, VCRUNTIME140*.dll, ucrtbase.dll, api-ms-win-*.dll",
     "Microsoft Software License Terms (redistributable)",
     "https://learn.microsoft.com/cpp/windows/latest-supported-vc-redist"),
    ("PyInstaller bootloader", "DICOM-RT-Toolkit.exe, dfm.exe (the loader part of both)",
     "GPL-2.0-or-later with the PyInstaller bootloader exception (a bundled program may use any licence)",
     "https://pyinstaller.org/"),
]
GNU_TEXTS = [
    ("GNU Lesser General Public License v3 (Qt 6, PySide6, shiboken6)", "LGPL-3.0.txt"),
    ("GNU General Public License v3 (referenced by the LGPL v3)", "GPL-3.0.txt"),
    ("GNU Lesser General Public License v2.1 (GEOS)", "LGPL-2.1.txt"),
    ("Apache License 2.0 (OpenSSL)", "Apache-2.0.txt"),
]


def _top_level(name: str) -> str:
    return re.split(r"[\\/.]", name, maxsplit=1)[0]


def bundled_names(*analyses) -> set:
    """Top-Level-Namen der Module, Binaries und Daten von PyInstaller-Analysen."""
    names = set()
    for a in analyses:
        for toc in (a.pure, a.binaries, a.datas):
            for entry in toc:
                names.add(_top_level(entry[0]))
    return names


def names_from_toc_dir(toc_dir: Path) -> set:
    """Dieselben Namen aus den TOC-Dateien eines fertigen Builds (build/dfm)."""
    names = set()
    for toc in list(toc_dir.glob("PYZ-*.toc")) + list(toc_dir.glob("COLLECT-*.toc")):
        text = toc.read_text(encoding="utf-8", errors="replace")
        for m in re.finditer(r"^\s*\('([^']+)'", text, re.M):
            names.add(_top_level(m.group(1)))
    return names


def bundled_distributions(names) -> list:
    """Installierte Distributionen, zu denen die Namen gehoeren (ohne das eigene Paket)."""
    mapping = packages_distributions()
    dist_names = set()
    for n in names:
        dist_names.update(mapping.get(n, []))
    out = []
    for dn in sorted(dist_names, key=str.lower):
        if dn.lower().replace("_", "-") in OWN_DISTS:
            continue
        out.append(distribution(dn))
    return out


def _licence(dist: Distribution) -> str:
    md = dist.metadata
    chosen = LGPL_CHOICE.get(md["Name"].lower())
    expr = md.get("License-Expression")
    if chosen:
        return f"{chosen} (chosen from: {expr or md.get('License')})"
    if expr:
        return expr
    classifiers = [c.split("::")[-1].strip() for c in md.get_all("Classifier", []) if c.startswith("License ::")]
    if classifiers:
        return "; ".join(dict.fromkeys(classifiers))
    lic = (md.get("License") or "").strip()
    if not lic:
        return "see licence text"
    first = lic.splitlines()[0].strip()
    return first[:70] + (" ..." if len(lic) > len(first) else "")


def _source(dist: Distribution) -> str:
    md = dist.metadata
    urls = md.get_all("Project-URL") or []
    for pat in ("source", "repository", "homepage", "home"):
        for u in urls:
            label, _, url = u.partition(",")
            if pat in label.lower():
                return url.strip()
    return md.get("Home-page") or ""


def licence_files(dist: Distribution) -> list:
    """(Pfad im dist-info, Text) der Lizenzdateien: PEP-639-Eintraege, sonst LICENSE*/COPYING*/NOTICE*."""
    declared = dist.metadata.get_all("License-File") or []
    out = []
    for f in dist.files or []:
        s = str(f).replace("\\", "/")
        if ".dist-info/" not in s:
            continue
        inside = s.split(".dist-info/", 1)[1]
        if declared:
            wanted = inside in declared or inside.removeprefix("licenses/") in declared
        else:
            wanted = re.search(r"(LICEN[CS]E|COPYING|NOTICE)", inside, re.I) is not None
        if not wanted or any(skip in inside for skip in SKIP_FILES):
            continue
        try:
            out.append((inside, Path(dist.locate_file(f)).read_text(encoding="utf-8", errors="replace")))
        except OSError:
            pass
    return out


def notices_text(names, app_version: str) -> str:
    dists = bundled_distributions(names)
    lines = []
    w = lines.append
    w(f"{APP_NAME} {app_version} - Third-party notices")
    w("=" * 72)
    w("")
    w(f"Generated {dt.date.today().isoformat()} from the build environment: "
      f"Python {platform.python_version()}, {platform.system()} {platform.machine()}.")
    w("")
    w(f"{APP_NAME} (Python package dicom-file-modifier) is licensed under the MIT License, see LICENSE.txt.")
    w("The program folder also contains the components listed below. Each is the property of its authors")
    w("and is redistributed under its own licence. The licence texts follow in section 3.")
    w("")
    w("LGPL components: Qt 6 and PySide6/shiboken6 are used under the GNU Lesser General Public License v3,")
    w("GEOS under the GNU Lesser General Public License v2.1. They are shipped unmodified as separate,")
    w("dynamically loaded DLLs in the _internal folder, so a user may replace them with other versions of")
    w("the same libraries (relinking as permitted by the LGPL). Source code: https://code.qt.io/ (Qt, PySide6),")
    w("https://libgeos.org/ (GEOS).")
    w("")
    w("1. Python packages")
    w("-" * 72)
    rows = [(d.metadata["Name"], d.version, _licence(d), _source(d)) for d in dists]
    name_w = max((len(r[0]) for r in rows), default=8)
    ver_w = max((len(r[1]) for r in rows), default=7)
    for name, ver, lic, src in rows:
        w(f"{name:<{name_w}}  {ver:<{ver_w}}  {lic}")
        if src:
            w(f"{'':<{name_w}}  {'':<{ver_w}}  {src}")
    w("")
    w("2. Native components")
    w("-" * 72)
    for name, files, lic, src in NATIVE:
        w(f"{name}")
        w(f"    files:   {files}")
        w(f"    licence: {lic}")
        w(f"    source:  {src}")
    w("")
    w("3. Licence texts")
    w("-" * 72)
    w("")
    w("3.1 Python packages (as shipped in each wheel's dist-info)")
    for d in dists:
        label = f"{d.metadata['Name']} {d.version}"
        files = licence_files(d)
        if not files:
            w("")
            w(f"==== {label}: no licence text in the wheel; licence {_licence(d)} ====")
            continue
        for inside, text in files:
            w("")
            w(f"==== {label}: {inside} ====")
            w(text.rstrip())
    w("")
    w("3.2 Python")
    py_licence = Path(sys.base_prefix) / "LICENSE.txt"
    w("")
    if py_licence.is_file():
        w(f"==== Python {platform.python_version()}: LICENSE.txt ====")
        w(py_licence.read_text(encoding="utf-8", errors="replace").rstrip())
    else:
        w(f"==== Python {platform.python_version()}: licence text not found in the build environment; "
          "see https://docs.python.org/3/license.html ====")
    w("")
    w("3.3 GNU and Apache licences")
    for title, fname in GNU_TEXTS:
        w("")
        w(f"==== {title} ====")
        w((LICENSES_DIR / fname).read_text(encoding="utf-8").rstrip())
    w("")
    return "\n".join(lines)


def write_notices(names, out: Path, app_version: str) -> str:
    text = notices_text(names, app_version)
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf-8")
    return text


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="THIRD-PARTY-NOTICES.txt aus den TOC-Dateien eines Builds erzeugen.")
    p.add_argument("--toc-dir", default=str(ROOT / "build" / "dfm"), help="PyInstaller-Arbeitsordner (build/dfm)")
    p.add_argument("-o", "--output", default=str(ROOT / "dist" / APP_NAME / "THIRD-PARTY-NOTICES.txt"))
    args = p.parse_args(argv)
    toc_dir = Path(args.toc_dir)
    if not toc_dir.is_dir():
        print(f"Kein Build-Ordner: {toc_dir}", file=sys.stderr)
        return 2
    sys.path.insert(0, str(ROOT))
    from dicom_file_modifier import __version__

    names = names_from_toc_dir(toc_dir)
    text = write_notices(names, Path(args.output), __version__)
    print(f"{args.output}: {len(text) / 1024:.0f} KB, {len(bundled_distributions(names))} Pakete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
