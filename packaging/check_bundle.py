"""
Pruefung des onedir-Bundles (Standard ``dist/DICOM-RT-Toolkit``); Exit 0 = ok.

  - Groesse unter MAX_MB, keine DICOM-Dateien im Bundle
  - dfm.exe: Versionen (alle Bibliotheken gefunden), alle Self-Tests,
    Demo-Fall, Dosisindizes mit Validierungsansicht (plotly, matplotlib)
  - DICOM-RT-Toolkit.exe --smoke-test: Strukturanalyse des Demo-Falls ueber
    den gefrorenen Worker (offscreen)

  python packaging/check_bundle.py [BUNDLE_DIR] [--zip]

Mit ``--zip`` wird ein bestandenes Bundle als Portable-ZIP neben den Ordner
gelegt (``DICOM-RT-Toolkit-<version>-win64.zip``).
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_BUNDLE = ROOT / "dist" / "DICOM-RT-Toolkit"
MAX_MB = 450


def dicom_files(bundle: Path) -> list:
    """``.dcm``-Dateien und Dateien mit DICOM-Praeambel (``DICM`` an Byte 128)."""
    found = []
    for p in bundle.rglob("*"):
        if not p.is_file():
            continue
        if p.suffix.lower() == ".dcm":
            found.append(p)
            continue
        with open(p, "rb") as f:
            if f.read(132)[128:132] == b"DICM":
                found.append(p)
    return found


def run(step: str, cmd: list, env=None, timeout: float = 900) -> subprocess.CompletedProcess:
    t0 = time.monotonic()
    proc = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, encoding="utf-8",
                          errors="replace", timeout=timeout, env=env, stdin=subprocess.DEVNULL)
    print(f"{'ok    ' if proc.returncode == 0 else 'FEHLER'} {step} ({time.monotonic() - t0:.1f} s)")
    if proc.returncode:
        print(proc.stdout[-2000:], proc.stderr[-4000:], sep="\n")
    return proc


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    sys.stdout.reconfigure(errors="replace")         # Ausgaben der EXEs auf der cp1252-Konsole
    make_zip = "--zip" in argv
    rest = [a for a in argv if a != "--zip"]
    bundle = Path(rest[0]) if rest else DEFAULT_BUNDLE
    version = "?"
    dfm, gui = bundle / "dfm.exe", bundle / "DICOM-RT-Toolkit.exe"
    missing_exe = [p.name for p in (dfm, gui) if not p.is_file()]
    if missing_exe:
        print(f"FEHLER Bundle unvollstaendig, es fehlt: {', '.join(missing_exe)}")
        return 1
    errors = []

    size_mb = sum(p.stat().st_size for p in bundle.rglob("*") if p.is_file()) / 2**20
    print(f"Groesse {size_mb:.0f} MB")
    if size_mb > MAX_MB:
        errors.append(f"Bundle groesser als {MAX_MB} MB")
    bad = dicom_files(bundle)
    if bad:
        errors.append("DICOM-Dateien im Bundle: " + ", ".join(str(p.relative_to(bundle)) for p in bad[:10]))

    tmp = Path(tempfile.mkdtemp(prefix="dfm_b_"))
    try:
        proc = run("dfm --version --json", [dfm, "--version", "--json"])
        if proc.returncode == 0:
            info = json.loads(proc.stdout)
            version = info["dfm"]
            missing = [k for k, v in info["packages"].items() if not v]
            if missing:
                errors.append(f"Versionen fehlen: {missing}")
        steps = [
            ("dfm selftest", [dfm, "selftest"]),
            ("dfm demo", [dfm, "demo", tmp / "demo"]),
            ("dfm dose-indices", [dfm, "dose-indices", tmp / "demo", "--output", tmp / "out"]),
        ]
        for step, cmd in steps:
            if run(step, cmd).returncode:
                errors.append(step)
        viz = sorted(p.name for p in (tmp / "out").rglob("*") if p.suffix in (".html", ".png"))
        if viz != ["dose_overview.png", "validation.html"]:
            errors.append(f"Validierungsansicht unvollstaendig: {viz}")
        env = {**os.environ, "QT_QPA_PLATFORM": "offscreen"}
        if run("App --smoke-test (Strukturanalyse ueber den Worker)",
               [gui, "--smoke-test", tmp / "demo"], env=env).returncode:
            errors.append("App-Smoke-Test")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    for e in errors:
        print("FEHLER", e)
    print("Bundle ok" if not errors else f"{len(errors)} Fehler")
    if errors:
        return 1
    if make_zip:
        base = bundle.parent / f"{bundle.name}-{version}-win64"
        archive = shutil.make_archive(str(base), "zip", bundle.parent, bundle.name)
        print(f"ZIP {archive} ({Path(archive).stat().st_size / 2**20:.0f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
