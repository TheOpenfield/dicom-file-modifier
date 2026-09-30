"""
cli.py - Einstiegspunkt ``dfm`` fuer alle Werkzeuge des Pakets (Plan P0.8).

  dfm analyze RS.dcm [...]          wie python -m dicom_file_modifier.analyzer
  dfm visualize RS.dcm [...]        ... visualizer
  dfm ct-transform CT_DIR [...]     ... modifier
  dfm case-transform CASE [...]     ... case_modifier
  dfm dose-indices CASE [...]       ... dose_indices
  dfm demo OUT [...]                ... demo (synthetischer Demo-Fall)
  dfm selftest                      alle Self-Tests; case_modifier auf einem frischen Demo-Fall
  dfm worker --job JOB.json         einen Job der Desktop-App ausfuehren (JSON-Lines)
  dfm --version [--json]

Die Werkzeug-Befehle leiten ihre Argumente unveraendert an ``main(argv)`` des
Moduls weiter: Optionen, Ausgaben und Exit-Codes sind dieselben.
"""

from __future__ import annotations

import importlib
import json
import sys
import tempfile
from pathlib import Path
from typing import Optional

from . import __version__

# Befehl -> (Modul, Kurzbeschreibung)
TOOLS = {
    "analyze": ("analyzer", "Geometrische Analyse eines RTSTRUCT (JSON + Konsole)"),
    "visualize": ("visualizer", "Plots und statistics.txt aus einem RTSTRUCT"),
    "ct-transform": ("modifier", "Starre Transformation einer CT-Serie"),
    "case-transform": ("case_modifier", "CT + RTSTRUCT gemeinsam transformieren"),
    "dose-indices": ("dose_indices", "Dosisindizes (CI, GI, HI) und Isodosen-RTSTRUCT"),
    "demo": ("demo", "Synthetischen Demo-Fall erzeugen (keine Patientendaten)"),
}
PACKAGES = ("numpy", "scipy", "pydicom", "shapely", "matplotlib", "scikit-image", "plotly")


def _usage() -> str:
    lines = [f"dfm {__version__} - DICOM-RT-Werkzeuge", "", "Befehle:"]
    lines += [f"  {cmd:<15} {desc}" for cmd, (_, desc) in TOOLS.items()]
    lines += [f"  {'selftest':<15} Alle Self-Tests (Installationspruefung)",
              f"  {'worker':<15} Job der Desktop-App ausfuehren (--job JOB.json)",
              "", "  dfm <befehl> --help    Optionen des Befehls",
              "  dfm --version [--json] Versionen"]
    return "\n".join(lines)


def versions() -> dict:
    """Paket-, Werkzeug- und Bibliotheksversionen (fuer --version, run.json, Diagnose)."""
    from importlib.metadata import PackageNotFoundError, version

    from .dose_constants import TOOL_NAME, TOOL_VERSION

    out = {"dfm": __version__, "python": sys.version.split()[0],
           "tools": {TOOL_NAME: TOOL_VERSION}, "packages": {}}
    for pkg in PACKAGES:
        try:
            out["packages"][pkg] = version(pkg)
        except PackageNotFoundError:
            out["packages"][pkg] = None
    return out


def _console_setup() -> None:
    """Nicht darstellbare Zeichen als '?' statt UnicodeEncodeError (cp1252-Konsole)."""
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(errors="replace")
        except (AttributeError, ValueError):     # None (pythonw) oder kein Textstrom
            pass


def selftest(argv: Optional[list] = None) -> int:
    """
    Installationspruefung: Self-Tests von analyzer, dose_indices und
    rtstruct_writer sowie der case_modifier-Self-Test auf einem frisch
    erzeugten Demo-Fall.  Exit 0 = alles bestanden.
    """
    import argparse

    p = argparse.ArgumentParser(prog="dfm selftest", description=selftest.__doc__)
    p.parse_args(argv)
    from . import analyzer, case_modifier, demo, dose_indices, rtstruct_writer

    results = []

    def run(name, fn):
        print(f"\n===== Self-Test: {name} =====")
        try:
            code = int(fn() or 0)
        except SystemExit as e:                     # argparse & Co.
            code = int(e.code or 0)
        except Exception as e:  # noqa: BLE001 - jeder Fehler ist ein FAIL der Pruefung
            print(f"  Abbruch: {type(e).__name__}: {e}")
            code = 1
        results.append((name, code))

    run("analyzer", lambda: analyzer.main(["--self-test"]))
    run("dose_indices", lambda: dose_indices.main(["--self-test"]))
    run("rtstruct_writer", lambda: rtstruct_writer.main(["--self-test"]))
    with tempfile.TemporaryDirectory(prefix="dfm_selftest_") as tmp:
        case = demo.make_demo_case(Path(tmp) / "demo")
        run("case_modifier (Demo-Fall)", lambda: case_modifier.main([str(case.root), "--self-test"]))

    print("\n===== Ergebnis =====")
    for name, code in results:
        print(f"  {'PASS' if code == 0 else 'FAIL'}  {name}")
    ok = all(code == 0 for _, code in results)
    print(f"\n  Gesamt: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def main(argv: Optional[list] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    _console_setup()
    if not argv or argv[0] in ("-h", "--help", "help"):
        print(_usage())
        return 0
    cmd, rest = argv[0], argv[1:]
    if cmd == "--version":
        if rest[:1] == ["--json"]:
            print(json.dumps(versions(), indent=2))
        else:
            print(f"dfm {__version__}")
        return 0
    if cmd in TOOLS:
        module = importlib.import_module(f"dicom_file_modifier.{TOOLS[cmd][0]}")
        return int(module.main(rest) or 0)
    if cmd == "selftest":
        return selftest(rest)
    if cmd == "worker":
        from .api import worker
        return worker.main(rest)
    print(f"Unbekannter Befehl {cmd!r}.\n\n{_usage()}", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
