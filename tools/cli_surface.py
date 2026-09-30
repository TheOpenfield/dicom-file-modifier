#!/usr/bin/env python
"""
cli_surface.py - Schnappschuss der argparse-Oberflaeche aller CLIs (Plan P0.1).

Faengt in jedem Modul den Parser ab, den ``main()`` baut (``parse_args`` wird
durch eine Ausnahme ersetzt), und schreibt je Option ``option_strings``,
``dest``, ``default``, ``choices``, ``type``, ``nargs``, ``required`` und
die Action-Klasse nach ``<ausgabe>/cli_surface.json``.  Hilfetexte bleiben
bewusst aussen vor (Formulierungen duerfen sich aendern, Flags nicht).

  python tools/cli_surface.py <ausgabeordner>
"""

from __future__ import annotations

import argparse
import importlib
import inspect
import json
import sys
from pathlib import Path

MODULES = ("analyzer", "visualizer", "modifier", "case_modifier", "dose_indices",
           "rtstruct_writer", "dose", "demo")


class _Captured(Exception):
    pass


def _fake_parse_args(self, args=None, namespace=None):
    raise _Captured(self)


def _describe(parser: argparse.ArgumentParser) -> list:
    out = []
    for a in parser._actions:
        typ = a.type
        out.append({
            "option_strings": list(a.option_strings),
            "dest": a.dest,
            "default": repr(a.default),
            "choices": [repr(c) for c in a.choices] if a.choices is not None else None,
            "type": getattr(typ, "__name__", repr(typ)) if typ is not None else None,
            "nargs": a.nargs if not isinstance(a.nargs, int) else int(a.nargs),
            "required": bool(a.required),
            "action": type(a).__name__.lstrip("_"),
        })
    return out


def main(argv=None) -> int:
    args = sys.argv[1:] if argv is None else argv
    out_dir = Path(args[0]) if args else Path(".")
    out_dir.mkdir(parents=True, exist_ok=True)
    argparse.ArgumentParser.parse_args = _fake_parse_args
    result = {}
    for name in MODULES:
        mod = importlib.import_module("dicom_file_modifier." + name)
        fn = getattr(mod, "main", None)
        if fn is None:
            result[name] = {"error": "kein main()"}
            continue
        try:
            if inspect.signature(fn).parameters:
                fn([])
            else:
                fn()
            result[name] = {"error": "main() hat nicht geparst"}
        except _Captured as c:
            result[name] = {"options": _describe(c.args[0])}
        except SystemExit as e:  # pragma: no cover - sollte nicht vorkommen
            result[name] = {"error": f"SystemExit {e.code}"}
    (out_dir / "cli_surface.json").write_text(json.dumps(result, indent=1), encoding="utf-8")
    print(f"CLI-Oberflaeche von {len(result)} Modulen geschrieben: {out_dir / 'cli_surface.json'}")
    return 0 if all("options" in v for v in result.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
