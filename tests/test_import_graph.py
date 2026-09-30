"""
Importgraph des Pakets (Plan P0.3, Schichten ab P0.8).

- Auf Modulebene gibt es keine Import-Zyklen.
- Kein Paketmodul importiert die CLI-Orchestratoren ``case_modifier`` oder
  ``dose_indices``, weder auf Modulebene noch in Funktionen.
- Schichten: Kein Kernmodul importiert ``api`` (ausser dem Einstiegspunkt
  ``cli``); ``api`` importiert nie ``gui``.
- Jedes Modul (auch in ``api/``) laesst sich in einem frischen Interpreter als
  erstes importieren.
- Die nach ``dicom_utils``/``dose_constants`` verschobenen Namen sind unter den
  alten Pfaden dieselben Objekte.
"""

from __future__ import annotations

import ast
import importlib
import subprocess
import sys
from pathlib import Path

import pytest

PKG_DIR = Path(__file__).resolve().parents[1] / "dicom_file_modifier"
PKG = PKG_DIR.name
MODULES = sorted(p.stem for p in PKG_DIR.glob("*.py") if p.stem not in ("__init__", "__main__"))
ORCHESTRATORS = {"case_modifier", "dose_indices"}
# Einstiegspunkte duerfen die Orchestratoren aufrufen (Unterpakete wie api/ werden nicht gescannt)
ENTRY_POINTS = {"cli"}


def _package_imports(node: ast.AST) -> set:
    """Namen der Paketmodule, die ein Import-Knoten laedt."""
    if isinstance(node, ast.ImportFrom):
        if node.level == 1:
            return {node.module.split(".")[0]} if node.module else {a.name for a in node.names}
        if node.module == PKG:
            return {a.name for a in node.names}
        if node.module and node.module.startswith(PKG + "."):
            return {node.module.split(".")[1]}
    elif isinstance(node, ast.Import):
        return {a.name.split(".")[1] for a in node.names if a.name.startswith(PKG + ".")}
    return set()


def _imports(module: str) -> tuple:
    """(beim Import ausgefuehrte, in Funktionen verzoegerte) Paketimporte eines Moduls."""
    tree = ast.parse((PKG_DIR / f"{module}.py").read_text(encoding="utf-8"))
    eager, lazy = set(), set()

    def visit(node: ast.AST, in_function: bool) -> None:
        for child in ast.iter_child_nodes(node):
            inner = in_function or isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda))
            (lazy if inner else eager).update(_package_imports(child))
            visit(child, inner)

    visit(tree, False)
    return eager & set(MODULES), lazy & set(MODULES)


GRAPH = {m: _imports(m) for m in MODULES}


def test_no_cycles_at_module_level():
    state: dict = {}
    cycles = []

    def dfs(m: str, path: list) -> None:
        state[m] = "open"
        for dep in sorted(GRAPH[m][0]):
            if state.get(dep) == "open":
                cycles.append(" -> ".join(path[path.index(dep):] + [dep]))
            elif dep not in state:
                dfs(dep, path + [dep])
        state[m] = "done"

    for m in MODULES:
        if m not in state:
            dfs(m, [m])
    assert not cycles, "Import-Zyklen auf Modulebene: " + "; ".join(cycles)


def test_no_module_imports_an_orchestrator():
    offenders = {m: sorted((eager | lazy) & ORCHESTRATORS)
                 for m, (eager, lazy) in GRAPH.items()
                 if m not in ORCHESTRATORS | ENTRY_POINTS and (eager | lazy) & ORCHESTRATORS}
    assert not offenders, f"Module importieren CLI-Orchestratoren: {offenders}"


API_DIR = PKG_DIR / "api"
API_MODULES = sorted(p.stem for p in API_DIR.glob("*.py") if p.stem != "__init__")


def _top_level_imports(path: Path, depth: int, own: str = "") -> set:
    """
    Direkte Kinder des Pakets (Module/Unterpakete), die eine Datei importiert,
    auf Modulebene und in Funktionen.  ``depth`` = Tiefe der Datei unter dem
    Paket (0 = Paketmodul, 1 = ``api/``); Importe innerhalb des eigenen
    Unterpakets zaehlen als ``own``.
    """
    out = set()
    for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
        if isinstance(node, ast.ImportFrom):
            if node.level == 0:
                if node.module == PKG:
                    out |= {a.name for a in node.names}
                elif node.module and node.module.startswith(PKG + "."):
                    out.add(node.module.split(".")[1])
            elif node.level - 1 == depth:                     # relativ zur Paketwurzel
                out |= {node.module.split(".")[0]} if node.module else {a.name for a in node.names}
            else:                                             # im eigenen Unterpaket
                out.add(own)
        elif isinstance(node, ast.Import):
            out |= {a.name.split(".")[1] for a in node.names if a.name.startswith(PKG + ".")}
    return out


def test_core_does_not_import_the_api_layer():
    offenders = sorted(m for m in MODULES if m not in ENTRY_POINTS
                       and {"api", "gui"} & _top_level_imports(PKG_DIR / f"{m}.py", 0))
    assert not offenders, f"Kernmodule importieren api/gui: {offenders}"
    assert "api" in _top_level_imports(PKG_DIR / "cli.py", 0)      # die Pruefung greift


def test_api_never_imports_the_gui():
    found = {m: _top_level_imports(API_DIR / f"{m}.py", 1, own="api")
             for m in API_MODULES + ["__init__"]}
    offenders = sorted(m for m, names in found.items() if "gui" in names)
    assert not offenders, f"api-Module importieren gui: {offenders}"
    assert "dose_indices" in found["dose"] and "api" in found["dose"]   # die Pruefung greift


@pytest.mark.parametrize("module", MODULES + [f"api.{m}" for m in API_MODULES])
def test_module_imports_first_in_fresh_interpreter(module):
    proc = subprocess.run([sys.executable, "-c", f"import {PKG}.{module}"],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=120)
    assert proc.returncode == 0, proc.stderr.decode("utf-8", "replace")[-2000:]


@pytest.mark.parametrize("old, new", [
    ("case_modifier.get_rs_frame_of_references", "dicom_utils.get_rs_frame_of_references"),
    ("case_modifier._label_with_suffix", "dicom_utils._label_with_suffix"),
    ("case_modifier._truncate", "dicom_utils._truncate"),
    ("modifier.set_sop_instance_uid", "dicom_utils.set_sop_instance_uid"),
    ("case_modifier.find_point_markers", "dicom_utils.find_point_markers"),
    ("case_modifier.validate_ct_geometry", "modifier.validate_ct_geometry"),
    ("case_modifier.load_ct_headers", "modifier.load_ct_headers"),
    ("dose_indices.TOOL_NAME", "dose_constants.TOOL_NAME"),
    ("dose_indices.TOOL_VERSION", "dose_constants.TOOL_VERSION"),
    ("dose_indices.LEVEL_COLORS", "dose_constants.LEVEL_COLORS"),
    ("dose_indices.HELPER_COLORS", "dose_constants.HELPER_COLORS"),
])
def test_old_import_paths_still_work(old, new):
    def resolve(dotted: str):
        mod, attr = dotted.split(".")
        return getattr(importlib.import_module(f"{PKG}.{mod}"), attr)

    assert resolve(old) is resolve(new)
