"""``dfm``-Einstiegspunkt (Plan P0.8): Weiterleitung, Versionen, Self-Test."""

from __future__ import annotations

import json
import subprocess
import sys
import tomllib
from importlib.metadata import version
from pathlib import Path

import dicom_file_modifier
from dicom_file_modifier import analyzer, cli
from dicom_file_modifier.demo import make_demo_case

ROOT = Path(__file__).resolve().parents[1]


def test_version_is_the_same_everywhere(capsys):
    pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    assert dicom_file_modifier.__version__ == pyproject["project"]["version"] == version("dicom-file-modifier")
    assert cli.main(["--version", "--json"]) == 0
    info = json.loads(capsys.readouterr().out)
    assert info["dfm"] == dicom_file_modifier.__version__
    assert set(info["packages"]) == set(cli.PACKAGES) and all(info["packages"].values())


def test_tool_commands_forward_argv_unchanged(tmp_path, capsys):
    case = make_demo_case(tmp_path / "demo")
    assert cli.main(["analyze", str(case.rs), "--list"]) == 0
    via_dfm = capsys.readouterr().out
    assert analyzer.main([str(case.rs), "--list"]) == 0
    assert via_dfm == capsys.readouterr().out and "PTV_1" in via_dfm


def test_unknown_command_is_exit_2(capsys):
    assert cli.main(["gibtsnicht"]) == 2
    assert "Unbekannter Befehl" in capsys.readouterr().err


def test_module_entry_point():
    out = subprocess.run([sys.executable, "-m", "dicom_file_modifier", "--help"],
                         capture_output=True, text=True, timeout=120)
    assert out.returncode == 0 and "dose-indices" in out.stdout


def test_selftest_passes(capsys):
    assert cli.main(["selftest"]) == 0
    assert "Gesamt: PASS" in capsys.readouterr().out
