"""
analyzer / visualizer / dose_viz (Plan P0.7): stiller Kern, RNG je Lauf,
Issues statt verschluckter Fehler, inspect_rtstruct, Namenslisten, Plot-Bericht,
kein pyplot.
"""

from __future__ import annotations

import subprocess
import sys

import pytest

from dicom_file_modifier import analyzer as ana
from dicom_file_modifier import visualizer as viz
from dicom_file_modifier.demo import DemoSpec, make_demo_case


@pytest.fixture(scope="module")
def demo(tmp_path_factory):
    return make_demo_case(tmp_path_factory.mktemp("demo") / "case", DemoSpec())


def _jsonable(results):
    return ana._results_to_jsonable(results)


def test_core_is_silent_and_matches_run_analysis(demo, capsys):
    results, info = ana.analyze_rtstruct(str(demo.rs))
    assert capsys.readouterr().out == ""
    printed = ana.run_analysis(str(demo.rs))
    assert "ZIELGEBIETE" in capsys.readouterr().out
    assert _jsonable(results) == _jsonable(printed)
    assert {s for s, _ in info["analyzed"]} == {"targets", "oars", "helpers"}
    assert info["issues"] == []


def test_each_run_starts_the_random_numbers_afresh(demo, monkeypatch):
    # Kleine Obergrenze erzwingt das Ausduennen (sonst erst ab 50 000 Punkten)
    monkeypatch.setattr(ana.pair_distances, "__defaults__", (50,))
    first, info = ana.analyze_rtstruct(str(demo.rs))
    ana._cap_points(ana.np.zeros((200, 3)), 10)          # Zufallszahlen ausserhalb eines Laufs verbrauchen
    second, _ = ana.analyze_rtstruct(str(demo.rs))
    assert first["distances"] == second["distances"]
    assert "ANA.POINTS_CAPPED" in {i.code for i in info["issues"]}


def test_swallowed_geometry_errors_become_issues(demo, monkeypatch):
    def broken_hull(*args, **kwargs):
        raise RuntimeError("QHull kaputt")

    monkeypatch.setattr(ana, "ConvexHull", broken_hull)
    results, info = ana.analyze_rtstruct(str(demo.rs))
    hull = [i for i in info["issues"] if i.code == "ANA.HULL_FAILED"]
    assert hull and all(i.level == "warning" and "QHull kaputt" in i.detail for i in hull)
    assert any(i.message_de.startswith("PTV_1:") for i in hull)
    assert results["targets"]["PTV_1"]["shape"]["max_diameter_mm"] == 0.0


def test_inspect_rtstruct(demo, capsys):
    info = ana.inspect_rtstruct(str(demo.rs))
    assert capsys.readouterr().out == ""
    by_name = {r["name"]: r for r in info["rois"]}
    assert len(info["rois"]) == 14 and len(info["frame_of_reference_uids"]) == 1
    assert by_name["PTV_1"]["category"] == ana.CAT_TARGET and by_name["PTV_1"]["volume_cm3"] > 4
    assert by_name["HS1"]["is_marker"] and by_name["Leer"]["n_contours"] == 0
    assert {m["name"] for m in info["markers"]} == {"HS1", "Iso"}
    assert ana.inspect_rtstruct(str(demo.rs), volumes=False)["rois"][0]["volume_cm3"] is None


def test_name_lists_are_trimmed():
    assert ana.parse_name_list(" PTV , GTV,") == ["PTV", "GTV"]
    assert ana.parse_name_list("") is None and ana.parse_name_list(None) is None


def test_run_visualization_reports_written_and_skipped(demo, tmp_path, capsys):
    results, _ = ana.analyze_rtstruct(str(demo.rs))
    report = viz.run_visualization(results, tmp_path / "a")
    assert set(report["written"]) >= {"plot_volumes", "plot_distances", "statistics"}
    assert all((tmp_path / "a" / p).exists() for p in report["written"].values())
    empty = {"targets": {}, "oars": {}, "helpers": {}, "distances": [], "meta": {}}
    report = viz.run_visualization(empty, tmp_path / "b")
    assert "plot_volumes" in report["skipped"] and "(-) plot_volumes uebersprungen" in capsys.readouterr().out


def test_no_module_imports_pyplot():
    code = ("import sys, tempfile, pathlib\n"
            "import dicom_file_modifier.visualizer as v, dicom_file_modifier.dose_viz, "
            "dicom_file_modifier.dose_indices, dicom_file_modifier.case_modifier\n"
            "fig, ax = v._subplots(figsize=(2, 2)); ax.plot([0, 1], [0, 1])\n"
            "fig.savefig(pathlib.Path(tempfile.mkdtemp()) / 'x.png')\n"
            "print('matplotlib.pyplot' in sys.modules)")
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=120)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == "False"
