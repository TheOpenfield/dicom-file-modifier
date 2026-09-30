"""
Desktop-App (Plan Phase 1/2), offscreen auf dem Demo-Fall: Seite
Strukturanalyse Ende-zu-Ende ueber den echten Worker, Abbruch ueber den
``JobRunner``, und der GUI-Prozess laedt keine Plot-Module.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication  # noqa: E402

from dicom_file_modifier.api import dose, jobs, selection  # noqa: E402
from dicom_file_modifier.api.outputs import OutputSpec  # noqa: E402
from dicom_file_modifier.demo import DemoSpec, make_demo_case  # noqa: E402
from dicom_file_modifier.gui.config import AppConfig  # noqa: E402
from dicom_file_modifier.gui.jobs import JobRunner  # noqa: E402
from dicom_file_modifier.gui.window import MainWindow  # noqa: E402

FORBIDDEN = ("matplotlib.pyplot", "plotly", "dicom_file_modifier.visualizer", "dicom_file_modifier.dose_viz")


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture(scope="module")
def demo(tmp_path_factory):
    return make_demo_case(tmp_path_factory.mktemp("demo") / "case", DemoSpec())


def wait_until(app, pred, timeout=240.0):
    end = time.monotonic() + timeout
    while not pred():
        if time.monotonic() > end:
            raise AssertionError("Zeitueberschreitung")
        app.processEvents()
        time.sleep(0.02)


def test_structures_page_end_to_end(app, demo, tmp_path):
    win = MainWindow(AppConfig(results_root=tmp_path / "results", jobs_dir=tmp_path / "jobs"))
    win.show()
    win.open_case(demo.root)
    page = win.structures
    assert page.rs_combo.count() == 1
    wait_until(app, page.start_button.isEnabled)
    assert page.roi_table.rowCount() == len(page.info.rtstruct["rois"]) > 0
    assert "PTV_1" in page.preview_label.text()
    assert page.out_label.text().endswith(f"{demo.root.name}_STRUCT")

    page.start_button.click()
    assert not page.start_button.isEnabled() and win.cancel_button.isVisible()
    wait_until(app, lambda: page.last_result is not None)
    res = page.last_result
    assert res["status"] in ("ok", "ok_warnings"), res["issues"]
    out = Path(res["output_dir"])
    assert out == tmp_path / "results" / f"{demo.root.name}_STRUCT" and (out / "run.json").is_file()
    assert page.gallery.count() == sum(1 for p in out.glob("*.png")) > 0
    assert "statistics" in res["outputs"] and page.stats_view.toPlainText().strip()
    assert page.start_button.isEnabled() and not win.cancel_button.isVisible()
    assert page.out_label.text().endswith("_STRUCT_2")          # naechster Lauf in neuen Ordner
    assert "RTSTRUCT Analyse" in win.log.toPlainText()
    assert sorted(p.name for p in (tmp_path / "results").iterdir()) == [out.name]
    assert not any((tmp_path / "jobs").iterdir())                # Job-Datei aufgeraeumt
    win.close()


def test_dose_page_end_to_end(app, demo, tmp_path):
    win = MainWindow(AppConfig(results_root=tmp_path / "results", jobs_dir=tmp_path / "jobs"))
    win.show()
    win.open_case(demo.root)
    page = win.dose
    wait_until(app, page.start_button.isEnabled)
    text = page.preview_label.text()
    assert "PTV_1" in text and "20.00 Gy (aus dem RTPLAN)" in text and "Feingitter: 0.25 mm" in text
    assert "RD:" in page.files_label.text() and "Dmax" in page.dose_label.text()

    # Eclipse-kompatibel sperrt die Rasterfelder und zeigt die effektiven Werte
    mode, grid = page.form.widget("eclipse_compat"), page.form.widget("grid_mm")
    mode.setCurrentIndex(mode.findData("high"))
    assert not grid.isEnabled() and "Gesperrt" in grid.toolTip()
    assert "am CT-Pixelraster" in page.preview_label.text() and page.start_button.isEnabled()
    mode.setCurrentIndex(0)
    assert grid.isEnabled()

    # ungueltige Isodosen: Feld rot umrandet, Start gesperrt
    iso = page.form.widget("isodose")
    iso.setText("100,abc")
    assert not page.start_button.isEnabled() and "border" in iso.styleSheet()
    page.form.set_settings(dose.Settings(viz=False))
    assert page.start_button.isEnabled() and iso.styleSheet() == ""

    page.start_button.click()
    wait_until(app, lambda: page.last_result is not None)
    res = page.last_result
    assert res["status"] in ("ok", "ok_warnings"), res["issues"]
    assert page.metrics.rowCount() == len(res["summary"]["targets"]) == 1
    ci = res["summary"]["targets"]["PTV_1"]["ci_paddick"]
    assert page.metrics.item(0, 0).text() == "PTV_1" and page.metrics.item(0, 4).text() == f"{ci:.3f}"
    assert "PTV_1" in page.report.toPlainText() and not page.viz_button.isEnabled()
    assert (Path(res["output_dir"]) / "run.json").is_file()
    win.close()


def test_job_runner_cancel_leaves_nothing(app, demo, tmp_path):
    runner = JobRunner(tmp_path / "jobs")
    got = {}

    def on_event(ev):
        if ev["type"] == "stage" and ev["key"] == "compute":
            runner.cancel()

    runner.event.connect(on_event)
    runner.finished.connect(lambda r: got.setdefault("r", r))
    job = jobs.new_job("dose", dose.Settings(grid_mm=0.1, viz=False, write_rs=False),
                       selection.for_dose(str(demo.root)), OutputSpec(str(tmp_path / "results")))
    runner.start(job)
    wait_until(app, lambda: "r" in got)
    assert got["r"]["status"] == "cancelled" and not runner.running
    assert not job.staging_dir.exists() and not any((tmp_path / "results").iterdir())
    assert not any((tmp_path / "jobs").iterdir())


def test_gui_process_never_loads_plot_modules(demo, tmp_path):
    code = (
        "import json, os, sys, time\n"
        "os.environ['QT_QPA_PLATFORM'] = 'offscreen'\n"
        "from pathlib import Path\n"
        "from PySide6.QtWidgets import QApplication\n"
        "from dicom_file_modifier.gui.config import AppConfig\n"
        "from dicom_file_modifier.gui.window import MainWindow\n"
        "app = QApplication([])\n"
        "tmp = Path(sys.argv[2])\n"
        "win = MainWindow(AppConfig(results_root=tmp / 'r', jobs_dir=tmp / 'j'))\n"
        "win.open_case(sys.argv[1])\n"
        "end = time.monotonic() + 120\n"
        "ready = lambda: all(pg.start_button.isEnabled() for pg in win.pages)\n"
        "while not ready() and time.monotonic() < end:\n"
        "    app.processEvents(); time.sleep(0.02)\n"
        f"bad = sorted(m for m in sys.modules if m in {FORBIDDEN!r})\n"
        "print(json.dumps({'ready': ready(), 'bad': bad}))\n")
    out = subprocess.run([sys.executable, "-c", code, str(demo.root), str(tmp_path)],
                         capture_output=True, text=True, timeout=180, stdin=subprocess.DEVNULL)
    assert out.returncode == 0, out.stderr
    assert json.loads(out.stdout.strip().splitlines()[-1]) == {"ready": True, "bad": []}


def test_smoke_test_option_runs_the_structure_page(demo):
    """``dfm gui --smoke-test CASE``: die Installationspruefung der (gefrorenen) App."""
    env = {**os.environ, "QT_QPA_PLATFORM": "offscreen"}
    out = subprocess.run([sys.executable, "-m", "dicom_file_modifier", "gui", "--smoke-test", str(demo.root)],
                         capture_output=True, text=True, timeout=300, env=env, stdin=subprocess.DEVNULL)
    assert out.returncode == 0, out.stderr
