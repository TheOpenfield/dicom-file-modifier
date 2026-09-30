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

from PySide6.QtGui import QValidator  # noqa: E402
from PySide6.QtWidgets import QApplication  # noqa: E402

from dicom_file_modifier.api import dose, jobs, selection  # noqa: E402
from dicom_file_modifier.api.outputs import OutputSpec  # noqa: E402
from dicom_file_modifier.demo import DemoSpec, make_demo_case  # noqa: E402
from dicom_file_modifier.gui.config import AppConfig  # noqa: E402
from dicom_file_modifier.gui.jobs import JobRunner  # noqa: E402
from dicom_file_modifier.gui.widgets import DecimalSpinBox, SettingsForm  # noqa: E402
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
    assert win.windowTitle().endswith(demo.root.name) and page.state.text.text().startswith("Bereit")
    roles = {page.roi_table.item(i, 0).text(): page.roi_table.item(i, 4).text()
             for i in range(page.roi_table.rowCount())}
    assert roles["PTV_1"] == "Zielvolumen" and roles["Hirnstamm"] == "Risikoorgan"

    page.start_button.click()
    assert not page.start_button.isEnabled() and win.cancel_button.isVisible()
    wait_until(app, lambda: page.last_result is not None)
    res = page.last_result
    assert res["status"] in ("ok", "ok_warnings"), res["issues"]
    out = Path(res["output_dir"])
    assert out == tmp_path / "results" / f"{demo.root.name}_STRUCT" and (out / "run.json").is_file()
    assert page.gallery.count() == sum(1 for p in out.glob("*.png")) > 0
    assert "statistics" in res["outputs"] and page.stats_view.toPlainText().strip()
    assert "Zielvolumen" in page.summary_label.text()
    assert page.more_button.isEnabled() and page.runjson_action.isEnabled()      # im Menu "…"
    page.command_action.trigger()
    assert QApplication.clipboard().text().startswith("dfm analyze")
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
    assert "RD:  Plan-Summe" in page.files_label.text() and "Dmax" in page.files_label.text()
    assert page.dose_label.text().startswith("DVH: ")

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
    assert not page.start_button.isEnabled() and page.form.error_fields() == {"isodose"}
    assert page.state.text.text() == "Nicht startbar"
    page.form.set_settings(dose.Settings(viz=False))
    assert page.start_button.isEnabled() and not page.form.error_fields()

    page.start_button.click()
    wait_until(app, lambda: page.last_result is not None)
    res = page.last_result
    assert res["status"] in ("ok", "ok_warnings"), res["issues"]
    # Kennzahlen transponiert: je Ziel Wert, Eclipse, Abweichung; CI ohne Fehlalarm (ganze Isodose wie Eclipse)
    m = page.metrics
    assert [m.horizontalHeaderItem(j).text() for j in range(m.columnCount())] == ["Kennzahl", "PTV_1", "Eclipse",
                                                                                   "Abw. %"]
    rows = {m.item(i, 0).text(): i for i in range(m.rowCount())}
    ci = res["summary"]["targets"]["PTV_1"]["ci_paddick"]
    assert m.item(rows["CI Paddick"], 1).text() == f"{ci:.3f}" and not m.item(rows["CI Paddick"], 2).text()
    g = rows["CI Paddick global"]                                 # PIV-Bereich component: eigene Vergleichszeile
    assert m.item(g, 2).text() and m.item(g, 3).icon().isNull() and "TV global [cm³]" not in rows
    assert "alle innerhalb der Toleranz" in page.metrics_note.text()
    assert not any(i["code"] == "DOSE.ECLIPSE_TOLERANCE" for i in res["issues"])
    assert "PTV_1" in page.report.toPlainText() and not page.viz_button.isEnabled()
    assert (Path(res["output_dir"]) / "run.json").is_file()
    win.close()


def test_transform_page_end_to_end(app, demo, tmp_path):
    win = MainWindow(AppConfig(results_root=tmp_path / "results", jobs_dir=tmp_path / "jobs"))
    win.show()
    win.open_case(demo.root)
    page = win.transform
    win.nav.setCurrentRow(win.pages.index(page))
    wait_until(app, page.start_button.isEnabled)
    assert "Volumenmitte" in page.case_label.text() and "keine Bewegung" in page.preview_label.text()

    # Drehpunkt: Marker aus dem RTSTRUCT; "Koordinate" zeigt das Eingabefeld mit der Volumenmitte
    combo, marker = page.center_combo, page.center_combo.findData("marker:HS1")
    assert marker > 0 and not page.center_edit.isVisible()
    combo.setCurrentIndex(combo.count() - 1)
    assert page.center_edit.isVisible() and page.form.settings().center.count(",") == 2
    combo.setCurrentIndex(marker)
    page.center_edit.setText("marker:hs1")                       # Kleinschreibung: derselbe Marker
    assert combo.currentIndex() == marker
    page.form.widget("tx").setValue(2.0)
    page.form.widget("rz").setValue(5.0)
    text = page.preview_label.text()
    assert "2 mm nach links · +5° um die Kopf-Fuss-Achse" in text and "Drehpunkt: Marker HS1 (" in text
    assert "FrameOfReference: beibehalten" in text

    page.start_button.click()
    wait_until(app, lambda: page.last_result is not None)
    res = page.last_result
    assert res["status"] in ("ok", "ok_warnings"), res["issues"]
    out = Path(res["output_dir"])
    assert out.name == f"{demo.root.name}_RB" and (out / "CT").is_dir() and (out / "RS_RB.dcm").is_file()
    assert page.overview.count() == page.displacement.count() == 1 and page.view3d_button.isEnabled()
    report = page.report.toPlainText()
    assert "Matrix T" in report and "Drehpunkt     Marker HS1" in report and "FoR           beibehalten" in report
    assert "2 mm nach links" in page.summary_label.text() and "Verschiebung X" not in page.summary_label.text()
    page.command_action.trigger()
    assert QApplication.clipboard().text().startswith("dfm case-transform")

    # nur CT: Drehpunkt fest auf der Volumenmitte, FoR gesperrt, Kennung nur fuer den Ordner
    page.ct_only.setChecked(True)
    wait_until(app, page.start_button.isEnabled)
    assert page.form.settings().center == "volume" and not page.center_combo.isEnabled()
    assert not page.form.widget("new_frame_of_reference").isEnabled() and page.form.widget("label").isEnabled()
    assert "nur das CT" in page.preview_label.text()
    page.ct_only.setChecked(False)                                # die Wahl von vorher kommt zurueck
    wait_until(app, page.start_button.isEnabled)
    assert page.form.settings().center == "marker:hs1"
    win.close()


def test_decimal_spinbox_accepts_a_comma_and_shows_a_point(app):
    w = DecimalSpinBox()
    w.setDecimals(2)
    w.setRange(-100.0, 100.0)
    assert w.textFromValue(2.5) == "2.50" and w.valueFromText("2,5") == 2.5
    assert w.validate("-2,5", 4)[0] == QValidator.State.Acceptable


def test_cli_options_in_core_texts_become_field_names(app):
    form = SettingsForm(dose.Settings)
    assert form.gui_text("weitere PTVs nicht ausgewertet (--target fuer alle).") == \
        "weitere PTVs nicht ausgewertet („Zielvolumen“ fuer alle)."
    assert form.gui_text("mit --no-rs") == f"mit „{form.label_of('write_rs')}“ aus"
    assert form.gui_text("mit --list anzeigen") == "mit --list anzeigen"      # keine Einstellung


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
