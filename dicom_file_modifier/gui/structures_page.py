"""
Seite Strukturanalyse: RTSTRUCT waehlen -> ROI-Tabelle -> Einstellungen ->
Pruefung (``structures.preview``) -> Start im Worker -> Plots und statistics.txt.
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QFontDatabase
from PySide6.QtWidgets import (QAbstractItemView, QComboBox, QGroupBox, QHBoxLayout, QLabel,
                               QPlainTextEdit, QPushButton, QSplitter, QTableWidget,
                               QTableWidgetItem, QVBoxLayout, QWidget)

from ..api import jobs as api_jobs
from ..api import selection, structures
from ..api.issues import Issue, issue_from_exception
from ..api.outputs import OutputSpec
from .jobs import run_in_background
from .widgets import STATUS_DE, Gallery, IssueList, SettingsForm, breakable, open_path

CATEGORY_DE = {"TARGET": "Zielvolumen", "OAR_SERIAL": "Risikoorgan (seriell)",
               "OAR_PARALLEL": "Risikoorgan (parallel)", "HELPER": "Hilfsstruktur",
               "EXTERNAL": "Außenkontur", "MARKER": "Marker"}


def _names(names: list, limit: int = 8) -> str:
    if not names:
        return "keine"
    more = f" … (+{len(names) - limit})" if len(names) > limit else ""
    return ", ".join(names[:limit]) + more


class StructuresPage(QWidget):
    title = "Strukturanalyse"

    def __init__(self, main):
        super().__init__()
        self.main = main                                    # MainWindow
        self.case = None
        self.info = None
        self.last_result = None
        self._token = 0
        self._running = False
        self._ok = False

        # -- links: Eingaben, Einstellungen, Pruefung, Start
        self.rs_combo = QComboBox()
        self.rs_combo.currentIndexChanged.connect(self._inspect)
        self.rs_label = QLabel()
        self.roi_table = QTableWidget(0, 4)
        self.roi_table.setHorizontalHeaderLabels(["Name", "DICOM-Typ", "Kategorie", "Volumen [cm³]"])
        self.roi_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.roi_table.verticalHeader().hide()
        self.roi_table.horizontalHeader().setStretchLastSection(True)
        inputs = QGroupBox("RTSTRUCT")
        box = QVBoxLayout(inputs)
        box.addWidget(self.rs_combo)
        box.addWidget(self.rs_label)
        box.addWidget(self.roi_table, 1)

        self.form = SettingsForm(structures.Settings)
        self.form.changed.connect(self._update_preview)
        settings_box = QGroupBox("Einstellungen")
        QVBoxLayout(settings_box).addWidget(self.form)

        self.preview_label = QLabel("Kein Fall geöffnet.")
        self.preview_label.setWordWrap(True)
        self.out_label = QLabel()
        self.out_label.setWordWrap(True)
        self.issues = IssueList()
        self.issues.setMaximumHeight(110)
        self.issues.hide()
        check = QGroupBox("Prüfung")
        box = QVBoxLayout(check)
        box.addWidget(self.preview_label)
        box.addWidget(self.issues)
        box.addWidget(self.out_label)

        self.start_button = QPushButton("Analyse starten")
        self.start_button.setEnabled(False)
        self.start_button.clicked.connect(self._start)

        left = QWidget()
        col = QVBoxLayout(left)
        col.addWidget(inputs, 1)
        col.addWidget(settings_box)
        col.addWidget(check)
        col.addWidget(self.start_button)

        # -- rechts: Ergebnis
        self.status_label = QLabel("Noch kein Ergebnis.")
        self.status_label.setWordWrap(True)
        self.open_button = QPushButton("Ordner öffnen")
        self.open_button.setEnabled(False)
        self.open_button.clicked.connect(lambda: open_path(self.last_result["output_dir"]))
        self.result_issues = IssueList()
        self.result_issues.setMaximumHeight(110)
        self.result_issues.hide()
        self.gallery = Gallery()
        self.stats_view = QPlainTextEdit()
        self.stats_view.setReadOnly(True)
        self.stats_view.setFont(QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont))
        result = QGroupBox("Ergebnis")
        box = QVBoxLayout(result)
        head = QHBoxLayout()
        head.addWidget(self.status_label, 1)
        head.addWidget(self.open_button)
        box.addLayout(head)
        box.addWidget(self.result_issues)
        split = QSplitter(Qt.Orientation.Vertical)
        split.addWidget(self.gallery)
        split.addWidget(self.stats_view)
        box.addWidget(split, 1)

        outer = QSplitter()
        outer.addWidget(left)
        outer.addWidget(result)
        outer.setStretchFactor(1, 1)
        outer.setSizes([560, 760])
        QVBoxLayout(self).addWidget(outer)

    # -- Fall und Inspektion ---------------------------------------------------------
    def set_case(self, case) -> None:
        self.case = case
        self.rs_combo.blockSignals(True)
        self.rs_combo.clear()
        for p in case.rs:
            self.rs_combo.addItem(p.name, str(p))
        self.rs_combo.blockSignals(False)
        self._inspect()

    def _inspect(self) -> None:
        self._token += 1
        token, self.info = self._token, None
        self.roi_table.setRowCount(0)
        self.rs_label.clear()
        rs = self.rs_combo.currentData()
        if not rs:
            self._show_check("Kein RTSTRUCT im Fallordner (RS*.dcm).", [], False)
            return
        self._show_check("RTSTRUCT wird gelesen …", [], False)
        case_id = self.case.folder.name
        run_in_background(lambda: structures.inspect(selection.from_rtstruct(rs, case_id=case_id)),
                          lambda info, err: self._on_inspected(token, info, err))

    def _on_inspected(self, token: int, info, err) -> None:
        if token != self._token:
            return                                          # veraltet: inzwischen anderes RS
        if err is not None:
            self._show_check("RTSTRUCT nicht lesbar.", [issue_from_exception(err)], False)
            return
        self.info = info
        rt = info.rtstruct
        if rt:
            self.rs_label.setText(f"{rt.get('structure_set_label') or '-'}   ·   {len(rt['rois'])} ROIs"
                                  f"   ·   {len(rt.get('markers', []))} Marker")
            self.roi_table.setRowCount(len(rt["rois"]))
            for i, r in enumerate(rt["rois"]):
                vol = "" if r.get("volume_cm3") is None else f"{r['volume_cm3']:.2f}"
                for j, text in enumerate((r["name"], r.get("rt_type", ""),
                                          CATEGORY_DE.get(r["category"], r["category"]), vol)):
                    self.roi_table.setItem(i, j, QTableWidgetItem(text))
            self.roi_table.resizeColumnsToContents()
        self.form.set_placeholder("targets", "automatisch: " + _names(info.auto_targets))
        self.form.set_placeholder("oars", "automatisch: " + _names(info.auto_oars))
        self._update_preview()

    # -- Pruefung --------------------------------------------------------------------
    def _output_spec(self) -> OutputSpec:
        return OutputSpec(str(self.main.config.results_root))

    def _update_preview(self) -> None:
        if self.info is None:
            return
        try:
            s = self.form.settings()
        except ValueError as exc:
            self._show_check("Eingabe prüfen.", [Issue("error", "INPUT.INVALID", str(exc))], False)
            return
        pv = structures.preview(self.info, s)
        text = (f"Zielvolumen: {_names(pv['targets'])}\nRisikoorgane: {_names(pv['oars'])}\n"
                f"Hilfsstrukturen: {len(pv['helpers'])}")
        out = self._output_spec()
        issues = pv["issues"] + out.validate()
        folder = out.target_dir(structures.default_folder(s, self.info.selection))
        self.out_label.setText(f"Ergebnisordner: {breakable(folder)}")
        self._show_check(text, issues, pv["ok"] and not any(i.level == "error" for i in issues))

    def _show_check(self, text: str, issues: list, ok: bool) -> None:
        self.preview_label.setText(text)
        self.issues.set_issues(issues)
        if not ok:
            self.out_label.clear()
        self._ok = ok
        self._update_start()

    def _update_start(self) -> None:
        self.start_button.setEnabled(self._ok and not self._running)
        self.start_button.setToolTip("" if self._ok else "Erst einen Fall mit RTSTRUCT öffnen "
                                                          "und die Prüfung bestehen.")

    def refresh(self) -> None:
        self._update_preview()

    def set_running(self, running: bool) -> None:
        self._running = running
        self._update_start()

    # -- Lauf und Ergebnis -----------------------------------------------------------
    def _start(self) -> None:
        job = api_jobs.new_job(structures.WORKFLOW, self.form.settings(), self.info.selection,
                               self._output_spec())
        self.last_result = None
        self.status_label.setText("Läuft …")
        self.open_button.setEnabled(False)
        self.result_issues.set_issues([])
        self.gallery.clear()
        self.stats_view.clear()
        self.main.start_job(self, job, self.title)

    def show_result(self, result: dict) -> None:
        self.last_result = result
        out = result.get("output_dir")
        self.status_label.setText(STATUS_DE.get(result["status"], result["status"])
                                  + (f":  {breakable(out)}" if out else ""))
        self.open_button.setEnabled(bool(out))
        self.result_issues.set_issues(result.get("issues", []))
        if out:
            outputs = result.get("outputs", {})
            self.gallery.set_images([Path(out) / rel for role, rel in outputs.items()
                                     if role.startswith("plot:") and str(rel).lower().endswith(".png")])
            stats = outputs.get("statistics")
            if stats and (Path(out) / stats).is_file():
                self.stats_view.setPlainText((Path(out) / stats).read_text(encoding="utf-8",
                                                                            errors="replace"))
        self._update_preview()                              # naechster Lauf bekaeme _2
