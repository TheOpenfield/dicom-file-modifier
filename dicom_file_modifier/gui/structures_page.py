"""
Seite Strukturanalyse: RTSTRUCT waehlen -> ROI-Tabelle -> Einstellungen ->
Pruefung (``structures.preview``) -> Start im Worker -> Plots und statistics.txt.
"""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QComboBox, QGroupBox, QLabel, QSplitter, QVBoxLayout

from ..api import selection, structures
from .page import WorkflowPage, names_text
from .widgets import Gallery, SettingsForm, fill_table, make_table, report_view

CATEGORY_DE = {"TARGET": "Zielvolumen", "OAR_SERIAL": "Risikoorgan (seriell)",
               "OAR_PARALLEL": "Risikoorgan (parallel)", "HELPER": "Hilfsstruktur",
               "EXTERNAL": "Außenkontur", "MARKER": "Marker"}


class StructuresPage(WorkflowPage):
    title = "Strukturanalyse"
    start_text = "Analyse starten"
    api = structures

    def __init__(self, main):
        super().__init__(main)
        self.rs_combo = QComboBox()
        self.rs_combo.currentIndexChanged.connect(lambda _i: self.inspect())
        self.rs_label = QLabel()
        self.roi_table = make_table(["Name", "DICOM-Typ", "Kategorie", "Volumen [cm³]"])
        inputs = QGroupBox("RTSTRUCT")
        box = QVBoxLayout(inputs)
        box.addWidget(self.rs_combo)
        box.addWidget(self.rs_label)
        box.addWidget(self.roi_table, 1)

        self.gallery = Gallery()
        self.stats_view = report_view()
        split = QSplitter(Qt.Orientation.Vertical)
        split.addWidget(self.gallery)
        split.addWidget(self.stats_view)
        self.build([(inputs, 1)], SettingsForm(structures.Settings), [(split, 1)])

    def set_case(self, case) -> None:
        self.rs_combo.blockSignals(True)
        self.rs_combo.clear()
        for p in case.rs:
            self.rs_combo.addItem(p.name, str(p))
        self.rs_combo.blockSignals(False)
        super().set_case(case)

    def selection(self, case):
        self.roi_table.setRowCount(0)
        self.rs_label.clear()
        rs = self.rs_combo.currentData()
        if not rs:
            self.show_check("Kein RTSTRUCT im Datensatz-Ordner (RS*.dcm).", [], False)
            return None
        return selection.from_rtstruct(rs, case_id=case.folder.name)

    def on_inspected(self, info) -> None:
        rt = info.rtstruct
        if rt:
            self.rs_label.setText(f"{rt.get('structure_set_label') or '-'}   ·   {len(rt['rois'])} ROIs"
                                  f"   ·   {len(rt.get('markers', []))} Marker")
            fill_table(self.roi_table, [
                [r["name"], r.get("rt_type", ""), CATEGORY_DE.get(r["category"], r["category"]),
                 "" if r.get("volume_cm3") is None else f"{r['volume_cm3']:.2f}"]
                for r in rt["rois"]])
        self.form.set_placeholder("targets", "automatisch: " + names_text(info.auto_targets))
        self.form.set_placeholder("oars", "automatisch: " + names_text(info.auto_oars))

    def check(self, settings) -> tuple:
        pv = structures.preview(self.info, settings)
        text = (f"Zielvolumen: {names_text(pv['targets'])}\nRisikoorgane: {names_text(pv['oars'])}\n"
                f"Hilfsstrukturen: {len(pv['helpers'])}")
        return text, pv["issues"], pv["ok"]

    def clear_outputs(self) -> None:
        self.gallery.clear()
        self.stats_view.clear()

    def show_outputs(self, result: dict, out) -> None:
        outputs = result.get("outputs", {})
        self.gallery.set_images([out / rel for role, rel in outputs.items()
                                 if role.startswith("plot:") and str(rel).lower().endswith(".png")])
        stats = outputs.get("statistics")
        if stats and (out / stats).is_file():
            self.stats_view.setPlainText((out / stats).read_text(encoding="utf-8", errors="replace"))
