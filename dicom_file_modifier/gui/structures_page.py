"""
Seite Strukturanalyse: RTSTRUCT waehlen -> ROI-Tabelle (mit der Rolle im
geplanten Lauf) -> Einstellungen -> Pruefung (``structures.preview``) -> Start
im Worker -> Plots und statistics.txt.
"""

from __future__ import annotations

from PySide6.QtGui import QColor, QIcon, QPixmap
from PySide6.QtWidgets import QComboBox, QGroupBox, QLabel, QTabWidget, QVBoxLayout

from ..api import selection, structures
from .page import WorkflowPage, names_text
from .widgets import ImageViewer, SettingsForm, fill_table, make_table, report_view

CATEGORY_DE = {"TARGET": "Zielvolumen", "OAR_SERIAL": "Risikoorgan (seriell)",
               "OAR_PARALLEL": "Risikoorgan (parallel)", "HELPER": "Hilfsstruktur",
               "EXTERNAL": "Außenkontur", "MARKER": "Marker"}
ROLE_COLUMN = 4


def _swatch(color) -> QIcon:
    """Farbfeld der ROI (``ROIDisplayColor``)."""
    pix = QPixmap(12, 12)
    try:
        pix.fill(QColor(*[int(c) for c in color][:3]))
    except (TypeError, ValueError):
        pix.fill(QColor(0, 0, 0, 0))
    return QIcon(pix)


class StructuresPage(WorkflowPage):
    title = "Strukturanalyse"
    start_text = "Analyse starten"
    api = structures

    def __init__(self, main):
        super().__init__(main)
        self.rs_combo = QComboBox()
        self.rs_combo.currentIndexChanged.connect(lambda _i: self.inspect())
        self.rs_label = QLabel()
        self.roi_table = make_table(["Name", "DICOM-Typ", "Kategorie", "Volumen [cm³]", "Rolle im Lauf"],
                                    sortable=True)
        inputs = QGroupBox("RTSTRUCT")
        box = QVBoxLayout(inputs)
        box.addWidget(self.rs_combo)
        box.addWidget(self.rs_label)
        box.addWidget(self.roi_table, 1)

        self.gallery = ImageViewer()
        self.stats_view = report_view()
        self.tabs = QTabWidget()
        self.tabs.addTab(self.gallery, "Plots")
        self.tabs.addTab(self.stats_view, "Statistik")
        self.build([(inputs, 1)], SettingsForm(structures.Settings), [(self.tabs, 1)])

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
            rois = rt["rois"]
            fill_table(self.roi_table, [
                [r["name"], r.get("rt_type", ""), CATEGORY_DE.get(r["category"], r["category"]),
                 ("", None) if r.get("volume_cm3") is None else (f"{r['volume_cm3']:.2f}", r["volume_cm3"]),
                 ""]
                for r in rois])
            colors = {r["name"]: r.get("color") for r in rois}
            for i in range(self.roi_table.rowCount()):
                item = self.roi_table.item(i, 0)
                item.setIcon(_swatch(colors.get(item.text())))
            self.roi_table.resizeColumnToContents(0)
        self.form.set_placeholder("targets", "automatisch: " + names_text(info.auto_targets))
        self.form.set_placeholder("oars", "automatisch: " + names_text(info.auto_oars))

    def _show_roles(self, pv: dict) -> None:
        """Spalte "Rolle im Lauf" aus der Pruefung."""
        role = {n: "Zielvolumen" for n in pv["targets"]}
        role.update({n: "Risikoorgan" for n in pv["oars"]})
        role.update({n: "Hilfsstruktur" for n in pv["helpers"]})
        sorting = self.roi_table.isSortingEnabled()
        self.roi_table.setSortingEnabled(False)
        for i in range(self.roi_table.rowCount()):
            name = self.roi_table.item(i, 0).text()
            self.roi_table.item(i, ROLE_COLUMN).setText(role.get(name, "nicht ausgewertet"))
        self.roi_table.setSortingEnabled(sorting)

    def check(self, settings) -> tuple:
        pv = structures.preview(self.info, settings)
        self._show_roles(pv)
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
        self.tabs.setCurrentIndex(0)

    def summary_text(self, result: dict) -> str:
        s = result.get("summary") or {}
        if not s:
            return ""
        text = (f"{s.get('n_targets', 0)} Zielvolumen  ·  {s.get('n_oars', 0)} Risikoorgane  ·  "
                f"{s.get('n_helpers', 0)} Hilfsstrukturen")
        # kleinster Abstand Zielvolumen - Risikoorgan (nicht Ziel - Ziel wie GTV im PTV)
        cat = {r["name"]: r["category"] for r in (self.info.rtstruct.get("rois", []) if self.info else [])}
        for c in s.get("closest") or []:
            kinds = {cat.get(c["a"], ""), cat.get(c["b"], "")}
            if "TARGET" in kinds and kinds & {"OAR_SERIAL", "OAR_PARALLEL"} and c.get("min_distance_mm") is not None:
                text += f"  ·  kleinster Abstand Ziel – Risikoorgan: {c['a']} – {c['b']} {c['min_distance_mm']:.1f} mm"
                break
        return text
