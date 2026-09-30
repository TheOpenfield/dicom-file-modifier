"""
Seite Dosisindizes: RS + RD (+ RP, CT, eclipse_ref.json) aus dem Datensatz ->
Ziele, Verschreibung, Isodosen, Modus -> Pruefung (``dose.preview``: effektive
Werte, Feingitter, ROI-Namen) -> Start im Worker -> Kennzahlen je Ziel,
Bericht, dose_overview.png und Validierungsansicht.
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import QGroupBox, QLabel, QPushButton, QTabWidget, QVBoxLayout

from ..api import dose, selection
from ..api.sysinfo import format_bytes
from .page import WorkflowPage, names_text
from .widgets import ImageViewer, SettingsForm, fill_table, make_table, open_path, report_view

RX_SOURCE_DE = {"rtplan": "aus dem RTPLAN", "cli": "manuell", "pct_of_max": "in % von Dmax"}
ECLIPSE_SOURCE_DE = {"dvh": "DVH der RTDOSE", "json": "eclipse_ref.json", "cli": "Eclipse-Werte",
                     "derived": "abgeleitet"}
METRICS = [("Ziel", None, None), ("TV [cm³]", "tv_cm3", 3), ("PIV [cm³]", "piv_cm3", 3),
           ("TV∩PIV [cm³]", "tv_piv_cm3", 3), ("CI Paddick", "ci_paddick", 3), ("GI", "gi", 2),
           ("HI ICRU 83", "hi_icru83", 3), ("D98 [Gy]", "d98_gy", 2), ("D50 [Gy]", "d50_gy", 2),
           ("D2 [Gy]", "d2_gy", 2), ("Eclipse-Abgleich", None, None)]


def _short(name: str, n: int = 44) -> str:
    return name if len(name) <= n else name[:n - 1] + "…"


def _num(v, digits: int):
    return ("–", None) if v is None else (f"{v:.{digits}f}", v)


def _eclipse_cell(t: dict) -> str:
    n = t.get("eclipse_n_compared")
    if not n:
        return "–"
    flagged = t.get("eclipse_n_flagged") or 0
    return f"ok ({n} Werte)" if not flagged else f"{flagged} von {n} außerhalb"


class DosePage(WorkflowPage):
    title = "Dosisindizes"
    start_text = "Berechnung starten"
    api = dose

    def __init__(self, main):
        super().__init__(main)
        self.files_label = QLabel()
        self.dose_label = QLabel()
        self.dose_label.setWordWrap(True)
        inputs = QGroupBox("Eingaben")
        box = QVBoxLayout(inputs)
        box.addWidget(self.files_label)
        box.addWidget(self.dose_label)

        metas = dose.Settings.field_meta()
        form = SettingsForm(dose.Settings, fields=[n for n, m in metas.items() if m.level == "basic"],
                            advanced=[n for n, m in metas.items() if m.level != "basic"])

        self.metrics = make_table([m[0] for m in METRICS])
        fill_table(self.metrics, [], fit=True)
        self.viz_button = QPushButton("Validierungsansicht")
        self.viz_button.setToolTip("validation.html im Browser öffnen")
        self.viz_button.setEnabled(False)
        self.viz_button.clicked.connect(lambda: self._viz and open_path(self._viz))
        self._viz = None
        self.gallery = ImageViewer()
        self.report = report_view()
        self.tabs = QTabWidget()
        self.tabs.addTab(self.gallery, "Dosisübersicht")
        self.tabs.addTab(self.report, "Bericht")
        self.build([(inputs, 0)], form, [(self.metrics, 0), (self.tabs, 1)], buttons=[self.viz_button])

    def selection(self, case):
        self.files_label.clear()
        self.dose_label.clear()
        return selection.for_dose(str(case.folder))

    def on_inspected(self, info) -> None:
        sel = info.selection

        def name(p) -> str:
            return _short(Path(p).name) if p else "–"

        if info.ct:
            ct = f"{info.ct['n_slices']} Schichten"
        else:
            ct = info.rs_export.get("reason") or "keins"
        lines = [f"RS:  {name(sel.rs)}", f"RD:  {name(sel.rd)}", f"RP:  {name(sel.rp)}", f"CT:  {ct}"]
        if sel.eclipse_ref:
            lines.append(f"Eclipse-Referenz:  {name(sel.eclipse_ref)}")
        self.files_label.setText("\n".join(lines))
        self.files_label.setToolTip("\n".join(str(p) for p in (sel.rs, sel.rd, sel.rp, sel.eclipse_ref) if p))
        d, dvh = info.dose, info.dvh
        if not d:
            return
        if not dvh.get("available"):
            dvh_text = "keine DVH in der RTDOSE"
        elif not dvh.get("trusted"):
            dvh_text = "DVH gehören zu einem anderen RTSTRUCT (nicht verwendet)"
        else:
            dvh_text = f"DVH für {len(dvh['rois'])} ROIs"
        spacing = " × ".join(f"{v:g}" for v in d["spacing_mm"])
        self.dose_label.setText(f"Dosis: Dmax {d['dmax_gy']:.2f} Gy  ·  Raster {spacing} mm  ·  {dvh_text}")
        self.form.set_placeholder("target", "automatisch: " + names_text(info.targets.get("default", [])))
        self.form.set_placeholder("rx_pct_of_max", "nicht verwendet")
        rx = info.prescriptions.get("default")
        self.form.set_placeholder("rx", f"automatisch: {rx['target_prescription_dose_gy']:.2f} Gy aus dem RTPLAN"
                                  if rx else "nicht im RTPLAN gefunden: bitte angeben")

    def check(self, settings) -> tuple:
        pv = dose.preview(self.info, settings)
        lines = []
        if pv.targets:
            lines.append(f"Zielvolumen: {names_text(pv.targets)}")
        if pv.rx:
            lines.append(f"Verschreibung: {pv.rx['gy']:.2f} Gy ({RX_SOURCE_DE.get(pv.rx['source'], pv.rx['source'])})")
        if pv.levels:
            lines.append("Isodosen: " + ", ".join(f"{lv['label']} = {lv['gy']:.2f} Gy" for lv in pv.levels))
        e = pv.effective
        if pv.grid:
            g = pv.grid
            aligned = " am CT-Pixelraster" if e.get("align") else ""
            lines.append(f"Feingitter: {e['grid_mm']:g} mm{aligned}  ·  {g['n_voxels'] / 1e6:.1f} M Voxel"
                         f"  ·  ca. {format_bytes(g['memory_bytes'])} Speicher")
        if settings.eclipse_compat and pv.grid:
            lines.append(f"Eclipse-kompatibel ({settings.eclipse_compat}): Volumenmodell {e['volume_model']}, "
                         f"PIV {e['piv_scope']}, Interpolation {e['dose_interp']}, Konturen {e['iso_contours']}")
        if pv.grid:
            lines.append(f"Isodosen-RTSTRUCT: {len(pv.roi_names)} ROIs" if e.get("write_rs")
                         else "Isodosen-RTSTRUCT: wird nicht geschrieben")
        sources = pv.eclipse.get("sources") or []
        if sources:
            lines.append("Eclipse-Referenz: " + ", ".join(ECLIPSE_SOURCE_DE.get(s, s) for s in sources))
        return "\n".join(lines) or "Prüfung nicht möglich.", pv.issues, pv.ok

    def clear_outputs(self) -> None:
        self.metrics.setRowCount(0)
        self.gallery.clear()
        self.report.clear()
        self._viz = None
        self.viz_button.setEnabled(False)

    def show_outputs(self, result: dict, out: Path) -> None:
        targets = (result.get("summary") or {}).get("targets", {})
        fill_table(self.metrics, [
            [name] + [_num(t.get(key), digits) for _, key, digits in METRICS[1:-1]] + [_eclipse_cell(t)]
            for name, t in targets.items()], fit=True)
        outputs = result.get("outputs", {})
        if outputs.get("viz_html"):
            self._viz = out / outputs["viz_html"]
            self.viz_button.setEnabled(True)
        self.gallery.set_images([out / outputs["viz_png"]] if outputs.get("viz_png") else [])
        txt = outputs.get("txt")
        if txt and (out / txt).is_file():
            self.report.setPlainText((out / txt).read_text(encoding="utf-8", errors="replace"))

    def summary_text(self, result: dict) -> str:
        s = result.get("summary") or {}
        if not s:
            return ""
        n = len(s.get("targets", {}))
        parts = [f"{n} Zielvolumen"]
        if s.get("rx_gy") is not None:
            parts.insert(0, f"Rx {s['rx_gy']:.2f} Gy")
        if s.get("n_voxels"):
            parts.append(f"Feingitter {s['n_voxels'] / 1e6:.1f} M Voxel")
        return "  ·  ".join(parts)
