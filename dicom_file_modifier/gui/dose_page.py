"""
Seite Dosisindizes: RS + RD (+ RP, CT, eclipse_ref.json) aus dem Datensatz ->
Ziele, Verschreibung, Isodosen, Modus -> Pruefung (``dose.preview``: effektive
Werte, Feingitter, ROI-Namen) -> Start im Worker -> Kennzahlen je Ziel mit
Eclipse-Abgleich, Bericht, dose_overview.png und Validierungsansicht.
"""

from __future__ import annotations

import json
from pathlib import Path

from PySide6.QtGui import QColor
from PySide6.QtWidgets import QGroupBox, QLabel, QPushButton, QTabWidget, QVBoxLayout, QWidget

from ..api import dose, selection
from ..api.sysinfo import format_bytes
from .page import WorkflowPage, names_text
from .widgets import ImageViewer, SettingsForm, fill_table, make_table, open_path, report_view, std_icon

RX_SOURCE_DE = {"rtplan": "aus dem RTPLAN", "cli": "manuell", "pct_of_max": "in % von Dmax"}
ECLIPSE_SOURCE_DE = {"dvh": "DVH der RTDOSE", "json": "eclipse_ref.json", "cli": "Eclipse-Werte",
                     "derived": "abgeleitet"}
SUMMATION_DE = {"PLAN": "Plan-Summe", "MULTI_PLAN": "Summe mehrerer Pläne", "BEAM": "Einzelfeld",
                "FRACTION": "eine Fraktion"}
# Zeilen der Kennzahltabelle in der Reihenfolge des Eclipse-Abgleichs: Anzeige, Schluessel, Stellen
METRICS = [("TV [cm³]", "tv_cm3", 3), ("PIV [cm³]", "piv_cm3", 3), ("TV∩PIV [cm³]", "tv_piv_cm3", 3),
           ("PIV50 [cm³]", "piv50_cm3", 3), ("CI Paddick", "ci_paddick", 3), ("GI", "gi", 2),
           ("HI ICRU 83", "hi_icru83", 3), ("D98 [Gy]", "d98_gy", 2), ("D50 [Gy]", "d50_gy", 2),
           ("D2 [Gy]", "d2_gy", 2), ("Dmean [Gy]", "dmean_gy", 2), ("Dmin [Gy]", "dmin_gy", 2),
           ("Dmax [Gy]", "dmax_gy", 2)]
FLAG_COLOR = QColor("#c62828")


def _short(name: str, n: int = 44) -> str:
    return name if len(name) <= n else name[:n - 1] + "…"


def _num(v, digits: int):
    return ("–", None) if v is None else (f"{v:.{digits}f}", v)


def _diff(pct):
    return ("–", None) if pct is None else (f"{pct:+.1f}".replace("+0.0", "0.0").replace("-0.0", "0.0"), pct)


def _global_label(label: str) -> str:
    """"PIV [cm³]" -> "PIV global [cm³]"."""
    return label.replace(" [", " global [") if " [" in label else label + " global"


def _voxels(n: int) -> str:
    return f"{n / 1e6:.1f} M Voxel" if n >= 1e6 else f"{n / 1e3:.0f} k Voxel" if n >= 1e3 else f"{n} Voxel"


def _target_blocks(result: dict, out: Path) -> dict:
    """Je Ziel ``(Werte, Eclipse-Block)`` aus der Ergebnis-JSON; ohne JSON nur die Kurzfassung."""
    rel = (result.get("outputs") or {}).get("json")
    if rel:
        try:
            report = json.loads((out / rel).read_text(encoding="utf-8"))
            return {name: ({**t["components"], **t["indices"], **t["dvh_stats"]}, t.get("eclipse") or {})
                    for name, t in report["targets"].items()}
        except (OSError, ValueError, KeyError, TypeError):
            pass
    return {name: (t, {}) for name, t in ((result.get("summary") or {}).get("targets") or {}).items()}


class DosePage(WorkflowPage):
    title = "Dosisindizes"
    start_text = "Berechnung starten"
    api = dose

    def __init__(self, main):
        super().__init__(main)
        self.files_label = QLabel()
        self.files_label.setWordWrap(True)
        self.dose_label = QLabel()
        self.dose_label.setWordWrap(True)
        inputs = QGroupBox("Eingaben")
        box = QVBoxLayout(inputs)
        box.addWidget(self.files_label)
        box.addWidget(self.dose_label)

        metas = dose.Settings.field_meta()
        form = SettingsForm(dose.Settings, fields=[n for n, m in metas.items() if m.level == "basic"],
                            advanced=[n for n, m in metas.items() if m.level != "basic"])

        self.metrics = make_table(["Kennzahl"])
        self.metrics.horizontalHeader().setStretchLastSection(False)
        self.metrics_note = QLabel()
        self.metrics_note.setWordWrap(True)
        self.metrics_note.hide()
        metrics_tab = QWidget()
        lay = QVBoxLayout(metrics_tab)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self.metrics, 1)
        lay.addWidget(self.metrics_note)
        self.viz_button = QPushButton("Validierungsansicht")
        self.viz_button.setToolTip("validation.html im Browser öffnen")
        self.viz_button.setEnabled(False)
        self.viz_button.clicked.connect(lambda: self._viz and open_path(self._viz))
        self._viz = None
        self.gallery = ImageViewer()
        self.report = report_view()
        self.tabs = QTabWidget()
        self.tabs.addTab(metrics_tab, "Kennzahlen")
        self.tabs.addTab(self.gallery, "Dosisübersicht")
        self.tabs.addTab(self.report, "Bericht")
        self.build([(inputs, 0)], form, [(self.tabs, 1)], buttons=[self.viz_button])

    def selection(self, case):
        self.files_label.clear()
        self.dose_label.clear()
        return selection.for_dose(str(case.folder))

    def on_inspected(self, info) -> None:
        sel, d, labels = info.selection, info.dose, info.labels

        def name(p) -> str:
            return _short(Path(p).name) if p else "–"

        rs = labels.get("rs") or name(sel.rs)
        if info.rois:
            rs += f"  ·  {len(info.rois)} ROIs"
        rd = name(sel.rd)
        if d:
            spacing = " × ".join(f"{v:g}" for v in d["spacing_mm"])
            rd = (f"{SUMMATION_DE.get(d.get('summation_type'), d.get('summation_type') or 'Dosis')}  ·  "
                  f"Dmax {d['dmax_gy']:.2f} Gy  ·  Raster {spacing} mm")
        rp = (labels.get("rp") or name(sel.rp)) if sel.rp else "–"
        ct = f"{info.ct['n_slices']} Schichten" if info.ct else (info.rs_export.get("reason") or "keins")
        lines = [f"RS:  {rs}", f"RD:  {rd}", f"RP:  {rp}", f"CT:  {ct}"]
        if sel.eclipse_ref:
            lines.append(f"Eclipse-Referenz:  {name(sel.eclipse_ref)}")
        self.files_label.setText("\n".join(lines))
        self.files_label.setToolTip("\n".join(str(p) for p in (sel.rs, sel.rd, sel.rp, sel.eclipse_ref) if p))
        if not d:
            return
        dvh = info.dvh
        if not dvh.get("available"):
            dvh_text = "keine in der RTDOSE"
        elif not dvh.get("trusted"):
            dvh_text = "gehören zu einem anderen RTSTRUCT (nicht verwendet)"
        else:
            dvh_text = f"{len(dvh['rois'])} ROIs" + ("" if dvh.get("body") else ", ohne Körper-DVH (kein Eclipse-PIV)")
        self.dose_label.setText(f"DVH: {dvh_text}")
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
            lines.append(f"Feingitter: {e['grid_mm']:g} mm{aligned}  ·  {_voxels(g['n_voxels'])}"
                         f"  ·  ca. {format_bytes(g['memory_bytes'])} Speicher")
        if settings.eclipse_compat and pv.grid:
            lines.append(f"Eclipse-kompatibel ({settings.eclipse_compat}): Volumenmodell {e['volume_model']}, "
                         f"PIV {e['piv_scope']}, Interpolation {e['dose_interp']}, Konturen {e['iso_contours']}")
        if pv.grid:
            lines.append(f"Isodosen-RTSTRUCT: {len(pv.roi_names)} ROIs" if e.get("write_rs")
                         else "Isodosen-RTSTRUCT: wird nicht geschrieben")
        sources = pv.eclipse.get("sources") or []
        if sources:
            named = [ECLIPSE_SOURCE_DE.get(s, s) for s in sources if s != "derived"]
            text = ", ".join(named) or ECLIPSE_SOURCE_DE["derived"]
            if named and "derived" in sources:
                text += " (übrige Werte abgeleitet)"
            lines.append("Eclipse-Referenz: " + text)
        return "\n".join(lines) or "Prüfung nicht möglich.", pv.issues, pv.ok

    def clear_outputs(self) -> None:
        self.metrics.setRowCount(0)
        self.metrics_note.hide()
        self.gallery.clear()
        self.report.clear()
        self._viz = None
        self.viz_button.setEnabled(False)

    def show_outputs(self, result: dict, out: Path) -> None:
        self._fill_metrics(_target_blocks(result, out))
        outputs = result.get("outputs", {})
        if outputs.get("viz_html"):
            self._viz = out / outputs["viz_html"]
            self.viz_button.setEnabled(True)
        self.gallery.set_images([out / outputs["viz_png"]] if outputs.get("viz_png") else [])
        txt = outputs.get("txt")
        if txt and (out / txt).is_file():
            self.report.setPlainText((out / txt).read_text(encoding="utf-8", errors="replace"))
        self.tabs.setCurrentIndex(0)

    def _fill_metrics(self, blocks: dict) -> None:
        """
        Zeilen = Kennzahlen, Spalten je Ziel: Wert, bei Eclipse-Referenz auch Eclipse und
        Abweichung.  Hat der Abgleich einen anderen Wert verglichen (ganze Isodose bei
        PIV-Bereich component), steht er in einer eigenen Zeile "... global".
        """
        headers, cols = ["Kennzahl"], []
        for name, (values, ec) in blocks.items():
            rows = {r["key"]: r for r in ec.get("rows", [])} if ec.get("n_compared") else None
            headers += [_short(name, 24)] + (["Eclipse", "Abw. %"] if rows is not None else [])
            cols.append((values, rows))
        self.metrics.setColumnCount(len(headers))
        self.metrics.setHorizontalHeaderLabels(headers)
        empty = ("", None)
        data, marks, split_any = [], [], False
        for label, key, digits in METRICS:
            main, other, pending = [label], [_global_label(label)], []
            for values, rows in cols:
                shown = _num(values.get(key), digits)
                main.append(shown)
                other.append(empty)
                if rows is None:
                    continue
                r = rows.get(key) or {}
                cells = [_num(r.get("eclipse"), digits), _diff(r.get("diff_pct"))]
                compared = _num(r.get("tool"), digits)
                col = len(main) - 1
                if r.get("tool") is not None and compared[0] != shown[0]:
                    other[col] = compared
                    main += [empty, empty]
                    other += cells
                    pending.append((1, col, r, digits))
                else:
                    main += cells
                    other += [empty, empty]
                    pending.append((0, col, r, digits))
            split = any(which for which, *_ in pending)
            split_any |= split
            marks += [(len(data) + which, col, r, d) for which, col, r, d in pending
                      if r.get("diff_abs") is not None]
            data.append(main)
            if split:
                data.append(other)
        fill_table(self.metrics, data)
        for row, col, r, digits in marks:
            self._mark(row, col, r, digits)
        self.metrics.resizeColumnsToContents()                  # Platz fuer die Warnsymbole

        notes = []
        for name, (_values, ec) in blocks.items():
            if ec.get("n_compared"):
                n, flagged, tol = ec["n_compared"], ec.get("n_flagged") or 0, ec.get("tol_pct", 5.0)
                notes.append(f"{name}: {n} Werte mit Eclipse verglichen, "
                             + (f"{flagged} außerhalb der Toleranz ({tol:g} %)" if flagged
                                else f"alle innerhalb der Toleranz ({tol:g} %)"))
        if split_any:
            notes.append(f"Zeilen „global“: ganze Isodose, so vergleicht Eclipse; die übrigen Werte des Tools gelten "
                         f"für „{self.form.label_of('piv_scope')}“ = component (nur Isodosen-Anteile am Ziel).")
        self.metrics_note.setText("\n".join(notes))
        self.metrics_note.setVisible(bool(notes))

    def _mark(self, row: int, col: int, r: dict, digits: int) -> None:
        """Vergleichszeile: Tooltip mit den verglichenen Werten, ausserhalb der Toleranz rot mit Symbol."""
        def fmt(v):
            return "–" if v is None else f"{v:.{digits}f}"

        src = ECLIPSE_SOURCE_DE.get(r.get("source"), r.get("source") or "–")
        tip = f"Tool {fmt(r.get('tool'))} gegen Eclipse {fmt(r.get('eclipse'))}  ·  Quelle: {src}"
        if r.get("note"):
            tip += f"  ·  {r['note']}"
        for c in (col + 1, col + 2):
            self.metrics.item(row, c).setToolTip(tip)
        if r.get("within_tol") is False:
            self.metrics.item(row, col + 2).setIcon(std_icon(self, "SP_MessageBoxWarning"))
            for c in (col, col + 1, col + 2):
                self.metrics.item(row, c).setForeground(FLAG_COLOR)

    def summary_text(self, result: dict) -> str:
        s = result.get("summary") or {}
        if not s:
            return ""
        n = len(s.get("targets", {}))
        parts = [f"{n} Zielvolumen"]
        if s.get("rx_gy") is not None:
            parts.insert(0, f"Rx {s['rx_gy']:.2f} Gy")
        if s.get("n_voxels"):
            parts.append(f"Feingitter {_voxels(s['n_voxels'])}")
        return "  ·  ".join(parts)
