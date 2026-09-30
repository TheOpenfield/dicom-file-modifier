"""
Seite Transformation: CT und RTSTRUCT gemeinsam (``case-transform``) oder nur das
CT (``ct-transform``) -> Verschiebung und Rotation an anatomisch beschrifteten
Achsen, Drehpunkt (Volumenmitte, Marker, Koordinate) -> Pruefung
(``transform.preview``: Bewegung im Klartext, Drehpunkt, Methode und FoR,
Clipping, Speicher) -> Start im Worker -> Vorher/Nachher, Verschiebung je ROI,
Bericht (Matrix, Clipping, Schwerpunkt-Pruefung), 3D-Ansicht.
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QCheckBox, QComboBox, QGridLayout, QGroupBox, QHBoxLayout, QLabel,
                               QPushButton, QTabWidget, QVBoxLayout, QWidget)

from ..api import selection, transform
from ..api.sysinfo import format_bytes
from .page import WorkflowPage
from .widgets import ImageViewer, SettingsForm, open_path, report_view

SHIFTS, ROTATIONS = ("tx", "ty", "tz"), ("rx", "ry", "rz")
AXES = ("X  (+ links)", "Y  (+ posterior)", "Z  (+ superior)")
MAX_SHIFT_MM = 1000.0
METHOD_DE = {"resample": "Neuabtastung", "metadata": "nur Lage-Tags (Pixel unverändert)"}
ORDER_DE = {0: "nächster Nachbar", 1: "linear", 3: "kubisch"}
FOR_DE = {"keep": "FoR beibehalten", "new": "neue FoR"}
FOR_SHORT = {"keep": "beibehalten", "new": "neu"}
FOR_KEPT_DE = "beibehalten (das TPS legt vorhandene Pläne und Dosen darüber)"


def _mm(v) -> str:
    return ", ".join(f"{x:.1f}" for x in v)


def _center_text(spec: str) -> str:
    """Drehpunkt-Angabe fuer die Anzeige: Volumenmitte, Marker NAME oder Koordinate."""
    s = (spec or "").strip()
    if s.lower() == "volume":
        return "Volumenmitte"
    if s.lower().startswith("marker:"):
        return f"Marker {s.split(':', 1)[1].strip()}"
    return "Koordinate"


class TransformPage(WorkflowPage):
    title = "Transformation"
    start_text = "Transformation starten"
    api = transform
    summary_fields = SHIFTS + ROTATIONS + ("center",)    # stehen schon in der Kurzfassung

    def __init__(self, main):
        super().__init__(main)
        self._center_before = ""
        self.rs_combo = QComboBox()
        self.rs_combo.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.rs_combo.setMinimumContentsLength(20)
        self.rs_combo.currentIndexChanged.connect(lambda _i: self.inspect())
        self.ct_only = QCheckBox("Nur das CT transformieren (ohne RTSTRUCT)")
        self.ct_only.toggled.connect(self._toggle_ct_only)
        self.case_label = QLabel()
        self.case_label.setWordWrap(True)
        inputs = QGroupBox("Eingaben")
        box = QVBoxLayout(inputs)
        for w in (self.rs_combo, self.ct_only, self.case_label):
            box.addWidget(w)

        metas = transform.Settings.field_meta()
        own = SHIFTS + ROTATIONS + ("center",)
        form = SettingsForm(transform.Settings,
                            fields=[n for n, m in metas.items() if m.level == "basic" and n not in own],
                            advanced=[n for n, m in metas.items() if m.level != "basic"], external=own)
        form.add_top(self._motion_box(form))

        self.overview = ImageViewer()
        self.displacement = ImageViewer()
        self.report = report_view()
        self.tabs = QTabWidget()
        self.tabs.addTab(self.overview, "Vorher/Nachher")
        self.tabs.addTab(self.displacement, "Verschiebung je ROI")
        self.tabs.addTab(self.report, "Bericht")
        self.view3d_button = QPushButton("3D-Ansicht")
        self.view3d_button.setToolTip("Interaktive Vorher/Nachher-Ansicht im Browser öffnen")
        self.view3d_button.setEnabled(False)
        self.view3d_button.clicked.connect(lambda: self._html and open_path(self._html))
        self._html = None
        self.build([(inputs, 0)], form, [(self.tabs, 1)], buttons=[self.view3d_button])

    def _motion_box(self, form: SettingsForm) -> QWidget:
        """Verschiebung und Rotation als Raster mit Achsen im Patientensystem, darunter der Drehpunkt."""
        box = QWidget()
        grid = QGridLayout(box)
        grid.setContentsMargins(0, 0, 0, 6)
        for j, axis in enumerate(AXES):
            head = QLabel(axis)
            head.setAlignment(Qt.AlignmentFlag.AlignCenter)
            grid.addWidget(head, 0, j + 1)
        for i, (names, text, suffix) in enumerate(((SHIFTS, "Verschiebung", " mm"), (ROTATIONS, "Rotation", " °"))):
            grid.addWidget(QLabel(text), i + 1, 0)
            for j, name in enumerate(names):
                w = form.widget(name)
                w.setSuffix(suffix)
                w.setSingleStep(0.5)
                if suffix == " mm":
                    w.setRange(-MAX_SHIFT_MM, MAX_SHIFT_MM)
                grid.addWidget(w, i + 1, j + 1)
        self.center_combo = QComboBox()
        self.center_combo.currentIndexChanged.connect(self._center_mode)
        self.center_edit = form.widget("center")
        self.center_edit.setPlaceholderText("x,y,z in mm (LPS)")
        self.center_edit.textChanged.connect(self._sync_center_combo)
        row = QHBoxLayout()
        row.addWidget(self.center_combo, 1)
        row.addWidget(self.center_edit, 1)
        grid.addWidget(QLabel("Drehpunkt"), 3, 0)
        grid.addLayout(row, 3, 1, 1, 3)
        self._set_centers([])
        return box

    # -- Drehpunkt: Auswahl schreibt die Angabe ("volume", "marker:NAME", "x,y,z") ins Feld --
    def _find_center(self, spec: str) -> int:
        """Eintrag der Angabe (Gross-/Kleinschreibung egal wie im Kern), sonst -1."""
        s = spec.strip().lower()
        return next((k for k in range(self.center_combo.count())
                     if (self.center_combo.itemData(k) or "").lower() == s), -1)

    def _center_display(self, spec: str) -> str:
        """Anzeige mit dem Markernamen, wie er im RTSTRUCT steht (auch bei getipptem "marker:hs1")."""
        i = self._find_center(spec)
        return _center_text((self.center_combo.itemData(i) if i >= 0 else None) or spec)

    def _set_centers(self, markers: list) -> None:
        self.center_combo.blockSignals(True)
        self.center_combo.clear()
        self.center_combo.addItem("Volumenmitte", "volume")
        for m in markers:
            self.center_combo.addItem(f"Marker {m['name']}  ({_mm(m['position_mm'])} mm)", f"marker:{m['name']}")
        self.center_combo.addItem("Koordinate …", None)
        self.center_combo.blockSignals(False)
        spec = self.center_edit.text().strip()
        if spec.lower().startswith("marker:") and self._find_center(spec) < 0:
            self.center_edit.setText("volume")                  # Marker gibt es in diesem Datensatz nicht
        self._sync_center_combo(self.center_edit.text())

    def _sync_center_combo(self, text: str) -> None:
        i = self._find_center(text)
        if i < 0:
            i = self.center_combo.count() - 1                   # "Koordinate …"
        self.center_combo.blockSignals(True)
        self.center_combo.setCurrentIndex(i)
        self.center_combo.blockSignals(False)
        self.center_edit.setVisible(self.center_combo.currentData() is None)

    def _center_mode(self, i: int) -> None:
        spec = self.center_combo.itemData(i)
        if spec is not None:
            self.center_edit.setText(spec)
        else:
            text = self.center_edit.text().strip().lower()
            if (text == "volume" or text.startswith("marker:")) and self.info is not None and self.info.ct:
                self.center_edit.setText(",".join(f"{v:.1f}" for v in self.info.ct["volume_center_mm"]))
            self.center_edit.show()
            self.center_edit.setFocus()

    def _toggle_ct_only(self, on: bool) -> None:
        if on:                                                  # ohne RTSTRUCT immer die Volumenmitte
            self._center_before = self.center_edit.text()
            self.center_edit.setText("volume")
        elif self._center_before:
            self.center_edit.setText(self._center_before)       # Wahl von vorher zurueck
        self.inspect()

    # -- Seite ---------------------------------------------------------------------------
    def set_case(self, case) -> None:
        self.rs_combo.blockSignals(True)
        self.rs_combo.clear()
        for p in case.rs:
            self.rs_combo.addItem(p.name, str(p))
        self.rs_combo.blockSignals(False)
        self.ct_only.blockSignals(True)
        self.ct_only.setChecked(not case.rs)
        self.ct_only.blockSignals(False)
        super().set_case(case)

    def selection(self, case):
        self.case_label.clear()
        only = self.ct_only.isChecked()
        self.rs_combo.setEnabled(not only)
        if only:
            return selection.from_ct_dir(str(case.folder / "CT"), case_id=case.folder.name)
        rs = self.rs_combo.currentData()
        if not rs:
            self.show_check("Kein RTSTRUCT im Datensatz-Ordner (RS*.dcm): „Nur das CT transformieren“ wählen.",
                            [], False)
            return None
        return selection.for_transform(str(case.folder), rs=rs)

    def on_inspected(self, info) -> None:
        lines = []
        ct = info.ct
        if ct:
            nz, ny, nx = ct["shape"]
            sp = ct["spacing_mm"]
            lines.append(f"CT:  {ct['n_slices']} Schichten  ·  {nx} × {ny} Pixel  ·  "
                         f"{sp['dc']:g} × {sp['dr']:g} × {sp['dz']:g} mm")
            lines.append(f"Volumenmitte:  {_mm(ct['volume_center_mm'])} mm")
        rt = info.rtstruct
        if info.has_rs and rt:
            lines.append(f"RS:  {rt.get('label') or rt['file']}  ·  {rt['n_rois']} ROIs  ·  {len(info.markers)} Marker")
        self.case_label.setText("\n".join(lines))
        self._set_centers(info.markers)

    def check(self, settings) -> tuple:
        pv = transform.preview(self.info, settings)
        self._sync_center_combo(self.center_edit.text())       # auch nach set_settings
        self.center_combo.setEnabled(self.center_edit.isEnabled())
        s, with_rs = settings, self.info.has_rs
        order = " (Drehreihenfolge X, Y, Z, intrinsisch)" if sum(1 for v in (s.rx, s.ry, s.rz) if v) > 1 else ""
        lines = [f"Bewegung: {pv.description}{order}"]
        if pv.center_mm:
            lines.append(f"Drehpunkt: {self._center_display(s.center) if with_rs else 'Volumenmitte'} "
                         f"({_mm(pv.center_mm)} mm)")
        if pv.drehpunkt_mm:
            lines.append(f"Neuer POINT-Marker „Drehpunkt“ im Ergebnis bei ({_mm(pv.drehpunkt_mm)} mm)")
        method = METHOD_DE.get(s.method, s.method) + (f", {ORDER_DE.get(s.order, s.order)}" if s.method == "resample"
                                                      else "")
        new_for = with_rs and s.new_frame_of_reference
        lines.append(f"Methode: {method}  ·  FrameOfReference: {'neu' if new_for else FOR_KEPT_DE}")
        if pv.clipping:
            top = sorted(pv.clipping, key=lambda c: -c["fraction"])
            lines.append("Clipping: " + ", ".join(f"{c['roi']} {100 * c['fraction']:.1f} %" for c in top[:3])
                         + (f" und {len(top) - 3} weitere" if len(top) > 3 else "")
                         + " der Konturpunkte außerhalb des CT")
        if pv.planned:
            what = "CT und RTSTRUCT gemeinsam" if pv.tool == transform.CASE_TOOL else "nur das CT"
            lines.append(f"Schreibt: {what} ({', '.join(pv.planned.values())})")
        if pv.memory_bytes:
            lines.append(f"Speicherbedarf: ca. {format_bytes(pv.memory_bytes)}")
        return "\n".join(lines), pv.issues, pv.ok

    def clear_outputs(self) -> None:
        self.overview.clear()
        self.displacement.clear()
        self.report.clear()
        self._html = None
        self.view3d_button.setEnabled(False)

    def show_outputs(self, result: dict, out: Path) -> None:
        outputs = result.get("outputs", {})
        views = [out / rel for role, rel in outputs.items() if (role.startswith("viz:") or role == "viz_html") and rel]
        pngs = [p for p in views if p.suffix.lower() == ".png"]
        with_rs = (result.get("summary") or {}).get("tool") == transform.CASE_TOOL
        for viewer, stem in ((self.overview, "overview"), (self.displacement, "displacement")):
            if with_rs:
                viewer.set_images([p for p in pngs if stem in p.stem])
            else:
                viewer.clear("Nur mit RTSTRUCT (Konturen vorher/nachher); siehe 3D-Ansicht und Bericht")
        self.report.setPlainText(self._report_text(result))
        self._html = next((p for p in views if p.suffix.lower() == ".html" and p.is_file()), None)
        self.view3d_button.setEnabled(self._html is not None)
        self.tabs.setCurrentWidget(self.overview if self.overview.count() else self.report)

    def _center_of(self, result: dict) -> str:
        s, js = result.get("summary") or {}, self._job_settings
        if s.get("tool") == transform.CASE_TOOL and js is not None:
            return self._center_display(js.center)
        return s.get("rotation_center_label") or "Volumenmitte"

    def _report_text(self, result: dict) -> str:
        """Bericht des Laufs aus der Kurzfassung und den Einstellungen beim Start."""
        s, js = result.get("summary") or {}, self._job_settings
        if not s:
            return ""
        motion = (s.get("description") or "").split(" · ")      # Bericht ohne Zeilenumbruch: ein Teil je Zeile
        L = [f"Bewegung      {motion[0]}"] + [f"              {m}" for m in motion[1:]]
        if s.get("rotation_center_mm"):
            L.append(f"Drehpunkt     {self._center_of(result)} ({_mm(s['rotation_center_mm'])} mm)")
        if s.get("drehpunkt_mm"):
            L.append(f"POINT-Marker  „Drehpunkt“ im Ergebnis bei ({_mm(s['drehpunkt_mm'])} mm)")
        if js is not None:
            L.append(f"Methode       {METHOD_DE.get(js.method, js.method)}"
                     + (f", Interpolation {ORDER_DE.get(js.order, js.order)}" if js.method == "resample" else ""))
        L.append(f"CT-Schichten  {s.get('n_slices', 0)}")
        L.append(f"FoR           {FOR_SHORT.get(s.get('for_strategy'), 'beibehalten')}")
        if s.get("T"):
            L += ["", "Matrix T (Patient -> Patient, LPS, mm):"]
            L += ["  " + "".join(f"{v:12.6f}" for v in row) for row in s["T"]]
        clip = s.get("clipping") or []
        if clip:
            L += ["", "Konturpunkte außerhalb des CT nach der Bewegung:", f"  {'ROI':<28}{'Punkte':>16}{'Anteil':>10}"]
            L += [f"  {c['roi']:<28}{c['n_outside']:>8} / {c['n_total']:<7}{100 * c['fraction']:>8.1f} %"
                  for c in sorted(clip, key=lambda c: -c["fraction"])]
        v = s.get("verify")
        if v:
            L += ["", f"Schwerpunkt-Prüfung: {v['checked']} ROIs, max. Abweichung {v['max_err_mm']:.2e} mm "
                      f"({v.get('worst_roi') or '-'})",
                  f"  Schwelle {v['threshold_mm']:g} mm -> {'bestanden' if v['passed'] else 'NICHT bestanden'}"]
        return "\n".join(L)

    def summary_text(self, result: dict) -> str:
        s = result.get("summary") or {}
        if not s:
            return ""
        parts = [s.get("description", ""), f"Drehpunkt {self._center_of(result)}", f"{s.get('n_slices', 0)} CT-Schichten"]
        if s.get("for_strategy"):
            parts.append(FOR_DE.get(s["for_strategy"], s["for_strategy"]))
        v = s.get("verify")
        if v:
            parts.append(f"Schwerpunkte {'ok' if v['passed'] else 'ABWEICHUNG'} (max {v['max_err_mm']:.1e} mm)")
        if s.get("n_clipped_rois"):
            parts.append(f"{s['n_clipped_rois']} ROIs ragen aus dem CT")
        return "  ·  ".join(p for p in parts if p)
