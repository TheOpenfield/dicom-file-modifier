"""
Grundgeruest der Workflow-Seiten: Eingaben -> Einstellungen -> Pruefung ->
Start -> Ergebnis.  Die Seite liefert ``api`` (Workflow-Modul), ``selection``,
``on_inspected``, ``check``, ``clear_outputs``/``show_outputs`` und optional
``summary_text``; Inspektion im Hintergrund, Pruefung nach jeder Eingabepause
(``CHECK_DELAY_MS``, bis dahin ist Start gesperrt), Lauf ueber das Hauptfenster
(Worker).
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QTimer
from PySide6.QtWidgets import (QApplication, QGroupBox, QHBoxLayout, QLabel, QMenu, QPushButton,
                               QScrollArea, QSplitter, QToolButton, QVBoxLayout, QWidget)

from ..api import jobs as api_jobs
from ..api.fields import command_string
from ..api.issues import Issue, has_errors, issue_from_exception
from ..api.outputs import OutputSpec
from .icons import page_icon
from .jobs import run_in_background
from .widgets import ImageViewer, STATUS_DE, IssueList, StateLine, breakable, open_path

# Laengster Dateiname unter dem Ergebnisordner (z.B. RS.<64-Zeichen-UID>_analysis.json)
# plus Trenner: laengere Ordnerpfade stossen an MAX_PATH (260) von Windows
LONGEST_FILE_NAME = 100
MAX_PATH = 259
RESULT_LEVEL = {"ok": "ok", "ok_warnings": "warning", "failed": "error", "cancelled": "info"}
CHECK_DELAY_MS = 200                    # Pruefung erst nach dieser Eingabepause: eine je Tipp- oder Klickfolge


def _value_text(v, unit: str = "") -> str:
    if isinstance(v, bool):
        return "ja" if v else "nein"
    if isinstance(v, (list, tuple)):
        return ", ".join(map(str, v))
    if v is None:
        return "aus"
    if isinstance(v, (int, float)) and unit:
        return f"{v:g} {unit}"
    return str(v)


def names_text(names: list, limit: int = 8) -> str:
    if not names:
        return "keine"
    more = f" … (+{len(names) - limit})" if len(names) > limit else ""
    return ", ".join(names[:limit]) + more


class WorkflowPage(QWidget):
    title = ""
    start_text = "Starten"
    api = None                          # Workflow-Modul: api.structures, api.dose, api.transform
    summary_fields = ()                 # Felder, die schon die Kurzfassung nennt (nicht in der Lauf-Zeile)

    def __init__(self, main):
        super().__init__()
        self.main = main                # MainWindow
        self.case = self.info = self.last_result = self.form = self._left = None
        self._token = 0
        self._running = self._ok = False
        self._job_settings = None
        self._job_disabled: set = set()
        self._check_timer = QTimer(self)
        self._check_timer.setSingleShot(True)
        self._check_timer.setInterval(CHECK_DELAY_MS)
        self._check_timer.timeout.connect(self._update_preview)
        self.state = StateLine()
        self.preview_label = QLabel("Kein Datensatz geöffnet.")
        self.preview_label.setWordWrap(True)
        self.out_label = QLabel()
        self.out_label.setWordWrap(True)
        self.issues = IssueList()
        self.issues.max_height = 120
        self.issues.hide()
        self.start_button = QPushButton(self.start_text)
        self.start_button.setEnabled(False)
        self.start_button.setMinimumHeight(34)
        self.start_button.setMinimumWidth(200)
        font = self.start_button.font()
        font.setBold(True)
        self.start_button.setFont(font)
        self.start_button.clicked.connect(self._start)

        self.result_state = StateLine()
        self.result_state.set_state(None, "Noch kein Ergebnis.")
        self.summary_label = QLabel()
        self.summary_label.setWordWrap(True)
        self.summary_label.hide()
        self.open_button = QPushButton("Ordner öffnen")
        self.open_button.clicked.connect(lambda: open_path(self.last_result["output_dir"]))
        menu = QMenu(self)                              # Nachvollziehbarkeit, bewusst unauffaellig
        menu.setToolTipsVisible(True)
        self.command_action = menu.addAction("dfm-Befehl kopieren")
        self.command_action.triggered.connect(self._copy_command)
        self.runjson_action = menu.addAction("run.json öffnen")
        self.runjson_action.setToolTip("Protokoll des Laufs: Versionen, Einstellungen, Eingaben mit UIDs")
        self.runjson_action.triggered.connect(
            lambda: open_path(Path(self.last_result["output_dir"]) / "run.json"))
        self.more_button = QToolButton()
        self.more_button.setText("…")
        self.more_button.setToolTip("Weitere Aktionen")
        self.more_button.setMenu(menu)
        self.more_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        self.more_button.setMinimumHeight(self.open_button.sizeHint().height())
        for b in (self.open_button, self.more_button):
            b.setEnabled(False)
        self.result_issues = IssueList()
        self.result_issues.hide()

    def build(self, inputs: list, form, results: list, buttons: list = ()) -> None:
        """Layout: links (scrollbar) Eingaben, Einstellungen, Pruefung, Start; rechts das Ergebnis."""
        self.form = form
        form.changed.connect(self._schedule_check)
        self.issues.text_map = self.result_issues.text_map = form.gui_text
        left = QWidget()
        col = QVBoxLayout(left)
        for widget, stretch in inputs:
            col.addWidget(widget, stretch)
        box = QGroupBox("Einstellungen")
        QVBoxLayout(box).addWidget(form)
        col.addWidget(box)
        check = QGroupBox("Prüfung vor dem Start")
        lay = QVBoxLayout(check)
        for w in (self.preview_label, self.issues, self.out_label):
            lay.addWidget(w)
        col.addWidget(check)
        if not any(stretch for _, stretch in inputs):
            col.addStretch(1)                           # Gruppen oben kompakt halten
        self._left = left
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(left)
        side = QWidget()                                # Zustand und Start immer sichtbar unter der Spalte
        lay = QVBoxLayout(side)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(scroll, 1)
        foot = QHBoxLayout()
        foot.addWidget(self.state, 1)
        foot.addWidget(self.start_button)
        lay.addLayout(foot)

        result = QGroupBox("Ergebnis")
        lay = QVBoxLayout(result)
        lay.addWidget(self.result_state)
        lay.addWidget(self.summary_label)
        row = QHBoxLayout()                             # eigene Zeile: der Text behaelt die volle Breite
        for b in list(buttons) + [self.open_button, self.more_button]:
            row.addWidget(b)
        row.addStretch(1)
        lay.addLayout(row)
        lay.addWidget(self.result_issues)
        for widget, stretch in results:
            lay.addWidget(widget, stretch)

        outer = QSplitter()
        outer.addWidget(side)
        outer.addWidget(result)
        outer.setStretchFactor(1, 1)
        outer.setSizes([580, 740])
        QVBoxLayout(self).addWidget(outer)
        for viewer in self.findChildren(ImageViewer):          # Leerzustand mit dem Symbol der Seite
            viewer.set_placeholder_icon(page_icon(self.api.WORKFLOW))

    # -- von der Seite -----------------------------------------------------------------
    def selection(self, case):
        """``CaseSelection`` fuer diese Seite (darf werfen; ``None``: nichts zu tun)."""
        raise NotImplementedError

    def on_inspected(self, info) -> None:
        """Eingabe-Widgets aus dem Inspektionsergebnis fuellen."""

    def check(self, settings) -> tuple:
        """``(text, issues, ok)`` der Pruefung (``api.preview``)."""
        raise NotImplementedError

    def clear_outputs(self) -> None:
        """Ergebnis-Widgets leeren."""

    def show_outputs(self, result: dict, out: Path) -> None:
        """Ergebnis-Widgets aus dem ``JobResult``-dict fuellen."""

    def summary_text(self, result: dict) -> str:
        """Kurzfassung des Ergebnisses (eine Zeile)."""
        return ""

    # -- Datensatz und Inspektion ------------------------------------------------------
    def set_case(self, case) -> None:
        self.case = case
        self._reset_result("Noch kein Ergebnis.")      # das Ergebnis gehoerte zum vorigen Datensatz
        self.inspect()

    def inspect(self) -> None:
        self._token += 1
        token, self.info = self._token, None
        try:
            sel = self.selection(self.case)
        except Exception as exc:  # noqa: BLE001 - Befund statt Absturz
            self.show_check("Datensatz für diese Seite nicht verwendbar.", [issue_from_exception(exc)], False)
            return
        if sel is None:
            return
        self.show_check("Datensatz wird gelesen …", [], False, pending=True)
        run_in_background(lambda: self.api.inspect(sel),
                          lambda info, err: self._on_inspected(token, info, err))

    def _on_inspected(self, token: int, info, err) -> None:
        if token != self._token:
            return                                      # veraltet: inzwischen andere Auswahl
        if err is not None:
            self.show_check("Datensatz nicht lesbar.", [issue_from_exception(err)], False)
            return
        self.info = info
        self.on_inspected(info)
        self._update_preview()

    # -- Pruefung ----------------------------------------------------------------------
    def output_spec(self) -> OutputSpec:
        return OutputSpec(str(self.main.config.results_root))

    def _schedule_check(self) -> None:
        """Nach einer Eingabe: sofort "wird geprueft" (Start gesperrt), die Pruefung folgt nach der Eingabepause."""
        if self.info is None:
            return                                      # die Inspektion prueft, wenn sie fertig ist
        self._ok = False
        self.state.set_state(None, "Wird geprüft …")
        self._update_start()
        self._check_timer.start()

    def flush_check(self) -> None:
        """Eine anstehende Pruefung sofort ausfuehren (vor dem Start, in Tests)."""
        if self._check_timer.isActive():
            self._update_preview()

    def _update_preview(self) -> None:
        self._check_timer.stop()
        if self.info is None:
            return
        try:
            s = self.form.settings()
        except ValueError as exc:
            self.show_check("Eingabe prüfen.", [Issue("error", "INPUT.INVALID", str(exc))], False)
            return
        self.form.set_disabled(s.disabled_fields(self.info))
        text, issues, ok = self.check(s)
        out = self.output_spec()
        issues = list(issues) + out.validate()
        folder = None
        if ok and not has_errors(issues):
            folder = out.target_dir(self.api.default_folder(s, self.info.selection))
            if len(str(folder)) + LONGEST_FILE_NAME > MAX_PATH:
                issues.append(Issue(
                    "warning", "OUTPUT.PATH_LONG",
                    f"Ergebnispfad sehr lang ({len(str(folder))} Zeichen): Windows erlaubt 260 Zeichen "
                    "je Datei, lange DICOM-Dateinamen können das Schreiben scheitern lassen.",
                    hint_de="Kürzeren Ergebnis-Stammordner wählen (Symbolleiste: Ergebnisordner …)."))
        ok = ok and not has_errors(issues)
        self.form.mark_issues(issues)
        self.show_check(text, issues, ok)
        if folder is not None and ok:
            self.out_label.setText(f"Ergebnisordner: {breakable(folder)}")

    def show_check(self, text: str, issues: list, ok: bool, pending: bool = False) -> None:
        if pending:
            self.state.set_state(None, "Wird geprüft …")
        elif ok:
            warned = any((i.level if hasattr(i, "level") else i["level"]) == "warning" for i in issues)
            self.state.set_state("warning" if warned else "ok",
                                 "Bereit, mit Hinweisen" if warned else "Bereit zum Start")
        else:
            self.state.set_state("error", "Nicht startbar")
        self.preview_label.setText(text)
        self.issues.set_issues(issues)
        self.out_label.clear()
        self._ok = ok
        self._update_start()

    def _update_start(self) -> None:
        self.start_button.setEnabled(self._ok and not self._running)
        self.start_button.setToolTip("" if self._ok else "Erst einen passenden Datensatz öffnen "
                                                          "und die Prüfung bestehen.")

    def refresh(self) -> None:
        self._update_preview()

    def set_running(self, running: bool) -> None:
        self._running = running
        if self._left is not None:
            self._left.setEnabled(not running)          # der Lauf rechnet mit den Werten beim Start
        self._update_start()

    # -- Lauf und Ergebnis -------------------------------------------------------------
    def _start(self) -> None:
        self.flush_check()
        if not self._ok or self._running:
            return
        self._job_settings = self.form.settings()
        self._job_disabled = set(self._job_settings.disabled_fields(self.info))    # wirkten nicht
        job = api_jobs.new_job(self.api.WORKFLOW, self._job_settings, self.info.selection,
                               self.output_spec())
        self._reset_result("Läuft …")
        self.main.start_job(self, job, self.title)

    def _reset_result(self, text: str) -> None:
        self.last_result = None
        self.result_state.set_state(None, text)
        self.result_state.setToolTip("")
        self.summary_label.hide()
        for b in (self.open_button, self.more_button):
            b.setEnabled(False)
        self.result_issues.set_issues([])
        self.clear_outputs()

    def _copy_command(self) -> None:
        QApplication.clipboard().setText(command_string(self.last_result.get("command") or []))

    def _run_line(self, result: dict) -> str:
        """Dauer und die vom Standard abweichenden Einstellungen des Laufs."""
        parts = []
        total = (result.get("timings") or {}).get("total_s")
        if total is not None:
            parts.append(f"Dauer {total:.1f} s")
        if self._job_settings is not None:
            metas = self.form.cls.field_meta()
            changed = {k: v for k, v in self._job_settings.non_default().items()
                       if k not in self._job_disabled and k not in self.summary_fields}
            if changed:
                parts.append("abweichend vom Standard: " + ", ".join(
                    f"{metas[k].label} = {_value_text(v, metas[k].unit)}" for k, v in changed.items()))
        return "  ·  ".join(parts)

    def show_result(self, result: dict) -> None:
        self.last_result = result
        out = result.get("output_dir")
        status = result["status"]
        self.result_state.set_state(RESULT_LEVEL.get(status, "error"), STATUS_DE.get(status, status)
                                    + (f"  ·  {Path(out).name}" if out else ""))
        self.result_state.setToolTip(str(out) if out else "")
        lines = [t for t in (self.summary_text(result) if out else "", self._run_line(result)) if t]
        self.summary_label.setText("\n".join(lines))
        self.summary_label.setVisible(bool(lines))
        self.open_button.setEnabled(bool(out))
        self.runjson_action.setEnabled(bool(out) and (Path(out) / "run.json").is_file())
        self.command_action.setEnabled(bool(result.get("command")))
        self.command_action.setToolTip(command_string(result.get("command") or []) or
                                       "Kein Befehl verfügbar")
        self.more_button.setEnabled(self.runjson_action.isEnabled() or self.command_action.isEnabled())
        self.result_issues.set_issues(result.get("issues", []))
        if out:
            self.show_outputs(result, Path(out))
        self._update_preview()                          # naechster Lauf bekaeme _2
