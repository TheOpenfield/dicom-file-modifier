"""
Grundgeruest der Workflow-Seiten: Eingaben -> Einstellungen -> Pruefung ->
Start -> Ergebnis.  Die Seite liefert ``api`` (Workflow-Modul), ``selection``,
``on_inspected``, ``check``, ``clear_outputs``/``show_outputs`` und optional
``summary_text``; Inspektion im Hintergrund, Pruefung bei jeder Aenderung,
Lauf ueber das Hauptfenster (Worker).
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (QApplication, QGroupBox, QHBoxLayout, QLabel, QPushButton,
                               QScrollArea, QSplitter, QVBoxLayout, QWidget)

from ..api import jobs as api_jobs
from ..api.fields import command_string
from ..api.issues import Issue, has_errors, issue_from_exception
from ..api.outputs import OutputSpec
from .jobs import run_in_background
from .widgets import STATUS_DE, IssueList, StateLine, breakable, open_path

# Laengster Dateiname unter dem Ergebnisordner (z.B. RS.<64-Zeichen-UID>_analysis.json)
# plus Trenner: laengere Ordnerpfade stossen an MAX_PATH (260) von Windows
LONGEST_FILE_NAME = 100
MAX_PATH = 259


def _value_text(v) -> str:
    if isinstance(v, bool):
        return "ja" if v else "nein"
    if isinstance(v, (list, tuple)):
        return ", ".join(map(str, v))
    return "aus" if v is None else str(v)


def names_text(names: list, limit: int = 8) -> str:
    if not names:
        return "keine"
    more = f" … (+{len(names) - limit})" if len(names) > limit else ""
    return ", ".join(names[:limit]) + more


class WorkflowPage(QWidget):
    title = ""
    start_text = "Starten"
    api = None                          # Workflow-Modul: api.structures, api.dose, api.transform

    def __init__(self, main):
        super().__init__()
        self.main = main                # MainWindow
        self.case = self.info = self.last_result = self.form = self._left = None
        self._token = 0
        self._running = self._ok = False
        self._job_settings = None
        self.state = StateLine()
        self.preview_label = QLabel("Kein Datensatz geöffnet.")
        self.preview_label.setWordWrap(True)
        self.out_label = QLabel()
        self.out_label.setWordWrap(True)
        self.issues = IssueList()
        self.issues.hide()
        self.start_button = QPushButton(self.start_text)
        self.start_button.setEnabled(False)
        self.start_button.setMinimumHeight(34)
        font = self.start_button.font()
        font.setBold(True)
        self.start_button.setFont(font)
        self.start_button.clicked.connect(self._start)

        self.status_label = QLabel("Noch kein Ergebnis.")
        self.status_label.setWordWrap(True)
        self.summary_label = QLabel()
        self.summary_label.setWordWrap(True)
        self.summary_label.hide()
        self.open_button = QPushButton("Ordner öffnen")
        self.open_button.clicked.connect(lambda: open_path(self.last_result["output_dir"]))
        self.command_button = QPushButton("Befehl kopieren")
        self.command_button.setToolTip("Gleichwertigen dfm-Befehl in die Zwischenablage kopieren")
        self.command_button.clicked.connect(self._copy_command)
        self.runjson_button = QPushButton("run.json")
        self.runjson_button.setToolTip("Protokoll des Laufs: Versionen, Einstellungen, Eingaben mit UIDs")
        self.runjson_button.clicked.connect(lambda: open_path(Path(self.last_result["output_dir"]) / "run.json"))
        for b in (self.open_button, self.command_button, self.runjson_button):
            b.setEnabled(False)
        self.result_issues = IssueList()
        self.result_issues.hide()

    def build(self, inputs: list, form, results: list, buttons: list = ()) -> None:
        """Layout: links (scrollbar) Eingaben, Einstellungen, Pruefung, Start; rechts das Ergebnis."""
        self.form = form
        form.changed.connect(self._update_preview)
        left = QWidget()
        col = QVBoxLayout(left)
        for widget, stretch in inputs:
            col.addWidget(widget, stretch)
        box = QGroupBox("Einstellungen")
        QVBoxLayout(box).addWidget(form)
        col.addWidget(box)
        check = QGroupBox("Prüfung vor dem Start")
        lay = QVBoxLayout(check)
        for w in (self.state, self.preview_label, self.issues, self.out_label):
            lay.addWidget(w)
        col.addWidget(check)
        col.addWidget(self.start_button)
        if not any(stretch for _, stretch in inputs):
            col.addStretch(1)                           # Gruppen oben kompakt halten
        self._left = left
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(left)

        result = QGroupBox("Ergebnis")
        lay = QVBoxLayout(result)
        head = QHBoxLayout()
        text = QVBoxLayout()
        text.addWidget(self.status_label)
        text.addWidget(self.summary_label)
        head.addLayout(text, 1)
        for b in list(buttons) + [self.command_button, self.runjson_button, self.open_button]:
            head.addWidget(b)
        lay.addLayout(head)
        lay.addWidget(self.result_issues)
        for widget, stretch in results:
            lay.addWidget(widget, stretch)

        outer = QSplitter()
        outer.addWidget(scroll)
        outer.addWidget(result)
        outer.setStretchFactor(1, 1)
        outer.setSizes([580, 740])
        QVBoxLayout(self).addWidget(outer)

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

    def _update_preview(self) -> None:
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
        self._job_settings = self.form.settings()
        job = api_jobs.new_job(self.api.WORKFLOW, self._job_settings, self.info.selection,
                               self.output_spec())
        self.last_result = None
        self.status_label.setText("Läuft …")
        self.summary_label.hide()
        for b in (self.open_button, self.command_button, self.runjson_button):
            b.setEnabled(False)
        self.result_issues.set_issues([])
        self.clear_outputs()
        self.main.start_job(self, job, self.title)

    def _copy_command(self) -> None:
        QApplication.clipboard().setText(command_string(self.last_result.get("command") or []))

    def _run_line(self, result: dict) -> str:
        """Dauer und die vom Standard abweichenden Einstellungen des Laufs."""
        parts = []
        total = (result.get("timings") or {}).get("total_s")
        if total is not None:
            parts.append(f"Dauer {total:.1f} s")
        if self._job_settings is not None:
            changed = self._job_settings.non_default()
            if changed:
                parts.append("abweichend vom Standard: " + ", ".join(
                    f"{self.form.label_of(k)} = {_value_text(v)}" for k, v in changed.items()))
        return "  ·  ".join(parts)

    def show_result(self, result: dict) -> None:
        self.last_result = result
        out = result.get("output_dir")
        self.status_label.setText(STATUS_DE.get(result["status"], result["status"])
                                  + (f":  {breakable(out)}" if out else ""))
        lines = [t for t in (self.summary_text(result) if out else "", self._run_line(result)) if t]
        self.summary_label.setText("\n".join(lines))
        self.summary_label.setVisible(bool(lines))
        self.open_button.setEnabled(bool(out))
        self.runjson_button.setEnabled(bool(out) and (Path(out) / "run.json").is_file())
        self.command_button.setEnabled(bool(result.get("command")))
        self.command_button.setToolTip(command_string(result.get("command") or []) or
                                       "Kein Befehl verfügbar")
        self.result_issues.set_issues(result.get("issues", []))
        if out:
            self.show_outputs(result, Path(out))
        self._update_preview()                          # naechster Lauf bekaeme _2
