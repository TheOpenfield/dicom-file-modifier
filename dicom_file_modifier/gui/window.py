"""Hauptfenster: Datensatzzeile, Seitenleiste mit Seiten, Fortschritt mit Abbrechen, Protokoll."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from PySide6.QtCore import QSize, Qt, Slot
from PySide6.QtGui import QPalette
from PySide6.QtWidgets import (QDockWidget, QFileDialog, QFrame, QHBoxLayout, QLabel, QListWidget,
                               QListWidgetItem,
                               QMainWindow, QMessageBox, QPlainTextEdit, QProgressBar,
                               QPushButton, QSizePolicy, QStackedWidget, QStyle, QToolBar,
                               QWidget)

from .config import APP_NAME, AppConfig
from .jobs import JobRunner
from .widgets import STATUS_DE

ROOT_LABEL_PX = 300                     # Ergebnis-Stammordner in der Statusleiste (Mitte gekuerzt)
ABOUT_HTML = (
    "<b>{app} {version}</b><br>Strukturanalyse, Dosisindizes und starre Transformation "
    "von DICOM-RT-Daten.<br><br>"
    "<b>Nur für Forschung und Lehre – kein Medizinprodukt.</b> Nicht für die klinische Anwendung "
    "validiert oder zertifiziert; nicht zum Erstellen, Ändern oder Prüfen von Daten für die "
    "Behandlung von Patientinnen und Patienten verwenden. Alle Ergebnisse vor jeder klinischen "
    "Verwendung unabhängig durch qualifizierte Medizinphysik prüfen. Nutzung auf eigene Gefahr."
    "<br><br>Lizenz: MIT")


@dataclass
class Case:
    """Geoeffneter Datensatz-Ordner; Dateien nach der Namenskonvention der CLI (``RS*``/``RD*``/``RP*``, ``CT/``)."""

    folder: Path
    rs: list = field(default_factory=list)
    rd: list = field(default_factory=list)
    rp: list = field(default_factory=list)
    n_ct: int = 0

    @classmethod
    def open(cls, folder) -> "Case":
        folder = Path(folder)
        if not folder.is_dir():
            raise FileNotFoundError(f"Ordner nicht gefunden: {folder}")

        def find(prefix: str) -> list:
            return sorted(folder.glob(f"{prefix}*.dcm"))

        ct = folder / "CT"
        n_ct = len(list(ct.glob("*.dcm"))) if ct.is_dir() else len(find("CT"))
        return cls(folder, find("RS"), find("RD"), find("RP"), n_ct)

    def summary(self) -> str:
        return (f"{self.folder.name}   ·   CT: {self.n_ct} Schichten   ·   RS: {len(self.rs)}"
                f"   ·   RP: {len(self.rp)}   ·   RD: {len(self.rd)}")


class MainWindow(QMainWindow):
    def __init__(self, config: AppConfig):
        super().__init__()
        from .dose_page import DosePage
        from .structures_page import StructuresPage

        self.config = config
        self.case = None
        self.runner = JobRunner(config.jobs_dir, self)
        self._job_page = None
        self.setWindowTitle(APP_NAME)
        self.setAcceptDrops(True)
        self.resize(1320, 860)

        bar = QToolBar("Datensatz", self)
        bar.setMovable(False)
        self.addToolBar(bar)
        self.open_case_button = QPushButton(
            self.style().standardIcon(QStyle.StandardPixmap.SP_DirOpenIcon), "Datensatz öffnen …")
        self.open_case_button.clicked.connect(self.choose_case)
        bar.addWidget(self.open_case_button)
        self.case_label = QLabel("  Kein Datensatz geöffnet: Ordner wählen oder hierher ziehen")
        bar.addWidget(self.case_label)
        spacer = QWidget()
        spacer.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        bar.addWidget(spacer)
        self.root_action = bar.addAction("Ergebnisordner …", self.choose_results_root)

        self.structures = StructuresPage(self)
        self.dose = DosePage(self)
        self.pages = [self.structures, self.dose]
        self.nav = QListWidget()
        self.nav.setFixedWidth(170)
        self.nav.setFrameShape(QFrame.Shape.NoFrame)
        palette = self.nav.palette()
        palette.setColor(QPalette.ColorRole.Base, palette.color(QPalette.ColorRole.Window))
        self.nav.setPalette(palette)
        self.stack = QStackedWidget()
        for page in self.pages:
            item = QListWidgetItem(page.title)
            item.setSizeHint(QSize(0, 34))
            self.nav.addItem(item)
            self.stack.addWidget(page)
        self.nav.currentRowChanged.connect(self.stack.setCurrentIndex)
        self.nav.setCurrentRow(0)
        central = QWidget()
        row = QHBoxLayout(central)
        row.addWidget(self.nav)
        row.addWidget(self.stack, 1)
        self.setCentralWidget(central)

        self.job_label = QLabel("Bereit")
        self.progress = QProgressBar()
        self.progress.setMaximumWidth(280)
        self.progress.hide()
        self.cancel_button = QPushButton("Abbrechen")
        self.cancel_button.hide()
        self.cancel_button.clicked.connect(self.cancel_job)
        self.root_label = QLabel()
        self.root_label.setEnabled(False)                 # grau: nur Information
        self.statusBar().addWidget(self.job_label, 1)
        self.statusBar().addPermanentWidget(self.root_label)
        self.statusBar().addPermanentWidget(self.progress)
        self.statusBar().addPermanentWidget(self.cancel_button)

        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumBlockCount(5000)
        dock = self.log_dock = QDockWidget("Protokoll", self)
        dock.setObjectName("log")
        dock.setWidget(self.log)
        self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, dock)
        dock.hide()
        bar.addAction(dock.toggleViewAction())
        bar.addAction("Info", self.show_about)

        self.runner.event.connect(self._on_event)
        self.runner.finished.connect(self._on_finished)
        self._update_root_tooltip()

    # -- Datensatz ----------------------------------------------------------------------------
    def show_about(self) -> None:
        from .. import __version__
        QMessageBox.about(self, f"Über {APP_NAME}", ABOUT_HTML.format(app=APP_NAME, version=__version__))

    def choose_case(self) -> None:
        start = str(self.case.folder.parent) if self.case else ""
        folder = QFileDialog.getExistingDirectory(self, "Datensatz-Ordner wählen", start)
        if folder:
            self.open_case(folder)

    def open_case(self, folder) -> None:
        try:
            case = Case.open(folder)
        except OSError as exc:
            QMessageBox.warning(self, APP_NAME, str(exc))
            return
        self.case = case
        self.case_label.setText("  " + case.summary())
        self.setWindowTitle(f"{APP_NAME} – {case.folder.name}")
        for page in self.pages:
            page.set_case(case)

    def dragEnterEvent(self, event) -> None:
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event) -> None:
        paths = [Path(u.toLocalFile()) for u in event.mimeData().urls() if u.isLocalFile()]
        if paths:
            self.open_case(paths[0] if paths[0].is_dir() else paths[0].parent)

    # -- Ergebnisordner ------------------------------------------------------------------
    def choose_results_root(self) -> None:
        folder = QFileDialog.getExistingDirectory(self, "Ergebnis-Stammordner wählen",
                                                  str(self.config.results_root))
        if folder:
            self.config.results_root = Path(folder)
            self.config.save()
            self._update_root_tooltip()
            for page in self.pages:
                page.refresh()

    def _update_root_tooltip(self) -> None:
        root = str(self.config.results_root)
        self.root_action.setToolTip(f"Ergebnisse unter: {root}")
        short = self.root_label.fontMetrics().elidedText(root, Qt.TextElideMode.ElideMiddle, ROOT_LABEL_PX)
        self.root_label.setText(f"Ergebnisse: {short}")
        self.root_label.setToolTip(root)

    # -- Jobs ----------------------------------------------------------------------------
    def start_job(self, page, job, title: str) -> None:
        self._job_page = page
        self.log.appendPlainText(f"===== {title} =====")
        self.runner.start(job)
        self.job_label.setText(f"{title}: startet …")
        self.progress.setRange(0, 0)
        self.progress.show()
        self.cancel_button.setEnabled(True)
        self.cancel_button.show()
        for p in self.pages:
            p.set_running(True)

    def cancel_job(self) -> None:
        self.cancel_button.setEnabled(False)
        self.job_label.setText("Wird abgebrochen …")
        self.runner.cancel()

    @Slot(dict)
    def _on_event(self, ev: dict) -> None:
        kind = ev.get("type")
        if kind == "log":
            msg = ev.get("msg", "")
            if self.runner.job is not None:              # der Staging-Ordner existiert nach dem Lauf nicht mehr
                msg = msg.replace(str(self.runner.job.staging_dir), "<Ergebnisordner>")
            self.log.appendPlainText(msg)
        elif kind == "issue":
            self.log.appendPlainText(f"[{ev.get('level')}] {ev.get('message_de')}")
        elif kind == "stage" and self.cancel_button.isEnabled():
            step = f"Schritt {ev['i']}/{ev['n']}: " if ev.get("i") and ev.get("n") else ""
            self.job_label.setText(step + str(ev.get("label", "")))
            self.progress.setRange(0, 0)
        elif kind == "progress" and (ev.get("total") or 0) > 0:
            self.progress.setRange(0, 1000)
            self.progress.setValue(int(1000 * ev["done"] / ev["total"]))

    @Slot(dict)
    def _on_finished(self, result: dict) -> None:
        self.progress.hide()
        self.cancel_button.hide()
        self.job_label.setText(STATUS_DE.get(result.get("status"), str(result.get("status"))))
        if result.get("status") == "failed":
            self.log_dock.show()
        for p in self.pages:
            p.set_running(False)
        page, self._job_page = self._job_page, None
        if page is not None:
            page.show_result(result)

    def closeEvent(self, event) -> None:
        if self.runner.running:
            answer = QMessageBox.question(self, APP_NAME, "Ein Lauf ist noch aktiv. Abbrechen und beenden?")
            if answer != QMessageBox.StandardButton.Yes:
                event.ignore()
                return
            self.runner.cancel()
            self.runner.kill()
        super().closeEvent(event)
