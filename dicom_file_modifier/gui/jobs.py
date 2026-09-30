"""
Jobs der App.  ``JobRunner`` startet ``dfm worker --job`` als Kindprozess und
meldet dessen JSON-Zeilen als Qt-Signale; ``run_in_background`` fuehrt kurze
Lesearbeiten (Inspektion) im Thread-Pool aus.

``subprocess`` statt QProcess: Der Worker braucht ``CREATE_NO_WINDOW``, sonst
oeffnet die fensterlose App je Job ein Konsolenfenster.  Abbruch: Datei
``<job>.cancel``, nach ``KILL_AFTER_MS`` hart beenden; der Staging-Ordner
wird danach immer geloescht.
"""

from __future__ import annotations

import collections
import json
import subprocess
import sys
import threading
from pathlib import Path
from typing import Callable

from PySide6.QtCore import QObject, QThreadPool, QTimer, Signal, Slot

from ..api import jobs as api_jobs
from ..api.issues import Issue
from ..api.results import JobResult
from ..api.worker import cancel_file

KILL_AFTER_MS = 5000


def worker_command(job_path) -> list:
    """``dfm.exe`` neben der App-EXE (gefroren), sonst dieses Python (python.exe statt pythonw.exe)."""
    if getattr(sys, "frozen", False):
        return [str(Path(sys.executable).with_name("dfm.exe")), "worker", "--job", str(job_path)]
    exe = Path(sys.executable)
    if exe.name.lower() == "pythonw.exe" and exe.with_name("python.exe").exists():
        exe = exe.with_name("python.exe")
    return [str(exe), "-m", "dicom_file_modifier", "worker", "--job", str(job_path)]


class JobRunner(QObject):
    """Ein Job zur Zeit; ``event`` je Protokollzeile, ``finished`` genau einmal (``JobResult`` als dict)."""

    event = Signal(dict)
    finished = Signal(dict)
    _line = Signal(dict)
    _exited = Signal(int)

    def __init__(self, jobs_dir, parent=None):
        super().__init__(parent)
        self.jobs_dir = Path(jobs_dir)
        self.job = self.job_path = self._proc = self._result = None
        self._cancel_requested = False
        self._stderr: collections.deque = collections.deque(maxlen=200)
        self._line.connect(self._on_line)
        self._exited.connect(self._on_exit)

    @property
    def running(self) -> bool:
        return self._proc is not None

    def start(self, job) -> None:
        if self.running:
            raise RuntimeError("Es laeuft bereits ein Job.")
        self.job, self._result, self._cancel_requested = job, None, False
        self._stderr.clear()
        self.job_path = job.save(self.jobs_dir / f"{job.job_id}.json")
        self._proc = subprocess.Popen(
            worker_command(self.job_path), stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, text=True, encoding="utf-8", errors="replace",
            creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
        err = threading.Thread(target=self._read_stderr, args=(self._proc,), daemon=True)
        err.start()
        threading.Thread(target=self._read_stdout, args=(self._proc, err, job, self.job_path),
                         daemon=True).start()

    def cancel(self) -> None:
        if not self.running or self._cancel_requested:
            return
        self._cancel_requested = True
        cancel_file(self.job_path).touch()
        QTimer.singleShot(KILL_AFTER_MS, self.kill)

    def kill(self) -> None:
        """Worker sofort hart beenden (der Staging-Ordner wird danach geloescht)."""
        if self._proc is not None and self._proc.poll() is None:
            self._proc.kill()

    # -- Lese-Threads ------------------------------------------------------------------
    def _read_stderr(self, proc) -> None:
        for raw in proc.stderr:
            self._stderr.append(raw.rstrip("\n"))

    def _read_stdout(self, proc, err: threading.Thread, job, job_path: Path) -> None:
        for raw in proc.stdout:
            try:
                ev = json.loads(raw)
            except ValueError:
                continue
            if isinstance(ev, dict):
                self._line.emit(ev)
        code = proc.wait()
        err.join(timeout=5)
        api_jobs.remove_tree(job.staging_dir)            # nach Kill oder Absturz; sonst schon weg
        for p in (job_path, cancel_file(job_path)):
            try:
                p.unlink(missing_ok=True)
            except OSError:
                pass
        self._exited.emit(code)

    # -- GUI-Thread --------------------------------------------------------------------
    @Slot(dict)
    def _on_line(self, ev: dict) -> None:
        if ev.get("type") == "result":
            self._result = ev.get("result")
        self.event.emit(ev)

    @Slot(int)
    def _on_exit(self, code: int) -> None:
        self._proc = None
        result = self._result
        if result is None:                               # ohne Ergebnis beendet: Kill oder Absturz
            if self._cancel_requested:
                res = JobResult(self.job.workflow, "cancelled", 3,
                                issues=[Issue("info", "JOB.CANCELLED", "Der Lauf wurde abgebrochen.")])
            else:
                res = JobResult(self.job.workflow, "failed", code or 1, issues=[Issue(
                    "error", "JOB.WORKER_DIED",
                    f"Der Rechenprozess wurde unerwartet beendet (Exit-Code {code}).",
                    hint_de="Details im Protokoll.", detail="\n".join(self._stderr))])
            result = res.to_dict()
        self.finished.emit(result)


class _Relay(QObject):
    done = Signal(object, object)

    def __init__(self, on_done: Callable):
        super().__init__()
        self.on_done = on_done
        self.done.connect(self._finish)

    @Slot(object, object)
    def _finish(self, value, error) -> None:
        _ACTIVE.discard(self)
        self.on_done(value, error)


_ACTIVE: set = set()


def run_in_background(fn: Callable, on_done: Callable) -> None:
    """``fn()`` im Thread-Pool; danach ``on_done(wert, fehler)`` im GUI-Thread."""
    relay = _Relay(on_done)
    _ACTIVE.add(relay)

    def task():
        try:
            value, error = fn(), None
        except Exception as exc:  # noqa: BLE001 - an den Aufrufer
            value, error = None, exc
        relay.done.emit(value, error)

    QThreadPool.globalInstance().start(task)
