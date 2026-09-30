"""
Worker der Desktop-App: ``dfm worker --job JOB.json`` fuehrt genau einen Job
aus (``jobs.execute_job``) und meldet sich per JSON-Lines auf stdout.

Worker -> App, eine JSON-Zeile je Ereignis (Feld ``type``):

  ``hello``      ``{v, pid, versions, threads}``; immer die erste Zeile
  ``stage``      ``{key, label, i, n}``; ``i``/``n`` aus ``planned_stages`` (sonst null)
  ``progress``   ``{done, total, text}``; hoechstens etwa 10 je Sekunde
  ``log``        ``{level, msg}``; jede Ausgabezeile des Kerns
  ``issue``      Felder eines ``Issue``, fuer Befunde vor dem Lauf (die uebrigen
                 stehen im Ergebnis)
  ``artifact``   ``{role, path}``; Ergebnisdateien, absolut, nach dem Commit
  ``heartbeat``  ``{t}``; alle 2 s
  ``result``     ``{status, result}``; die letzte Zeile, ``result`` ist ``JobResult.to_dict()``

App -> Worker: Abbruch ueber die Datei ``<job>.cancel`` neben der Job-Datei
(``cancel_file``); der Lauf endet am naechsten Checkpoint.  stdin liest der
Worker bewusst nicht: Unter Windows blockiert ein wartender Lese-Thread auf
einer Pipe jedes ``GetFileType`` auf stdin, und das rufen DLLs beim Laden
auf (numpy/OpenBLAS) - der Worker hinge schon beim ersten Import, sogar
``os._exit`` (Loader-Lock).

Exit-Codes: 0 ok, 2 Eingabefehler, 3 abgebrochen, 1 interner Fehler
(Traceback und faulthandler auf stderr).

stdout-Schutz: Das Protokoll laeuft ueber ein Duplikat des urspruenglichen
stdout.  fd 1 zeigt danach auf stderr, ``sys.stdout`` ist ein
``LineForwarder``: ``print`` im Kern wird zu ``log``, und Ausgaben von
C-Bibliotheken koennen das Protokoll nicht zerstoeren.  Vor dem ersten
numpy-Import begrenzt der Worker die BLAS-/OpenMP-Threads, weil jeder
Thread Speicher reserviert (ca. 1,5 GB Commit je Prozess bei 24 Kernen).
"""

from __future__ import annotations

import argparse
import faulthandler
import json
import os
import sys
import threading
import time
import traceback
from pathlib import Path
from typing import Optional

from .. import _runtime

PROTOCOL_VERSION = 1
HEARTBEAT_S = 2.0
PROGRESS_INTERVAL_S = 0.1
CANCEL_POLL_S = 0.2
THREAD_VARS = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS")


def cancel_file(job_path) -> Path:
    """Abbruch-Signal eines Jobs: ``<job>.cancel`` neben der Job-Datei (die App legt sie an)."""
    return Path(job_path).with_suffix(".cancel")


class Emitter:
    """Schreibt Ereignisse als JSON-Zeilen (thread-sicher, sofort geleert)."""

    def __init__(self, stream):
        self.stream = stream
        self._lock = threading.Lock()

    def __call__(self, type_: str, **fields) -> None:
        from .results import jsonable
        line = json.dumps({"type": type_, **jsonable(fields)}, ensure_ascii=False, allow_nan=False)
        with self._lock:
            self.stream.write(line + "\n")
            self.stream.flush()


class WorkerContext(_runtime.JobContext):
    """``JobContext`` des Workers: Stufen und Fortschritt als Ereignisse, Abbruch per Datei."""

    def __init__(self, emit: Emitter, cancel: Optional[Path] = None, planned: Optional[list] = None):
        self.emit, self.cancel, self.planned = emit, cancel, list(planned or [])
        self._last_progress = self._last_poll = 0.0
        self._cancelled = False

    def stage(self, key: str, label: str) -> None:
        i = self.planned.index(key) + 1 if key in self.planned else None
        self.emit("stage", key=key, label=label, i=i, n=len(self.planned) or None)

    def progress(self, done: float, total: float, text: str = "") -> None:
        now = time.monotonic()
        if done < total and now - self._last_progress < PROGRESS_INTERVAL_S:
            return
        self._last_progress = now
        self.emit("progress", done=done, total=total, text=text)

    def log(self, message: str, level: str = "info") -> None:
        """Nichts: Protokollzeilen kommen ueber den ``LogSink`` (ein Ereignis je Zeile)."""

    def check_cancel(self) -> None:
        now = time.monotonic()
        if not self._cancelled and self.cancel is not None and now - self._last_poll >= CANCEL_POLL_S:
            self._last_poll = now
            self._cancelled = self.cancel.exists()
        if self._cancelled:
            raise _runtime.JobCancelled()


def limit_threads(n: Optional[int] = None) -> None:
    """BLAS-/OpenMP-Threads begrenzen, sofern nicht schon gesetzt (vor dem numpy-Import)."""
    n = n or min(4, os.cpu_count() or 1)
    for var in THREAD_VARS:
        os.environ.setdefault(var, str(n))


def protect_stdout():
    """Protokoll-Stream auf einem Duplikat von fd 1; fd 1 danach auf stderr umgelenkt."""
    try:
        sys.stdout.flush()
    except (AttributeError, ValueError):
        pass
    proto = os.fdopen(os.dup(1), "w", encoding="utf-8", newline="\n")
    os.dup2(2, 1)
    return proto


def _heartbeat(emit: Emitter, stop: threading.Event) -> None:
    while not stop.wait(HEARTBEAT_S):
        emit("heartbeat", t=round(time.time(), 3))


def _planned(job) -> list:
    from . import workflow
    from .selection import CaseSelection
    try:
        wf = workflow(job.workflow)
        return wf.planned_stages(wf.Settings.from_dict(job.settings),
                                 CaseSelection.from_dict(job.selection))
    except Exception:  # noqa: BLE001 - der Lauf meldet den Fehler selbst
        return []


def run_worker(job_path: str, proto) -> int:
    """Job ausfuehren und protokollieren; ``proto`` ist der Protokoll-Stream."""
    from ..cli import versions
    from .issues import exit_code_for, issue_from_exception
    from .jobs import Job, LineForwarder, LogSink, execute_job
    from .results import JobResult

    emit = Emitter(proto)
    emit("hello", v=PROTOCOL_VERSION, pid=os.getpid(), versions=versions(),
         threads={var: os.environ.get(var) for var in THREAD_VARS})
    stop = threading.Event()
    beat = threading.Thread(target=_heartbeat, args=(emit, stop), daemon=True)
    beat.start()
    try:
        try:
            job = Job.load(job_path)
        except Exception as exc:  # noqa: BLE001 - kaputte oder fehlende Job-Datei
            result = JobResult("?", "failed", exit_code_for(exc), issues=[issue_from_exception(exc)])
        else:
            sink = LogSink(on_line=lambda level, text: level == "stage" or emit("log", level=level, msg=text))
            saved, sys.stdout = sys.stdout, LineForwarder(sink)
            try:
                with _runtime.use(WorkerContext(emit, cancel_file(job_path), _planned(job))):
                    result = execute_job(job, sink, on_issue=lambda i: emit("issue", **i.to_dict()))
            except Exception as exc:  # noqa: BLE001 - Fehler im Job-Rahmen selbst
                traceback.print_exc(file=sys.stderr)
                result = JobResult(job.workflow, "failed", 1, issues=[issue_from_exception(exc)])
            finally:
                sys.stdout.flush()
                sys.stdout = saved
        if result.output_dir:
            for role, rel in result.outputs.items():
                emit("artifact", role=role, path=str(Path(result.output_dir) / rel))
        stop.set()
        beat.join(timeout=5)                 # "result" ist garantiert die letzte Zeile
        emit("result", status=result.status, result=result.to_dict())
        return result.exit_code
    finally:
        stop.set()


def main(argv: Optional[list] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="dfm worker",
        description="Einen Job der Desktop-App ausfuehren; Protokoll als JSON-Zeilen auf stdout, "
                    "die Datei <job>.cancel bricht ab.")
    parser.add_argument("--job", required=True, metavar="JOB.json", help="Job-Datei")
    args = parser.parse_args(argv)
    limit_threads()
    faulthandler.enable(file=sys.stderr)
    proto = protect_stdout()
    try:
        return run_worker(args.job, proto)
    finally:
        proto.close()
