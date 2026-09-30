"""
Jobs: ein vollstaendig beschriebener Lauf als JSON-Datei, ausgefuehrt vom
Worker (``dfm worker --job``) oder im laufenden Prozess (``run_job_inprocess``).

Ablauf von ``execute_job``:

1. Pruefen (Job-Datei, Ausgabe, Einstellungen, Kopfzeile der Sammel-CSV);
   bei Fehlern ``failed`` (Exit 2), ohne etwas zu schreiben.
2. Im Staging-Ordner ``<root>/.stg<job_id>`` rechnen (``workflow.run``);
   ``run.log`` (alle Ausgaben) und ``run.json`` (Versionen, Einstellungen,
   Eingaben mit UIDs, Ergebnis, aequivalenter Befehl) liegen mit darin.
3. Bei Erfolg den Staging-Ordner in den Ergebnisordner umbenennen, mit
   Wiederholung, falls ein Virenscanner Dateien sperrt; ``overwrite`` tauscht
   einen vorhandenen Ordner erst nach dem Umbenennen aus.
4. Erst danach: Zeilen an die Sammel-CSV anhaengen (Workflow-Hook).

Bei Fehler oder Abbruch wird der Staging-Ordner geloescht; bleibt einer nach
einem harten Abbruch liegen, raeumt ``cleanup_staging`` ihn beim naechsten
Start weg.  Das Modul importiert die Workflows erst bei Bedarf.
"""

from __future__ import annotations

import datetime as _dt
import io
import json
import os
import shutil
import stat
import threading
import time
import uuid
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable, Optional

from .. import _runtime
from . import workflow as _workflow
from .fields import command_string
from .issues import Issue, exit_code_for, has_errors, issue_from_exception
from .outputs import OutputSpec
from .results import RUN_JSON, RUN_LOG, JobResult, jsonable
from .selection import CaseSelection

JOB_VERSION = 1
STAGING_PREFIX = ".stg"
OLD_PREFIX = ".old"


def _now() -> str:
    return _dt.datetime.now().isoformat(timespec="seconds")


@dataclass
class Job:
    job_id: str
    workflow: str
    settings: dict
    selection: dict
    output: dict
    staging: Optional[str] = None       # Default: <root>/.stg<job_id>
    v: int = JOB_VERSION

    @property
    def staging_dir(self) -> Path:
        if self.staging:
            return Path(self.staging)
        return Path(self.output["root"]) / f"{STAGING_PREFIX}{self.job_id}"

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "Job":
        if d.get("v", JOB_VERSION) != JOB_VERSION:
            raise ValueError(f"Job-Version {d.get('v')!r} wird nicht unterstuetzt (erwartet {JOB_VERSION}).")
        unknown = sorted(set(d) - set(cls.__dataclass_fields__))
        if unknown:
            raise ValueError(f"Unbekannte Felder in der Job-Datei: {', '.join(unknown)}")
        return cls(**d)

    @classmethod
    def load(cls, path) -> "Job":
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))

    def save(self, path) -> Path:
        """Atomar schreiben (``.tmp`` + ``os.replace``)."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_name(path.name + ".tmp")
        tmp.write_text(json.dumps(self.to_dict(), indent=2, ensure_ascii=False), encoding="utf-8")
        os.replace(tmp, path)
        return path


def new_job(workflow: str, settings, selection: CaseSelection, output: OutputSpec,
            job_id: Optional[str] = None) -> Job:
    return Job(job_id=job_id or uuid.uuid4().hex[:12], workflow=workflow,
               settings=settings.to_dict(), selection=selection.to_dict(), output=output.to_dict())


# ---------------------------------------------------------------------------
# Dateisystem: Wiederholen, Loeschen, Umbenennen
# ---------------------------------------------------------------------------

def _retry(fn: Callable, attempts: int = 7, delay: float = 0.1):
    """``fn()`` mit Wiederholung bei ``OSError`` (gesperrte Dateien unter Windows)."""
    for i in range(attempts):
        try:
            return fn()
        except OSError:
            if i == attempts - 1:
                raise
            time.sleep(delay * 2 ** i)


def _clear_readonly(func, path, _exc):
    os.chmod(path, stat.S_IWRITE)
    func(path)


def remove_tree(path, attempts: int = 7) -> bool:
    """Ordner loeschen, mit Wiederholung; ``True``, wenn er danach nicht mehr existiert."""
    path = Path(path)
    if not path.exists():
        return True
    try:
        _retry(lambda: shutil.rmtree(path, onexc=_clear_readonly), attempts)
    except OSError:
        return False
    return not path.exists()


def cleanup_staging(root, keep: tuple = ()) -> list:
    """Verwaiste ``.stg*``/``.old*``-Ordner unter ``root`` loeschen (App-Start); die geloeschten."""
    root = Path(root)
    if not root.is_dir():
        return []
    keep = {Path(k).resolve() for k in keep}
    removed = []
    for p in root.iterdir():
        if (p.is_dir() and p.name.startswith((STAGING_PREFIX, OLD_PREFIX))
                and p.resolve() not in keep and remove_tree(p)):
            removed.append(p)
    return removed


def _free_name(final: Path) -> Path:
    n = 2
    while final.with_name(f"{final.name}_{n}").exists():
        n += 1
    return final.with_name(f"{final.name}_{n}")


def commit_staging(staging: Path, final: Path, policy: str = "suffix") -> Path:
    """
    Staging-Ordner zum Ergebnisordner machen und dessen Pfad liefern.  Existiert
    ``final`` inzwischen: ``suffix`` nimmt den naechsten freien Namen,
    ``overwrite`` benennt den alten Ordner zuerst um und loescht ihn erst,
    wenn das neue Ergebnis an seinem Platz ist.
    """
    staging, final = Path(staging), Path(final)
    final.parent.mkdir(parents=True, exist_ok=True)
    old = None
    if final.exists():
        if policy == "overwrite":
            old = final.with_name(f"{OLD_PREFIX}{uuid.uuid4().hex[:8]}_{final.name}")
            _retry(lambda: os.replace(final, old))
        else:
            final = _free_name(final)
    try:
        _retry(lambda: os.replace(staging, final))
    except OSError:
        if old is not None:
            os.replace(old, final)             # altes Ergebnis zurueck
        raise
    if old is not None:
        remove_tree(old)
    return final


# ---------------------------------------------------------------------------
# Protokoll: Zeilen von print() und JobContext.log nach run.log (und zum Worker)
# ---------------------------------------------------------------------------

class LogSink:
    """Sammelt Protokollzeilen: ``run.log`` (falls offen) und optional ``on_line(level, text)``."""

    def __init__(self, on_line: Optional[Callable] = None):
        self.on_line = on_line
        self.file = None
        self._lock = threading.Lock()

    def line(self, text: str, level: str = "info") -> None:
        with self._lock:
            if self.file is not None:
                self.file.write(text + "\n")
                self.file.flush()
        if self.on_line is not None:
            self.on_line(level, text)

    def open(self, path: Path) -> None:
        self.close()
        self.file = open(path, "w", encoding="utf-8")

    def close(self) -> None:
        with self._lock:
            if self.file is not None:
                self.file.close()
                self.file = None


class LineForwarder(io.TextIOBase):
    """Text-Stream fuer ``sys.stdout``: jede vollstaendige Zeile geht an ``LogSink``."""

    def __init__(self, sink: LogSink):
        self.sink = sink
        self._buf = ""
        self._lock = threading.Lock()

    def writable(self) -> bool:
        return True

    def write(self, s: str) -> int:
        with self._lock:
            self._buf += s.replace("\r\n", "\n").replace("\r", "\n")
            *lines, self._buf = self._buf.split("\n")
        for ln in lines:
            self.sink.line(ln)
        return len(s)

    def flush(self) -> None:
        with self._lock:
            rest, self._buf = self._buf, ""
        if rest:
            self.sink.line(rest)

    @property
    def encoding(self) -> str:        # manche Bibliotheken fragen danach
        return "utf-8"


class _LoggingContext(_runtime.JobContext):
    """Leitet an ``parent`` weiter und schreibt Stufen und ``log`` zusaetzlich ins Protokoll."""

    def __init__(self, parent: _runtime.JobContext, sink: LogSink):
        self.parent, self.sink = parent, sink

    def stage(self, key: str, label: str) -> None:
        self.sink.line(f"== {label}", "stage")
        self.parent.stage(key, label)

    def progress(self, done: float, total: float, text: str = "") -> None:
        self.parent.progress(done, total, text)

    def log(self, message: str, level: str = "info") -> None:
        self.sink.line(message, level)
        self.parent.log(message, level)

    def check_cancel(self) -> None:
        self.parent.check_cancel()


# ---------------------------------------------------------------------------
# run.json
# ---------------------------------------------------------------------------

def describe_inputs(sel: CaseSelection) -> list:
    """Eingaben mit UIDs fuer ``run.json`` (nur Header gelesen)."""
    import pydicom

    out = []
    for role in ("rs", "rd", "rp"):
        path = getattr(sel, role)
        if not path:
            continue
        entry = {"role": role, "path": str(path)}
        try:
            ds = pydicom.dcmread(str(path), stop_before_pixels=True,
                                 specific_tags=["SOPInstanceUID", "Modality"])
            entry.update(modality=str(ds.get("Modality", "")),
                         sop_instance_uid=str(ds.get("SOPInstanceUID", "")))
        except Exception as exc:  # noqa: BLE001 - nur Beschreibung
            entry["error"] = f"{type(exc).__name__}: {exc}"
        out.append(entry)
    if sel.eclipse_ref:
        out.append({"role": "eclipse_ref", "path": str(sel.eclipse_ref)})
    if sel.ct_files:
        entry = {"role": "ct", "n_files": len(sel.ct_files),
                 "dir": str(sel.ct_dir) if sel.ct_dir else None}
        try:
            ds = pydicom.dcmread(str(sel.ct_files[0]), stop_before_pixels=True,
                                 specific_tags=["SeriesInstanceUID", "FrameOfReferenceUID"])
            entry.update(series_instance_uid=str(ds.get("SeriesInstanceUID", "")),
                         frame_of_reference_uid=str(ds.get("FrameOfReferenceUID", "")))
        except Exception as exc:  # noqa: BLE001
            entry["error"] = f"{type(exc).__name__}: {exc}"
        out.append(entry)
    return out


def write_run_json(folder: Path, job: Job, result: JobResult, started: str, extra: dict) -> Path:
    from ..cli import versions

    doc = {
        "v": JOB_VERSION, "job_id": job.job_id, "workflow": job.workflow,
        "started": started, "finished": _now(), "versions": versions(),
        "settings": job.settings, "selection": job.selection, "output": job.output,
        "command": result.command, "command_text": command_string(result.command),
        "result": result.to_dict(), **extra,
    }
    path = Path(folder) / RUN_JSON
    path.write_text(json.dumps(jsonable(doc), indent=2, ensure_ascii=False, allow_nan=False),
                    encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# Ausfuehrung
# ---------------------------------------------------------------------------

def _failed(job: Job, exc: BaseException, issues: Optional[list] = None) -> JobResult:
    return JobResult(job.workflow, "cancelled" if isinstance(exc, _runtime.JobCancelled) else "failed",
                     exit_code_for(exc), issues=list(issues or []) + [issue_from_exception(exc)])


def execute_job(job: Job, sink: Optional[LogSink] = None,
                on_issue: Optional[Callable] = None) -> JobResult:
    """
    Job ausfuehren (Schritte siehe Modulkopf) und das ``JobResult`` liefern;
    ``output_dir`` ist der endgueltige Ergebnisordner (bei Fehlern ``None``).
    ``sink`` nimmt die Protokollzeilen auf, ``on_issue`` meldet Befunde, die vor
    dem Lauf entstehen (z.B. zu wenig Arbeitsspeicher).
    """
    sink = sink or LogSink()
    started = _now()
    try:
        wf = _workflow(job.workflow)
        settings = wf.Settings.from_dict(job.settings)
        sel = CaseSelection.from_dict(job.selection)
        out = OutputSpec.from_dict(job.output)
    except Exception as exc:  # noqa: BLE001 - kaputte Job-Datei
        return _failed(job, exc)

    hooks = getattr(wf, "job_hooks", None)
    run_settings, hook_issues, after_commit = hooks(settings) if hooks else (settings, [], None)
    final = out.target_dir(wf.default_folder(settings, sel))
    cmd = wf.command(settings, sel, final)
    pre = out.validate() + hook_issues
    if has_errors(pre + settings.validate()):
        return JobResult(job.workflow, "failed", 2, issues=pre + settings.validate(), command=cmd)
    for issue in getattr(wf, "job_checks", lambda s, sel: [])(settings, sel):
        pre.append(issue)
        if on_issue is not None:
            on_issue(issue)

    stg = job.staging_dir
    if stg.exists():
        remove_tree(stg)
    try:
        stg.mkdir(parents=True)
        sink.open(stg / RUN_LOG)
    except OSError as exc:
        sink.close()
        return _failed(job, exc, pre)

    try:
        with _runtime.use(_LoggingContext(_runtime.current(), sink)):
            result = wf.run(run_settings, sel, stg)
    finally:
        sink.close()
    result.issues = pre + result.issues
    result.command = cmd
    if not result.ok:
        remove_tree(stg)
        result.output_dir = None
        return result

    result.output_dir = str(final)
    try:
        write_run_json(stg, job, result, started, {"inputs": describe_inputs(sel)})
        committed = commit_staging(stg, final, out.policy)
    except Exception as exc:  # noqa: BLE001 - z.B. Ordner dauerhaft gesperrt
        remove_tree(stg)
        failed = _failed(job, exc, result.issues)
        failed.issues[-1] = Issue("error", "JOB.COMMIT_FAILED",
                                  "Das Ergebnis konnte nicht in den Ergebnisordner uebernommen werden.",
                                  hint_de="Ist der Ordner in einem anderen Programm geoeffnet? Erneut starten.",
                                  detail=failed.issues[-1].message_de)
        failed.command = cmd
        return failed

    changed = committed != final
    result.output_dir = str(committed)
    if after_commit is not None:
        try:
            result.outputs.update(after_commit(committed))
        except Exception as exc:  # noqa: BLE001 - Ergebnis ist gesichert, nur die Sammel-CSV fehlt
            result.issues.append(Issue("warning", "JOB.AFTER_COMMIT_FAILED",
                                       f"Nachbearbeitung fehlgeschlagen: {issue_from_exception(exc).message_de}",
                                       hint_de="Die Ergebnisdateien sind vollstaendig; nur der Nachtrag "
                                               "(z.B. die Sammel-CSV) fehlt.",
                                       detail=f"{type(exc).__name__}: {exc}"))
            if result.status == "ok":
                result.status = "ok_warnings"
        changed = True
    if changed:
        result.command = wf.command(settings, sel, committed)
        write_run_json(committed, job, result, started, {"inputs": describe_inputs(sel)})
    return result


def run_job_inprocess(job: Job, ctx: Optional[_runtime.JobContext] = None) -> JobResult:
    """
    Job im laufenden Prozess (Tests, Skripte): wie der Worker, aber ohne
    Protokoll auf stdout; ``print`` des Kerns landet in ``run.log``.
    """
    import contextlib

    sink = LogSink()
    fwd = LineForwarder(sink)
    with contextlib.redirect_stdout(fwd), _runtime.use(ctx or _runtime.current()):
        try:
            return execute_job(job, sink)
        finally:
            fwd.flush()
