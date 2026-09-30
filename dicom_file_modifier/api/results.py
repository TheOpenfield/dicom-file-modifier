"""
``JobResult`` und die gemeinsame Ausfuehrungshuelle der ``run``-Funktionen.

``run_guarded`` fuehrt den Rumpf eines Workflows unter einem ``StageTimer``
aus (misst die Stufen und leitet Fortschritt/Abbruch an den umgebenden
``JobContext`` weiter), faengt jede Ausnahme als ``Issue`` ab und baut das
Ergebnis: Status, Exit-Code, Ausgaben relativ zum Ergebnisordner, Manifest.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional

from .. import _runtime
from .issues import Issue, exit_code_for, issue_from_exception

STATUSES = ("ok", "ok_warnings", "failed", "cancelled")


@dataclass
class JobResult:
    workflow: str
    status: str                      # ok | ok_warnings | failed | cancelled
    exit_code: int                   # 0 ok, 2 Eingabefehler, 3 abgebrochen, 1 interner Fehler
    output_dir: Optional[str] = None
    outputs: dict = field(default_factory=dict)    # Rolle -> Pfad (relativ zu output_dir, sonst absolut)
    issues: list = field(default_factory=list)     # Issue
    summary: dict = field(default_factory=dict)    # kleine Kennzahlen fuer die Anzeige
    timings: dict = field(default_factory=dict)    # {"total_s", "stages": [{key, label, seconds}]}
    manifest: list = field(default_factory=list)   # [{"path", "bytes"}] aller Dateien in output_dir
    command: list = field(default_factory=list)    # aequivalente dfm-Befehle (argv-Listen)

    @property
    def ok(self) -> bool:
        return self.status in ("ok", "ok_warnings")

    def to_dict(self) -> dict:
        d = dict(self.__dict__)
        d["issues"] = [i.to_dict() for i in self.issues]
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "JobResult":
        d = dict(d)
        d["issues"] = [Issue(**i) for i in d.get("issues", [])]
        return cls(**d)


class StageTimer(_runtime.JobContext):
    """Leitet an ``parent`` weiter und misst die Dauer jeder Stufe."""

    def __init__(self, parent: _runtime.JobContext):
        self.parent = parent
        self.stages: list = []
        self._t0 = time.perf_counter()
        self._current = None

    def stage(self, key: str, label: str) -> None:
        self._close()
        self._current = (key, label, time.perf_counter())
        self.parent.stage(key, label)

    def progress(self, done: float, total: float, text: str = "") -> None:
        self.parent.progress(done, total, text)

    def log(self, message: str, level: str = "info") -> None:
        self.parent.log(message, level)

    def check_cancel(self) -> None:
        self.parent.check_cancel()

    def _close(self) -> None:
        if self._current is not None:
            key, label, t = self._current
            self.stages.append({"key": key, "label": label,
                                "seconds": round(time.perf_counter() - t, 3)})
            self._current = None

    def timings(self) -> dict:
        self._close()
        return {"total_s": round(time.perf_counter() - self._t0, 3), "stages": list(self.stages)}


def build_manifest(out_dir) -> list:
    """Alle Dateien unter ``out_dir`` mit Groesse, Pfade relativ (``/``), sortiert."""
    root = Path(out_dir)
    if not root.is_dir():
        return []
    return [{"path": p.relative_to(root).as_posix(), "bytes": p.stat().st_size}
            for p in sorted(root.rglob("*")) if p.is_file()]


def relative_outputs(outputs: dict, out_dir) -> dict:
    """Pfade unter ``out_dir`` relativ (``/``), andere absolut; ``None``-Eintraege entfallen."""
    root = Path(out_dir).resolve()
    out = {}
    for role, path in outputs.items():
        if not path:
            continue
        p = Path(path).resolve()
        try:
            out[role] = p.relative_to(root).as_posix()
        except ValueError:
            out[role] = str(p)
    return out


def status_for(issues: list) -> str:
    return "ok_warnings" if any(i.level in ("warning", "error") for i in issues) else "ok"


def run_guarded(workflow: str, out_dir, body: Callable[[], tuple],
                pre_issues: Optional[list] = None, command: Optional[list] = None) -> JobResult:
    """
    ``body()`` liefert ``(outputs, issues, summary)``.  Ausnahmen werden zu
    ``failed`` (Exit 2 bei Eingabefehlern, sonst 1) bzw. ``cancelled`` (Exit 3);
    ``pre_issues`` (z.B. Warnungen der Einstellungspruefung) stehen vorn.
    """
    timer = StageTimer(_runtime.current())
    pre = list(pre_issues or [])
    try:
        with _runtime.use(timer):
            outputs, issues, summary = body()
    except Exception as exc:  # noqa: BLE001 - jeder Fehler wird ein Issue im Ergebnis
        status = "cancelled" if isinstance(exc, _runtime.JobCancelled) else "failed"
        return JobResult(workflow, status, exit_code_for(exc), str(out_dir),
                         issues=pre + [issue_from_exception(exc)], timings=timer.timings(),
                         manifest=build_manifest(out_dir), command=list(command or []))
    issues = pre + list(issues)
    return JobResult(workflow, status_for(issues), 0, str(out_dir),
                     outputs=relative_outputs(outputs, out_dir), issues=issues,
                     summary=summary, timings=timer.timings(), manifest=build_manifest(out_dir),
                     command=list(command or []))


def blocked(workflow: str, out_dir, issues: list, command: Optional[list] = None) -> JobResult:
    """Ergebnis ohne Lauf: die Pruefung vor dem Start hat Fehler gefunden (Exit 2)."""
    return JobResult(workflow, "failed", 2, str(out_dir), issues=list(issues),
                     command=list(command or []))
