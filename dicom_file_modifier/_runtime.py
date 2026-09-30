"""
_runtime.py - Laufzeitkontext fuer lange Rechnungen: Stufen, Fortschritt, Abbruch.

Der Kern fragt den Kontext in seinen Schleifen mit ``current()`` ab, statt ihn
als Parameter durchzureichen; so bleiben die Signaturen der Kernfunktionen
unveraendert (Plan, Kernbausteine).  Ohne gesetzten Kontext liefert
``current()`` einen ``JobContext``, der nichts tut: CLI-Laeufe verhalten sich
wie bisher.  Der Worker der Desktop-App (Plan P0.8) setzt mit ``use()`` einen
eigenen Kontext, der Fortschritt meldet und ``check_cancel()`` scharf schaltet.

Das Modul importiert nichts aus dem Paket.
"""

from __future__ import annotations

import contextlib
import contextvars
from typing import Iterator


class JobCancelled(Exception):
    """Abbruch auf Wunsch des Nutzers, ausgeloest an einem Checkpoint."""


class JobContext:
    """Basisklasse und Null-Kontext: alle Methoden sind No-ops."""

    def stage(self, key: str, label: str) -> None:
        """Beginn einer Stufe (``key`` maschinenlesbar, ``label`` deutscher Klartext)."""

    def progress(self, done: float, total: float, text: str = "") -> None:
        """Fortschritt innerhalb der aktuellen Stufe."""

    def log(self, message: str, level: str = "info") -> None:
        """Protokollmeldung (neuer Code; bestehende ``print``-Ausgaben bleiben)."""

    def check_cancel(self) -> None:
        """Wirft ``JobCancelled``, wenn ein Abbruch angefordert wurde."""


NullContext = JobContext

_CURRENT: contextvars.ContextVar = contextvars.ContextVar("dfm_job_context", default=JobContext())


def current() -> JobContext:
    """Der fuer den laufenden Thread gesetzte Kontext (Default: No-op)."""
    return _CURRENT.get()


@contextlib.contextmanager
def use(ctx: JobContext) -> Iterator[JobContext]:
    """``with use(ctx): ...`` setzt ``ctx`` fuer die Dauer des Blocks."""
    token = _CURRENT.set(ctx)
    try:
        yield ctx
    finally:
        _CURRENT.reset(token)
