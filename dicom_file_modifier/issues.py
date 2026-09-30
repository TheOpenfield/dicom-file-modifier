"""
issues.py - Strukturierte Befunde (Fehler, Warnungen, Hinweise) fuer CLI und GUI.

Ein ``Issue`` traegt einen maschinenlesbaren Code (z.B.
``CT.ORIENTATION_UNSUPPORTED``), einen deutschen Klartext, einen Hinweis
"Was tun?", optional das betroffene Eingabefeld und aufklappbare Details.
``UserInputError`` ist ein ``ValueError`` mit Issue: Die CLIs fangen
``ValueError`` bereits ab (Exit 2) und drucken den Klartext; eine Oberflaeche
kann ``err.issue`` auswerten.  Die bestehenden Konsolenausgaben bleiben
unveraendert; Issues kommen zusaetzlich in die Rueckgabewerte.

Das Modul importiert nichts aus dem Paket (Plan: Kernbausteine; die spaetere
``api/``-Schicht re-exportiert es).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

LEVELS = ("error", "warning", "info")


@dataclass(frozen=True)
class Issue:
    level: str               # "error" | "warning" | "info"
    code: str                # z.B. "CASE.FOR_KEPT"
    message_de: str          # Klartext
    hint_de: str = ""        # "Was tun?"
    field: str = ""          # betroffene Einstellung (z.B. "label")
    detail: str = ""         # aufklappbare Details

    def __post_init__(self):
        if self.level not in LEVELS:
            raise ValueError(f"Issue.level muss einer von {LEVELS} sein, nicht {self.level!r}")

    def to_dict(self) -> dict:
        return asdict(self)


class UserInputError(ValueError):
    """Eingabefehler mit ``Issue``; ``str(err)`` ist der Klartext (CLI-Exit 2)."""

    def __init__(self, issue: Issue):
        super().__init__(issue.message_de)
        self.issue = issue
