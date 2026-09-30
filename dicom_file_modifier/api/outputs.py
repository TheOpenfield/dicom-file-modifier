"""
``OutputSpec``: wohin ein Lauf schreibt (Plan: Auswahl, Ausgabe, Ergebnis).

Ergebnisse liegen in ``<root>/<folder>``; ``root`` ist absolut (der
Ergebnis-Stammordner der App), ``folder`` Default = Vorschlag des Workflows
(z.B. ``<case_id>_IDX``).  Existiert der Ordner, legt ``policy="suffix"`` einen
neuen ``<folder>_2``, ``_3`` ... an; ``overwrite`` ersetzt ihn.  Der Worker
schreibt zuerst in einen Staging-Ordner und benennt erst am Ende um.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional

from .issues import Issue

POLICIES = ("suffix", "overwrite")


def folder_name_problem(name: str) -> Optional[str]:
    """Grund, warum Windows ``name`` als Datei-/Ordnernamen ablehnt, sonst ``None``."""
    from ..case_modifier import _WIN_FORBIDDEN, _WIN_RESERVED

    if not name or not name.strip():
        return "Der Name ist leer."
    bad = sorted({c for c in name if c in _WIN_FORBIDDEN or ord(c) < 32})
    if bad:
        return "Unzulaessige Zeichen: " + " ".join(repr(c) for c in bad)
    if name.split(".")[0].upper() in _WIN_RESERVED:
        return f"{name!r} ist unter Windows ein reservierter Name."
    if name.endswith((".", " ")):
        return "Punkt oder Leerzeichen am Ende ist unter Windows nicht erlaubt."
    return None


@dataclass
class OutputSpec:
    root: str
    folder: Optional[str] = None
    policy: str = "suffix"

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "OutputSpec":
        unknown = sorted(set(d) - set(cls.__dataclass_fields__))
        if unknown:
            raise ValueError(f"Unbekannte Felder der Ausgabe: {', '.join(unknown)}")
        return cls(**d)

    def validate(self) -> list:
        issues = []
        if not self.root or not str(self.root).strip():
            issues.append(Issue("error", "OUTPUT.ROOT_MISSING", "Kein Ergebnis-Stammordner gewaehlt.",
                                field="root"))
        elif not Path(self.root).is_absolute():
            issues.append(Issue("error", "OUTPUT.ROOT_RELATIVE",
                                f"Ergebnis-Stammordner {self.root!r} ist kein absoluter Pfad.",
                                hint_de="Einen vollstaendigen Pfad waehlen, z.B. C:\\Ergebnisse.",
                                field="root"))
        if self.policy not in POLICIES:
            issues.append(Issue("error", "OUTPUT.POLICY",
                                f"Unbekannte Regel {self.policy!r} (erlaubt: {', '.join(POLICIES)}).",
                                field="policy"))
        if self.folder is not None:
            problem = folder_name_problem(self.folder)
            if problem:
                issues.append(Issue("error", "OUTPUT.FOLDER_INVALID",
                                    f"Ergebnisordner {self.folder!r}: {problem}", field="folder"))
        return issues

    def target_dir(self, default_folder: str) -> Path:
        """Endgueltiger Ergebnisordner nach ``policy`` (bei ``suffix`` der erste freie Name)."""
        name = self.folder or default_folder
        path = Path(self.root) / name
        if self.policy == "overwrite" or not path.exists():
            return path
        n = 2
        while (Path(self.root) / f"{name}_{n}").exists():
            n += 1
        return Path(self.root) / f"{name}_{n}"
