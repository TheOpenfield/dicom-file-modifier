"""
``CaseSelection``: welche Dateien ein Workflow verwendet (Plan: Auswahl, Ausgabe, Ergebnis).

Eine Auswahl nennt die Dateien explizit (CT-Schichten, RS, RD, RP,
Eclipse-Referenz) und ist damit unabhaengig vom Ordnerlayout.  Die
``from_*``-Funktionen treffen genau die Wahl der CLI im Standardlayout
(``discover_case``/``discover_dose_case``); der DICOM-Scanner der App erzeugt
Auswahlen auch fuer andere Layouts (flacher Eclipse-Export).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Optional


@dataclass
class CaseSelection:
    case_id: str = ""                    # Name in Ergebnisordner und Dateinamen
    case_dir: Optional[str] = None       # Fallordner, falls Standardlayout (fuer den CLI-Befehl)
    ct_files: list = field(default_factory=list)   # Schichten genau einer CT-Serie
    rs: Optional[str] = None
    rd: Optional[str] = None
    rp: Optional[str] = None
    eclipse_ref: Optional[str] = None    # eclipse_ref.json
    related: list = field(default_factory=list)    # weitere RP/RD des Falls (nur Hinweis)

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "CaseSelection":
        unknown = sorted(set(d) - set(cls.__dataclass_fields__))
        if unknown:
            raise ValueError(f"Unbekannte Felder der Auswahl: {', '.join(unknown)}")
        return cls(**d)

    @property
    def ct_dir(self) -> Optional[Path]:
        """Gemeinsamer Ordner der CT-Dateien (``None`` ohne CT oder bei mehreren Ordnern)."""
        parents = {Path(f).parent for f in self.ct_files}
        return parents.pop() if len(parents) == 1 else None


def _case_id(path) -> str:
    return Path(path).resolve().name


def from_rtstruct(rs: str, case_id: Optional[str] = None) -> CaseSelection:
    """Strukturanalyse: nur das RTSTRUCT (``dfm analyze RS.dcm``)."""
    p = Path(rs)
    if not p.is_file():
        raise FileNotFoundError(f"RS-Datei nicht gefunden: {rs!r}")
    return CaseSelection(case_id=case_id or _case_id(p.parent), rs=str(p))


def from_ct_dir(ct_dir: str, case_id: Optional[str] = None) -> CaseSelection:
    """Transformation nur des CT (``dfm ct-transform CT_DIR``); alle ``*.dcm`` des Ordners."""
    from .. import modifier as mod

    d = Path(ct_dir)
    if not d.is_dir():
        raise FileNotFoundError(f"CT-Ordner nicht gefunden: {ct_dir!r}")
    base = d.parent if d.name.upper() == "CT" else d
    return CaseSelection(case_id=case_id or _case_id(base), ct_files=mod.ct_dir_files(d))


def for_transform(case_dir: str, rs: Optional[str] = None) -> CaseSelection:
    """Wie ``dfm case-transform CASE [--rs RS]``: ``<case>/CT/`` + genau ein ``RS*.dcm``."""
    from .. import case_modifier as cm
    from .. import modifier as mod

    ct_dir, rs_path, siblings = cm.discover_case(case_dir, rs_override=rs, return_siblings=True)
    return CaseSelection(case_id=_case_id(case_dir), case_dir=str(case_dir),
                         ct_files=mod.ct_dir_files(ct_dir), rs=str(rs_path),
                         related=[str(p) for p in siblings])


def for_dose(case_dir: Optional[str] = None, rs: Optional[str] = None, rd: Optional[str] = None,
             rp: Optional[str] = None, eclipse_ref: Optional[str] = None) -> CaseSelection:
    """Wie ``dfm dose-indices [CASE] [--rs/--rd/--rp/--eclipse-ref]`` (PLAN-RD bei mehreren)."""
    from .. import dose_indices as di
    from .. import modifier as mod

    files = di.discover_dose_case(case_dir, rs, rd, rp, eclipse_ref)
    return CaseSelection(
        case_id=files["case_id"], case_dir=str(case_dir) if case_dir is not None else None,
        ct_files=mod.ct_dir_files(files["ct_dir"]) if files["ct_dir"] else [],
        rs=str(files["rs"]), rd=str(files["rd"]),
        rp=str(files["rp"]) if files["rp"] else None,
        eclipse_ref=str(files["eclipse_ref"]) if files["eclipse_ref"] else None,
    )
