"""
dose_indices.py - Automatische Indexberechnung aus RTSTRUCT + RTDOSE (+ RTPLAN)

Berechnet je Zielvolumen:
  - Paddick Conformity Index   CI = (TV&PIV)^2 / (TV * PIV)
  - Coverage, Selektivitaet, RTOG-CI (PIV/TV), Dice
  - Gradient Index             GI = PIV50 / PIV100
  - Gradient Measure           GM = r_eq(PIV50) - r_eq(PIV100)  [cm]
  - Homogenitaetsindex ICRU83  HI = (D2 - D98) / D50
  - Dmin/Dmax/Dmean, D95, V95/V100 des Ziels
und schreibt eine separate RTSTRUCT mit Isodosen-ROIs (100 %, 50 %, ...) und
den fuer den Index benutzten Hilfskonturen (Schnitt, Unterdosierung, Spill).

Abgleich Eclipse: die ``DVHSequence`` der RTDOSE (Eclipse-DVHs) wird immer
gelesen (TV, V_Rx, D98/D50/D2 der Ziele; PIV/PIV50 nur bei einem Body-DVH);
weitere Eclipse-Komponentenwerte kommen aus ``<case>/eclipse_ref.json``
(``--eclipse-ref``) oder ``--eclipse-values TV=..,VRX=..,PIV=..,V50=..``.
Abweichungen ueber ``--eclipse-tol-pct`` (Default 5 %) werden markiert.

Validierungsansicht (``dose_viz``): validation.html (3D, Schichtbrowser, DVH,
Indextabelle; offline) und dose_overview.png; ``--no-viz`` / ``--no-viz-ct``.

Verwendung:
  python -m dicom_file_modifier.dose_indices data/<case-id> [Optionen]
  python -m dicom_file_modifier.dose_indices --rs RS.dcm --rd RD.dcm [--rp RP.dcm] [Optionen]
  python -m dicom_file_modifier.dose_indices --self-test

Ausgaben in output/<case-id><label>/:
  RS_<case-id><label>.dcm, <case-id>_indices.json, indices.txt, indices.csv,
  validation.html, dose_overview.png

Alle Volumina in cm3, Dosen in Gy, LPS-Koordinaten in mm.  Konsolenausgabe
ist ASCII (Windows-Konsole), Dateien sind UTF-8.
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import json
import math
import os
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import pydicom
from pydicom.valuerep import format_number_as_ds
from scipy.integrate import trapezoid

from . import _runtime
from . import analyzer as ana
from . import dose as dm
from . import rtstruct_writer as rw
from ._compat import PYDICOM_MAJOR
# Seit P0.3 in dose_constants (loest den Zyklus mit rtstruct_writer/dose_viz);
# hier weiter importierbar (alte Importpfade)
from .dose_constants import HELPER_COLORS, LEVEL_COLORS, TOOL_NAME, TOOL_VERSION  # noqa: F401

DEFAULT_ISODOSE = "100,50"
GRID_CHOICES = (1.0, 0.5, 0.25, 0.1)
INTERP_ORDER = {"linear": 1, "cubic": 3}

CSV_COLUMNS = [
    "case_id", "target", "rx_gy", "rx_source", "grid_mm", "dose_interp", "volume_model",
    "piv_scope", "tv_cm3", "piv_cm3", "tv_piv_cm3", "piv50_cm3", "ci_paddick", "ci_rtog",
    "coverage", "dice", "gi", "gm_cm", "hi_icru83", "d2_gy", "d50_gy", "d98_gy", "dmin_gy",
    "dmax_gy", "dmean_gy", "d95_gy", "v95_pct", "v100_pct", "rs_file", "rd_sop_uid",
    "run_timestamp",
    # Abgleich Eclipse (leer ohne Referenz)
    "ecl_source", "ecl_tv_cm3", "ecl_tv_piv_cm3", "ecl_piv_cm3", "ecl_piv50_cm3",
    "ecl_ci_paddick", "ecl_gi", "ecl_hi_icru83", "ecl_d98_gy", "ecl_d50_gy", "ecl_d2_gy",
    "d_tv_pct", "d_tv_piv_pct", "d_piv_pct", "d_ci_paddick_pct", "d_gi_pct",
    "d_hi_icru83_pct", "d_d98_pct", "ecl_tol_pct", "ecl_n_flagged",
]


# ---------------------------------------------------------------------------
# 1. Discovery, Zielauswahl, Verschreibung, Level
# ---------------------------------------------------------------------------

def discover_dose_case(case_dir: Optional[str], rs_override: Optional[str] = None,
                       rd_override: Optional[str] = None,
                       rp_override: Optional[str] = None,
                       eclipse_ref_override: Optional[str] = None,
                       need_rd: bool = True) -> dict:
    """
    ``{'case_id', 'case_dir', 'rs', 'rd', 'rp'|None, 'ct_dir'|None, 'eclipse_ref'|None}``.
    Sucht ``RS*.dcm``/``RD*.dcm``/``RP*.dcm`` im Case-Ordner; Overrides haben
    Vorrang.  Bei mehreren RD-Kandidaten wird die PLAN-Summendosis bevorzugt,
    die den RP referenziert; bleibt es mehrdeutig -> ``ValueError``.
    ``eclipse_ref`` = ``--eclipse-ref`` (muss existieren) oder
    ``<case>/eclipse_ref.json``, falls vorhanden.  ``need_rd=False`` (``--list``)
    sucht kein RD (``rd`` ist dann None bzw. der Override).
    """
    def _pick(kind: str, override: Optional[str], required: bool):
        if override is not None:
            p = Path(override)
            if not p.is_file():
                raise FileNotFoundError(f"{kind}-Datei nicht gefunden: {override!r}")
            return p
        if case_dir is None:
            if required:
                raise ValueError(f"Ohne Case-Ordner muss --{kind.lower()} angegeben werden.")
            return None
        cands = sorted(Path(case_dir).glob(f"{kind}*.dcm"))
        if not cands:
            if required:
                raise FileNotFoundError(
                    f"Keine '{kind}*.dcm'-Datei in {case_dir!r} gefunden "
                    f"(anderes Praefix? dann --{kind.lower()} <pfad>)."
                )
            return None
        if len(cands) == 1:
            return cands[0]
        if kind == "RD":
            plans = []
            for c in cands:
                try:
                    h = pydicom.dcmread(str(c), stop_before_pixels=True)
                except Exception:
                    continue
                if str(h.get("DoseSummationType", "")).upper() == "PLAN":
                    plans.append(c)
            if len(plans) == 1:
                return plans[0]
        joined = "\n  ".join(str(p) for p in cands)
        raise ValueError(
            f"Mehrere '{kind}*.dcm'-Kandidaten in {case_dir!r}:\n  {joined}\n"
            f"Bitte mit --{kind.lower()} <pfad> explizit auswaehlen."
        )

    if case_dir is not None and not Path(case_dir).is_dir():
        raise FileNotFoundError(f"Case-Ordner nicht gefunden: {case_dir!r}")
    rs = _pick("RS", rs_override, True)
    rd = _pick("RD", rd_override, True) if (need_rd or rd_override is not None) else None
    rp = _pick("RP", rp_override, False)
    base = Path(case_dir) if case_dir is not None else (rd or rs).parent
    ct_dir = base / "CT"
    if eclipse_ref_override is not None:
        eclipse_ref = Path(eclipse_ref_override)
        if not eclipse_ref.is_file():
            raise FileNotFoundError(f"Eclipse-Referenzdatei nicht gefunden: {eclipse_ref_override!r}")
    else:
        cand = base / "eclipse_ref.json"
        eclipse_ref = cand if cand.is_file() else None
    return {
        "case_id": base.resolve().name,
        "case_dir": base,
        "rs": rs, "rd": rd, "rp": rp,
        "ct_dir": ct_dir if ct_dir.is_dir() else None,
        "eclipse_ref": eclipse_ref,
    }


def dose_files(rs, rd, rp=None, ct=None, eclipse_ref=None,
               case_id: Optional[str] = None) -> dict:
    """
    Wie ``discover_dose_case``, aber fuer bereits gewaehlte Dateien (API):
    ``ct`` ist ein CT-Ordner oder eine Liste von CT-Dateien (``None`` = ohne CT),
    ``eclipse_ref`` wird nicht automatisch gesucht.  Fehlende Dateien ->
    ``FileNotFoundError``; ``case_id`` Default = Ordnername des RD.
    """
    def _file(kind: str, path):
        if path is None:
            return None
        p = Path(path)
        if not p.is_file():
            raise FileNotFoundError(f"{kind}-Datei nicht gefunden: {str(path)!r}")
        return p

    rs_p, rd_p = _file("RS", rs), _file("RD", rd)
    if rs_p is None or rd_p is None:
        raise ValueError("RTSTRUCT und RTDOSE muessen angegeben werden.")
    ct_dir = ct_files = None
    if isinstance(ct, (str, os.PathLike)):
        ct_dir = Path(ct)
        if not ct_dir.is_dir():
            raise FileNotFoundError(f"CT-Ordner nicht gefunden: {str(ct)!r}")
    elif ct:
        ct_files = [str(f) for f in ct]
        missing = [f for f in ct_files if not Path(f).is_file()]
        if missing:
            raise FileNotFoundError(f"{len(missing)} CT-Datei(en) nicht gefunden, z.B. {missing[0]!r}")
        parents = {Path(f).parent for f in ct_files}
        ct_dir = parents.pop() if len(parents) == 1 else None
    return {
        "case_id": case_id or rd_p.parent.resolve().name,
        "case_dir": rd_p.parent,
        "rs": rs_p, "rd": rd_p, "rp": _file("RP", rp),
        "ct_dir": ct_dir, "ct_files": ct_files,
        "eclipse_ref": _file("Eclipse-Referenz", eclipse_ref),
    }


def roi_table(rs_ds: pydicom.Dataset) -> list:
    """``[(roi_number, name, rt_type, category)]`` fuer alle ROIs (Kategorie aus ``analyzer``)."""
    names = ana.get_structure_names(rs_ds)
    types = ana.get_structure_type(rs_ds)
    geoms = ana.get_structure_geom_types(rs_ds)
    out = []
    for num, name in names.items():
        cat = ana.classify_structure(name, types.get(num, ""), geoms.get(num, set()))
        out.append((int(num), name, types.get(num, ""), cat))
    return out


_roi_table = roi_table        # alter Name (bis P0.4 privat)


def target_candidates(rs_ds: pydicom.Dataset, rp_refs: Optional[list] = None) -> dict:
    """
    Zielvorschlag ohne Abbruch, fuer die Auswahl in einer Oberflaeche:
    ``{'rois': roi_table, 'targets': TARGET-ROIs, 'ptvs': PTVs darunter,
    'default': [(roi_number, name), ...], 'reason': Begruendung, 'notes': [...]}``.
    ``default`` ist die Auto-Auswahl von ``select_targets``: alle PTVs; passt
    eine RTPLAN-``DoseReferenceDescription`` (SH, 16 Zeichen) als Praefix auf
    genau eines von mehreren, nur dieses.  Ohne PTV ist ``default`` leer.
    """
    table = roi_table(rs_ds)
    targets = [r for r in table if r[3] == ana.CAT_TARGET]
    ptvs = [r for r in targets if r[1].upper().startswith("PTV")]
    chosen, notes = ptvs, []
    reason = "alle PTVs (Kategorie TARGET, Name beginnt mit PTV)" if ptvs else "kein PTV gefunden"
    for ref in rp_refs or []:
        desc = (ref.get("description") or "").strip()
        if not desc:
            continue
        hits = [r for r in ptvs if r[1].lower().startswith(desc.lower())]
        if len(hits) == 1 and len(ptvs) > 1:
            chosen = hits
            rest = ", ".join(repr(r[1]) for r in ptvs if r not in hits)
            notes.append(
                f"RTPLAN-Verschreibung '{desc}' passt auf {hits[0][1]!r}; "
                f"weitere PTVs ({rest}) nicht ausgewertet (--target fuer alle)."
            )
            reason = f"RTPLAN-Verschreibung '{desc}'"
            break
    return {"rois": table, "targets": targets, "ptvs": ptvs,
            "default": [(r[0], r[1]) for r in chosen], "reason": reason, "notes": notes}


def select_targets(rs_ds: pydicom.Dataset, target_arg: Optional[str],
                   rp_refs: Optional[list] = None) -> tuple:
    """
    Liefert ``([(roi_number, name), ...], hinweise)``.
    ``--target``: exakter Name, sonst eindeutiger case-insensitiver Teilstring.
    Auto: ``target_candidates()['default']``; ohne PTV -> ``ValueError``.
    """
    if target_arg:
        table = roi_table(rs_ds)
        chosen = []
        for token in [t.strip() for t in target_arg.split(",") if t.strip()]:
            exact = [r for r in table if r[1] == token]
            if not exact:
                exact = [r for r in table if r[1].lower() == token.lower()]
            if not exact:
                exact = [r for r in table if token.lower() in r[1].lower()]
            if len(exact) != 1:
                cands = ", ".join(repr(r[1]) for r in exact) if exact else "keine"
                raise ValueError(
                    f"--target {token!r} ist nicht eindeutig (Kandidaten: {cands}). "
                    "Verfuegbare ROIs mit --list anzeigen."
                )
            if exact[0] not in chosen:
                chosen.append(exact[0])
        return [(r[0], r[1]) for r in chosen], []

    cand = target_candidates(rs_ds, rp_refs)
    if not cand["ptvs"]:
        listing = ", ".join(f"{r[1]!r}" for r in cand["targets"]) or "keine"
        raise ValueError(
            "Kein PTV gefunden. Zielvolumen mit --target NAME waehlen "
            f"(TARGET-Kandidaten: {listing})."
        )
    return cand["default"], cand["notes"]


def rx_candidates(rp_refs: Optional[list], dose: Optional[dm.DoseGrid] = None,
                  target_names: Optional[list] = None) -> dict:
    """
    Verschreibungsvorschlaege ohne Abbruch, fuer die Auswahl in einer
    Oberflaeche: ``{'refs': RTPLAN-DoseReferences mit TargetPrescriptionDose
    (Typ TARGET oder leer), 'default': ref | None, 'reason': Begruendung,
    'dmax_gy': Dmax fuer "% von Dmax" | None}``.  ``default`` ist die Wahl von
    ``resolve_prescription`` ohne ``--rx``: die einzige Verschreibung oder die,
    deren Beschreibung Praefix genau eines Zielnamens ist; sonst ``None``.
    """
    refs = [r for r in (rp_refs or [])
            if r.get("target_prescription_dose_gy") is not None
            and (r.get("reference_type", "").upper() in ("TARGET", ""))]
    default, reason = None, "keine Verschreibung im RTPLAN"
    if len(refs) == 1:
        default, reason = refs[0], "einzige Verschreibung im RTPLAN"
    elif len(refs) > 1:
        lowered = [n.lower() for n in (target_names or [])]
        pref = [r for r in refs if r.get("description")
                and any(n.startswith(r["description"].lower()) for n in lowered)]
        if len(pref) == 1:
            default, reason = pref[0], f"Beschreibung '{pref[0]['description']}' passt zum Ziel"
        else:
            reason = "mehrere Verschreibungen, keine passt eindeutig zu den Zielen"
    return {"refs": refs, "default": default, "reason": reason,
            "dmax_gy": float(dose.dmax) if dose is not None else None}


def resolve_prescription(rx_cli: Optional[float], rx_pct_of_max: Optional[float],
                         rp_refs: Optional[list], dose: dm.DoseGrid,
                         target_names: list) -> tuple:
    """(rx_gy, source in {'cli','pct_of_max','rtplan'}, detail_text)."""
    if rx_cli is not None and rx_pct_of_max is not None:
        raise ValueError("--rx und --rx-pct-of-max schliessen sich aus.")
    if rx_cli is not None:
        if rx_cli <= 0:
            raise ValueError("--rx muss > 0 sein.")
        return float(rx_cli), "cli", f"{rx_cli:.2f} Gy per --rx"
    if rx_pct_of_max is not None:
        if not (0 < rx_pct_of_max <= 100):
            raise ValueError("--rx-pct-of-max muss in (0, 100] liegen.")
        rx = rx_pct_of_max / 100.0 * dose.dmax
        return float(rx), "pct_of_max", f"{rx_pct_of_max:g} % von Dmax {dose.dmax:.2f} Gy"
    cand = rx_candidates(rp_refs, dose, target_names)
    refs = cand["refs"]
    if not refs:
        raise ValueError(
            "Keine Verschreibung gefunden (RTPLAN fehlt oder ohne TargetPrescriptionDose). "
            "Bitte --rx <Gy> oder --rx-pct-of-max <Prozent> angeben."
        )
    r = cand["default"]
    if r is None:
        listing = "; ".join(f"{x['description']!r}: {x['target_prescription_dose_gy']:g} Gy"
                            for x in refs)
        raise ValueError(
            f"Mehrere Verschreibungen im RTPLAN ({listing}). Bitte --rx angeben."
        )
    return (float(r["target_prescription_dose_gy"]), "rtplan",
            f"DoseReferenceSequence '{r['description']}'")


def parse_isodose_levels(spec: str, rx_gy: float) -> tuple:
    """
    ``'100,50,80,12Gy'`` -> ``[{'key','label','pct','gy'}, ...]`` absteigend
    nach Gy; 100 und 50 werden bei Bedarf ergaenzt (Hinweis in ``notes``).
    Ein Level, das auf Rx bzw. Rx/2 faellt (auch in Gy angegeben, z.B. ``20Gy``
    bei Rx 20 Gy), bekommt den Schluessel ``'100'`` bzw. ``'50'``; die Indizes
    greifen ueber diese Schluessel zu.
    """
    levels, notes = [], []
    for tok in [t.strip() for t in (spec or "").split(",") if t.strip()]:
        low = tok.lower()
        try:
            if low.endswith("gy"):
                gy = float(low[:-2])
                pct = gy / rx_gy * 100.0
                key, label = f"{gy:g}Gy", f"{gy:g} Gy"
            else:
                pct = float(low.rstrip("%"))
                gy = pct / 100.0 * rx_gy
                key, label = f"{pct:g}", f"{pct:g}%"
        except ValueError:
            raise ValueError(f"Ungueltiges Isodosen-Level: {tok!r} (z.B. 100,50,80 oder 12Gy)")
        if gy <= 0:
            raise ValueError(f"Isodosen-Level muss > 0 sein: {tok!r}")
        for need, canon in ((100.0, "100"), (50.0, "50")):
            if abs(pct - need) < 1e-6:
                key = canon
        levels.append({"key": key, "label": label, "pct": pct, "gy": gy})
    for need, lab in ((100.0, "100"), (50.0, "50")):
        if not any(abs(lv["pct"] - need) < 1e-6 for lv in levels):
            levels.append({"key": lab, "label": f"{lab}%", "pct": need, "gy": need / 100.0 * rx_gy})
            notes.append(f"Isodosen-Level {lab} % ergaenzt (fuer CI/GI erforderlich).")
    # Duplikate (gleiches Gy) entfernen, absteigend sortieren
    uniq = {}
    for lv in levels:
        uniq.setdefault(round(lv["gy"], 6), lv)
    out = sorted(uniq.values(), key=lambda lv: -lv["gy"])
    return out, notes


# ---------------------------------------------------------------------------
# 2. Reine Indexfunktionen (cm3; NaN bei Nenner 0)
# ---------------------------------------------------------------------------

def _safe_div(a: float, b: float) -> float:
    return float(a) / float(b) if b else float("nan")


def paddick_ci(tv: float, piv: float, tv_piv: float) -> float:
    return _safe_div(tv_piv * tv_piv, tv * piv)


def coverage(tv: float, tv_piv: float) -> float:
    return _safe_div(tv_piv, tv)


def selectivity(piv: float, tv_piv: float) -> float:
    return _safe_div(tv_piv, piv)


def rtog_ci(tv: float, piv: float) -> float:
    return _safe_div(piv, tv)


def dice_index(tv: float, piv: float, tv_piv: float) -> float:
    return _safe_div(2.0 * tv_piv, tv + piv)


def gradient_index(piv50: float, piv100: float) -> float:
    return _safe_div(piv50, piv100)


def equivalent_radius_cm(v_cm3: float) -> float:
    return (3.0 * v_cm3 / (4.0 * math.pi)) ** (1.0 / 3.0) if v_cm3 > 0 else float("nan")


def gradient_measure_cm(piv50: float, piv100: float) -> float:
    return equivalent_radius_cm(piv50) - equivalent_radius_cm(piv100)


def homogeneity_index(d2: float, d98: float, d50: float) -> float:
    return _safe_div(d2 - d98, d50)


# ---------------------------------------------------------------------------
# 2b. Eclipse-Referenzwerte (DVHSequence, eclipse_ref.json, --eclipse-values)
# ---------------------------------------------------------------------------

ECLIPSE_KEYS = ("tv_cm3", "tv_piv_cm3", "piv_cm3", "piv50_cm3", "ci_paddick", "gi", "hi_icru83",
                "d98_gy", "d50_gy", "d2_gy", "dmean_gy", "dmin_gy", "dmax_gy")
ECLIPSE_ALIASES = {
    "tv_cm3": ("tv", "volume", "volume_cm3"),
    "tv_piv_cm3": ("vrx", "v_rx", "v_rx_cm3", "tv_piv", "tvpiv"),
    "piv_cm3": ("piv", "piv100"),
    "piv50_cm3": ("v50", "piv50", "v_half_rx", "v_half_rx_cm3"),
    "ci_paddick": ("ci", "paddick"),
    "gi": ("gradient_index",),
    "hi_icru83": ("hi", "hi_icru"),
    "d98_gy": ("d98",), "d50_gy": ("d50",), "d2_gy": ("d2",),
    "dmean_gy": ("dmean", "mean"), "dmin_gy": ("dmin", "min"), "dmax_gy": ("dmax", "max"),
}
ECLIPSE_LABELS = {
    "tv_cm3": "TV [cm3]", "tv_piv_cm3": "TV&PIV [cm3]", "piv_cm3": "PIV [cm3]",
    "piv50_cm3": "PIV50 [cm3]", "ci_paddick": "CI Paddick", "gi": "GI", "hi_icru83": "HI ICRU83",
    "d98_gy": "D98 [Gy]", "d50_gy": "D50 [Gy]", "d2_gy": "D2 [Gy]", "dmean_gy": "Dmean [Gy]",
    "dmin_gy": "Dmin [Gy]", "dmax_gy": "Dmax [Gy]",
}
_VGY_RE = re.compile(r"^v([0-9]+(?:[.,][0-9]+)?)gy$")


def normalize_eclipse_key(key: str, rx_gy: Optional[float] = None) -> str:
    """
    Alias -> kanonischer Schluessel (case-insensitiv).  ``V<n>Gy`` wird ueber
    den Wert aufgeloest: n == Rx -> ``tv_piv_cm3``, n == Rx/2 -> ``piv50_cm3``,
    sonst ``ValueError``.  Unbekannte Schluessel -> ``ValueError`` mit Liste.
    """
    k = str(key).strip().lower().replace(" ", "").replace("-", "_")
    if k in ECLIPSE_KEYS:
        return k
    for canon, aliases in ECLIPSE_ALIASES.items():
        if k in aliases:
            return canon
    m = _VGY_RE.match(k)
    if m:
        gy = float(m.group(1).replace(",", "."))
        if rx_gy is None:
            raise ValueError(f"Eclipse-Wert {key!r}: V<n>Gy braucht die Verschreibung (Rx unbekannt).")
        if abs(gy - rx_gy) < 1e-3:
            return "tv_piv_cm3"
        if abs(gy - 0.5 * rx_gy) < 1e-3:
            return "piv50_cm3"
        raise ValueError(f"Eclipse-Wert {key!r} passt weder zu Rx {rx_gy:.2f} Gy noch zu "
                         f"Rx/2 {0.5 * rx_gy:.2f} Gy.")
    known = ", ".join(f"{c} ({'/'.join(a)})" for c, a in ECLIPSE_ALIASES.items())
    raise ValueError(f"Unbekannter Eclipse-Schluessel {key!r}. Bekannt: {known}, V<Rx>Gy, V<Rx/2>Gy.")


@dataclass
class EclipseReference:
    """Eclipse-Referenzwerte je Ziel: ``values[ziel][key] = {'value', 'source'}``."""
    values: dict = field(default_factory=dict)
    meta: dict = field(default_factory=dict)
    notes: list = field(default_factory=list)

    def get(self, target: str, key: str) -> Optional[dict]:
        return self.values.get(target, {}).get(key)

    def value(self, target: str, key: str) -> Optional[float]:
        e = self.get(target, key)
        return e["value"] if e else None

    def set(self, target: str, key: str, value: float, source: str) -> None:
        self.values.setdefault(target, {})[key] = {"value": float(value), "source": source}

    def merge(self, other: "EclipseReference") -> None:
        """``other`` gewinnt bei gleichem Ziel/Schluessel."""
        for t, vals in other.values.items():
            for k, v in vals.items():
                self.values.setdefault(t, {})[k] = dict(v)
        self.meta.update(other.meta)
        self.notes.extend(other.notes)

    def is_empty(self) -> bool:
        return not any(self.values.values())

    def sources(self) -> list:
        return sorted({v["source"] for vals in self.values.values() for v in vals.values()})


def match_target_name(token: str, target_names: list, allow_unknown: bool = False) -> Optional[str]:
    """
    exakt -> case-insensitiv -> eindeutiger Praefix; ``*``/leer nur bei genau
    einem Ziel.  ``allow_unknown``: passt der Name auf gar kein Ziel, ``None``
    statt ``ValueError`` (mehrdeutige Namen bleiben ein Fehler).
    """
    tok = (token or "").strip()
    if tok in ("", "*"):
        if len(target_names) == 1:
            return target_names[0]
        raise ValueError(
            "Eclipse-Referenz: ohne Zielname ('*') nur bei genau einem Ziel erlaubt; bei mehreren "
            "Zielen NAME:KEY=WERT (--eclipse-values) bzw. den Zielnamen als JSON-Schluessel verwenden."
        )
    if tok in target_names:
        return tok
    ci = [n for n in target_names if n.lower() == tok.lower()]
    if len(ci) == 1:
        return ci[0]
    pre = [n for n in target_names if n.lower().startswith(tok.lower())]
    if len(pre) == 1:
        return pre[0]
    if not pre and allow_unknown:
        return None
    cands = ", ".join(repr(n) for n in (pre or target_names)) or "keine"
    raise ValueError(f"Eclipse-Referenz: Ziel {token!r} nicht eindeutig (Kandidaten: {cands}).")


def _to_float(v, where: str) -> Optional[float]:
    if v is None:
        return None
    if isinstance(v, str):
        v = v.strip().replace(",", ".")
    try:
        f = float(v)
    except (TypeError, ValueError):
        raise ValueError(f"Eclipse-Referenz: Wert {v!r} fuer {where} ist nicht numerisch.")
    return None if math.isnan(f) else f


def eclipse_reference_from_dict(d: dict, target_names: list, rx_gy: Optional[float],
                                source: str) -> EclipseReference:
    """
    ``{"<Ziel>": {alias: wert, ...}, "_meta": {...}}`` -> ``EclipseReference``.
    ``"*"`` als Zielschluessel bei genau einem Ziel; ``null`` wird uebersprungen.
    Ziele, die in diesem Lauf nicht ausgewertet werden (die Datei darf mehr
    Ziele enthalten), ergeben einen Hinweis statt eines Fehlers.
    """
    ref = EclipseReference()
    if not isinstance(d, dict):
        raise ValueError("Eclipse-Referenz: erwartet ein JSON-Objekt {Ziel: {Schluessel: Wert}}.")
    for tkey, vals in d.items():
        if tkey == "_meta":
            if isinstance(vals, dict):
                ref.meta.update(vals)
            continue
        if not isinstance(vals, dict):
            raise ValueError(f"Eclipse-Referenz: Eintrag {tkey!r} muss ein Objekt mit Schluessel/Wert-Paaren sein.")
        target = match_target_name(tkey, target_names, allow_unknown=True)
        if target is None:
            evaluated = ", ".join(repr(n) for n in target_names) or "keine"
            ref.notes.append(f"Eclipse-Referenz: Ziel {tkey!r} wird in diesem Lauf nicht ausgewertet "
                             f"(ausgewertet: {evaluated}); Werte ignoriert.")
            continue
        for k, v in vals.items():
            fv = _to_float(v, f"{tkey}/{k}")
            if fv is None:
                continue
            ref.set(target, normalize_eclipse_key(k, rx_gy), fv, source)
    return ref


def load_eclipse_reference_json(path, target_names: list, rx_gy: Optional[float]) -> EclipseReference:
    """``eclipse_ref.json`` (UTF-8) lesen; Quelle ``json``."""
    p = Path(path)
    with open(p, encoding="utf-8") as fh:
        try:
            d = json.load(fh)
        except json.JSONDecodeError as e:
            raise ValueError(f"Eclipse-Referenz {p}: ungueltiges JSON ({e}).")
    ref = eclipse_reference_from_dict(d, target_names, rx_gy, "json")
    ref.meta["json_path"] = str(p)
    return ref


def parse_eclipse_values(spec: str, target_names: list, rx_gy: Optional[float]) -> EclipseReference:
    """``--eclipse-values "[ZIEL:]KEY=WERT,..."`` (Dezimalpunkt) -> Quelle ``cli``."""
    ref = EclipseReference()
    for tok in [t.strip() for t in (spec or "").split(",") if t.strip()]:
        if "=" not in tok:
            raise ValueError(f"--eclipse-values: {tok!r} hat kein '=' (Form [ZIEL:]SCHLUESSEL=WERT).")
        left, val = tok.split("=", 1)
        tname, key = left.rsplit(":", 1) if ":" in left else ("*", left)
        target = match_target_name(tname, target_names)
        fv = _to_float(val, key.strip())
        if fv is None:
            raise ValueError(f"--eclipse-values: Wert fuer {key.strip()!r} fehlt.")
        ref.set(target, normalize_eclipse_key(key, rx_gy), fv, "cli")
    ref.meta["cli"] = spec
    return ref


def eclipse_reference_from_dvh(dvh_map: dict, rs_ds: Optional[pydicom.Dataset],
                               targets: list, rx_gy: float) -> EclipseReference:
    """
    Eclipse-DVHs -> Referenz: je Ziel mit DVH ``tv_cm3``, ``tv_piv_cm3`` (V(Rx)),
    D98/D50/D2/Dmean/Dmin/Dmax; ein DVH der EXTERNAL-ROI (Body) liefert fuer
    alle Ziele ``piv_cm3`` = V(Rx) und ``piv50_cm3`` = V(Rx/2).  Quelle ``dvh``.
    """
    ref = EclipseReference()
    if not dvh_map:
        return ref
    names = {int(n): nm for n, nm in targets}
    body = None
    if rs_ds is not None:
        for num, name, _rt, cat in roi_table(rs_ds):
            if cat == ana.CAT_EXTERNAL and num in dvh_map:
                body = (num, name)
                break
    used = []
    for roi, name in names.items():
        dvh = dvh_map.get(roi)
        if dvh is None:
            continue
        st = dm.dvh_statistics(dvh, rx_gy)
        ref.set(name, "tv_cm3", st["total_cm3"], "dvh")
        ref.set(name, "tv_piv_cm3", st["v_rx_cm3"], "dvh")
        for k in ("d98_gy", "d50_gy", "d2_gy", "dmean_gy", "dmin_gy", "dmax_gy"):
            if st[k] is not None:
                ref.set(name, k, st[k], "dvh")
        used.append(roi)
    if body is not None:
        st = dm.dvh_statistics(dvh_map[body[0]], rx_gy)
        for name in names.values():
            ref.set(name, "piv_cm3", st["v_rx_cm3"], "dvh")
            ref.set(name, "piv50_cm3", st["v_half_rx_cm3"], "dvh")
        used.append(body[0])
    ref.meta["dvh_rois"] = sorted(set(used))
    ref.meta["body_roi"] = body[0] if body else None
    missing = [nm for roi, nm in names.items() if roi not in dvh_map]
    if missing:
        ref.notes.append("Kein Eclipse-DVH in der RTDOSE fuer: " + ", ".join(repr(m) for m in missing) + ".")
    if body is None:
        ref.notes.append("Kein Body-DVH in der RTDOSE; Eclipse-PIV/PIV50 manuell angeben "
                         "(--eclipse-values PIV=..,V50=.. oder eclipse_ref.json).")
    return ref


def derive_eclipse_values(ref: EclipseReference, target: str) -> None:
    """
    Fuellt nur Luecken (Quelle ``derived``): PIV aus CI, CI aus den
    Komponenten, GI = PIV50/PIV, HI = (D2-D98)/D50.  Sind CI und PIV beide
    gegeben und > 1 % inkonsistent, gibt es einen Hinweis.
    """
    v = lambda k: ref.value(target, k)  # noqa: E731
    tv, tvp, piv, piv50, ci = v("tv_cm3"), v("tv_piv_cm3"), v("piv_cm3"), v("piv50_cm3"), v("ci_paddick")
    if piv is None and ci and tv and tvp is not None:
        piv = tvp * tvp / (tv * ci)
        ref.set(target, "piv_cm3", piv, "derived")
        ref.notes.append(f"{target}: Eclipse-PIV {piv:.3f} cm3 aus CI {ci:.3f} abgeleitet "
                         "(PIV = (TV&PIV)^2 / (TV*CI)).")
    if tv and piv and tvp is not None:
        ci_c = tvp * tvp / (tv * piv)
        if ci is None:
            ref.set(target, "ci_paddick", ci_c, "derived")
        elif ref.get(target, "piv_cm3")["source"] != "derived" and abs(ci_c - ci) > 0.01 * abs(ci):
            ref.notes.append(f"{target}: Eclipse-CI {ci:.3f} vs. aus Komponenten {ci_c:.3f}: "
                             "inkonsistent (> 1 %).")
    if v("gi") is None and piv and piv50 is not None:
        ref.set(target, "gi", piv50 / piv, "derived")
    d2, d98, d50 = v("d2_gy"), v("d98_gy"), v("d50_gy")
    if v("hi_icru83") is None and None not in (d2, d98, d50) and d50:
        ref.set(target, "hi_icru83", (d2 - d98) / d50, "derived")


def build_eclipse_reference(rd_ds: Optional[pydicom.Dataset], rs_ds: Optional[pydicom.Dataset],
                            targets: list, rx_gy: float, json_path=None,
                            cli_string: Optional[str] = None, use_dvh: bool = True,
                            dose: Optional[dm.DoseGrid] = None) -> tuple:
    """
    Alle Quellen zusammenfuehren (Vorrang ``cli > json > dvh``), dann Luecken
    ableiten.  Liefert ``(EclipseReference, dvh_map, hinweise)``; ``dvh_map``
    ist leer, wenn die DVHs nicht zu diesem RTSTRUCT gehoeren (dann zeigt auch
    die Validierungsansicht kein Eclipse-DVH).
    """
    ref = EclipseReference()
    dvh_map = {}
    names = [n for _, n in targets]
    if use_dvh and rd_ds is not None:
        dvh_map, d_notes = dm.read_dvh_sequence(rd_ds)
        ref.notes.extend(d_notes)
        trusted = True
        if dose is not None and rs_ds is not None:
            rs_uid = str(rs_ds.get("SOPInstanceUID", ""))
            if dose.referenced_rtstruct_uid and rs_uid and dose.referenced_rtstruct_uid != rs_uid:
                trusted = False
                ref.notes.append("RTDOSE-DVHs referenzieren ein anderes Structure Set; ROI-Nummern "
                                 "nicht uebertragbar, DVH-Quelle uebersprungen.")
        if dvh_map and trusted:
            ref.merge(eclipse_reference_from_dvh(dvh_map, rs_ds, targets, rx_gy))
        elif not dvh_map:
            ref.notes.append("Keine DVHSequence in der RTDOSE (kein automatischer Eclipse-DVH-Abgleich).")
        if not trusted:
            dvh_map = {}
    if json_path:
        ref.merge(load_eclipse_reference_json(json_path, names, rx_gy))
    if cli_string:
        ref.merge(parse_eclipse_values(cli_string, names, rx_gy))
    for name in names:
        derive_eclipse_values(ref, name)
    return ref, dvh_map, list(ref.notes)


def _isnan(v) -> bool:
    try:
        return v is None or math.isnan(float(v))
    except (TypeError, ValueError):
        return True


def compare_with_eclipse(art: "DoseIndexArtifacts", ref: EclipseReference,
                         tol_pct: float = 5.0) -> dict:
    """
    Je Ziel mit Referenzwerten: Zeilen ``{key, label, tool, tool_key, eclipse,
    diff_abs, diff_pct, source, within_tol, note}`` in ``ECLIPSE_KEYS``-
    Reihenfolge (``diff = tool - eclipse``).  PIV/PIV50 werden immer gegen die
    globalen Isodosenvolumina verglichen (Eclipse-PIV = ganze Isodose).
    Zeilen ausserhalb der Toleranz werden in ``result['warnings']`` des Ziels
    eingetragen (Exit-Code bleibt unveraendert).
    """
    out = {}
    for name, tm in art.targets.items():
        vals = ref.values.get(name, {})
        if not vals:
            continue
        r = tm.result
        c, ix, dv = r["components"], r["indices"], r["dvh_stats"]
        tool = {
            "tv_cm3": ("components.tv_cm3", c["tv_cm3"]),
            "tv_piv_cm3": ("components.tv_piv_cm3", c["tv_piv_cm3"]),
            "piv_cm3": ("components.piv_global_cm3", c["piv_global_cm3"]),
            "piv50_cm3": ("components.piv50_global_cm3", c["piv50_global_cm3"]),
            "ci_paddick": ("indices.ci_paddick", ix["ci_paddick"]),
            "gi": ("indices.gi", ix["gi"]),
            "hi_icru83": ("indices.hi_icru83", ix["hi_icru83"]),
            "d98_gy": ("dvh_stats.d98_gy", dv["d98_gy"]),
            "d50_gy": ("dvh_stats.d50_gy", dv["d50_gy"]),
            "d2_gy": ("dvh_stats.d2_gy", dv["d2_gy"]),
            "dmean_gy": ("dvh_stats.dmean_gy", dv["dmean_gy"]),
            "dmin_gy": ("dvh_stats.dmin_gy", dv["dmin_gy"]),
            "dmax_gy": ("dvh_stats.dmax_gy", dv["dmax_gy"]),
        }
        rows, flagged = [], []
        for key in ECLIPSE_KEYS:
            tool_key, tool_v = tool[key]
            e = vals.get(key)
            ev = e["value"] if e else None
            src = e["source"] if e else None
            note = ""
            if ev is None and key in ("piv_cm3", "piv50_cm3"):
                note = "kein Body-DVH; manuell angeben"
            diff = diff_pct = within = None
            if ev is not None and not _isnan(tool_v):
                diff = float(tool_v) - ev
                if ev != 0:
                    diff_pct = 100.0 * diff / ev
                    within = bool(abs(diff_pct) <= tol_pct)
            rows.append({"key": key, "label": ECLIPSE_LABELS[key], "tool": None if _isnan(tool_v) else float(tool_v),
                         "tool_key": tool_key, "eclipse": ev, "diff_abs": diff, "diff_pct": diff_pct,
                         "source": src, "within_tol": within, "note": note})
            if within is False:
                flagged.append(key)
                nd = 2 if key.endswith("_gy") else 3
                r["warnings"].append(
                    f"Abgleich Eclipse: {ECLIPSE_LABELS[key]} {_fmt(tool_v, nd)} vs {_fmt(ev, nd)} "
                    f"({diff_pct:+.1f} %, Toleranz {tol_pct:g} %)"
                )
        scope_note = None
        if r["piv"].get("scope") == "component" and abs(c["piv_cm3"] - c["piv_global_cm3"]) > 1e-9:
            scope_note = ("Tool-CI/GI mit PIV-Scope component (PIV-Zeile zeigt das globale PIV); "
                          "fuer den Abgleich --piv-scope global oder --eclipse-compat")
        out[name] = {
            "rows": rows, "tol_pct": float(tol_pct),
            "n_compared": sum(1 for row in rows if row["diff_abs"] is not None),
            "n_flagged": len(flagged), "flagged": flagged,
            "sources": sorted({row["source"] for row in rows if row["source"]}),
            "piv_scope_note": scope_note,
        }
    return out


# ---------------------------------------------------------------------------
# 3. Ergebnisobjekte und Auswertung auf dem Feingitter
# ---------------------------------------------------------------------------

@dataclass
class IsodoseLevel:
    key: str
    label: str
    pct: float
    gy: float
    mask: np.ndarray
    labels: np.ndarray
    n_components: int
    volume_cm3: float
    color: tuple
    roi_name: Optional[str] = None


@dataclass
class TargetMasks:
    name: str
    roi_number: int
    color: tuple
    rt_type: str
    contours: list
    structure: dm.StructureMask
    structure_alt: dm.StructureMask
    piv: np.ndarray
    piv50: np.ndarray
    inter: np.ndarray
    under: np.ndarray
    spill: np.ndarray
    dose_samples: np.ndarray
    sample_weights: np.ndarray
    piv_info: dict
    result: dict
    helper_names: dict = field(default_factory=dict)


@dataclass
class DoseIndexArtifacts:
    grid: dm.FineGrid
    dose_fine: np.ndarray
    dose: dm.DoseGrid
    rs_ds: pydicom.Dataset
    rx_gy: float
    rx_source: str
    levels: dict                      # key -> IsodoseLevel (absteigend nach Gy)
    targets: dict                     # name -> TargetMasks
    global_result: Optional[dict]
    settings: dict
    warnings: list
    results: dict = field(default_factory=dict)
    eclipse_dvh: dict = field(default_factory=dict)   # roi_number -> dm.EclipseDVH
    eclipse: dict = field(default_factory=dict)       # name -> compare_with_eclipse-Block


def _roi_color(rs_ds: pydicom.Dataset, roi_number: int) -> tuple:
    for rc in rs_ds.ROIContourSequence:
        if int(rc.ReferencedROINumber) == int(roi_number):
            col = rc.get("ROIDisplayColor")
            if col is not None and len(col) == 3:
                return tuple(int(v) for v in col)
    return (255, 0, 0)


def _level_color(pct: float) -> tuple:
    best = min(LEVEL_COLORS, key=lambda p: abs(p - pct))
    return LEVEL_COLORS[best] if abs(best - pct) < 5 else (160, 160, 160)


def evaluate_target(name: str, roi_number: int, contours: list, color: tuple, rt_type: str,
                    grid: dm.FineGrid, dose_fine: np.ndarray, rx_gy: float, levels: dict,
                    volume_model: str, piv_scope: str) -> TargetMasks:
    """Komponenten, Indizes, DVH-Statistik und Hilfsmasken fuer EIN Ziel."""
    alt_model = "eclipse" if volume_model == "slab" else "slab"
    sm = dm.rasterize_structure(contours, grid, volume_model, name, roi_number)
    sm_alt = dm.rasterize_structure(contours, grid, alt_model, name, roi_number)
    warnings = []

    lv100, lv50 = levels["100"], levels["50"]
    if piv_scope == "component":
        piv, info = dm.scope_mask_to_structure(lv100.mask, sm.mask, lv100.labels, grid)
        piv50, info50 = dm.scope_mask_to_structure(lv50.mask, sm.mask, lv50.labels, grid)
        if info["fallback"]:
            warnings.append(
                f"Keine 100%-Isodosen-Komponente ueberlappt {name!r}; naechste Komponente "
                f"({info['fallback_distance_mm']:.1f} mm) verwendet."
            )
        info["piv50_components_used"] = info50["components_used"]
    else:
        piv, piv50 = lv100.mask, lv50.mask
        info = {"n_components_total": lv100.n_components,
                "components_used": list(range(1, lv100.n_components + 1)),
                "fallback": False, "fallback_distance_mm": None,
                "piv50_components_used": list(range(1, lv50.n_components + 1))}
    info["scope"] = piv_scope

    w = sm.weighted()
    w_alt = sm_alt.weighted()
    tv = dm.mask_volume_cm3(w, grid)
    tv_alt = dm.mask_volume_cm3(w_alt, grid)
    piv_v = dm.mask_volume_cm3(piv, grid)
    tv_piv = dm.intersection_volume_cm3(w, piv, grid)
    tv_piv_alt = dm.intersection_volume_cm3(w_alt, piv, grid)
    piv50_v = dm.mask_volume_cm3(piv50, grid)

    inter = sm.mask & piv
    under = sm.mask & ~piv
    spill = piv & ~sm.mask

    k_idx = np.nonzero(sm.mask)[0]
    samples = dose_fine[sm.mask]
    stats = dm.weighted_dose_statistics(samples, sm.slab_w[k_idx], rx_gy)
    stats_alt = dm.weighted_dose_statistics(samples, sm_alt.slab_w[k_idx], rx_gy)

    if piv_v <= 0:
        warnings.append(f"PIV leer: Rx {rx_gy:.2f} Gy liegt ueber der Maximaldosis im Gitter.")
    elif tv_piv <= 0:
        warnings.append(f"Kein Ueberlapp zwischen {name!r} und der Rx-Isodose.")
    if stats["outside_fraction"] > 0:
        warnings.append(
            f"{stats['outside_fraction'] * 100:.1f} % des Zielvolumens liegen ausserhalb "
            "des Dosisgitters (ohne Dosiswerte)."
        )
    if sm.n_gaps:
        warnings.append(f"{sm.n_gaps} z-Luecke(n) in den Konturen von {name!r} (nicht ueberbrueckt).")

    def _indices(tv_, tv_piv_, st):
        return {
            "ci_paddick": paddick_ci(tv_, piv_v, tv_piv_),
            "coverage": coverage(tv_, tv_piv_),
            "selectivity": selectivity(piv_v, tv_piv_),
            "ci_rtog": rtog_ci(tv_, piv_v),
            "dice": dice_index(tv_, piv_v, tv_piv_),
            "gi": gradient_index(piv50_v, piv_v),
            "gm_cm": gradient_measure_cm(piv50_v, piv_v),
            "hi_icru83": homogeneity_index(st["d2"], st["d98"], st["d50"]),
        }

    indices = _indices(tv, tv_piv, stats)
    indices_alt = _indices(tv_alt, tv_piv_alt, stats_alt)
    v100_cm3 = stats["v100_pct"] / 100.0 * tv if not math.isnan(stats["v100_pct"]) else float("nan")

    result = {
        "roi_number": int(roi_number),
        "n_contours": sm.n_contours,
        "n_slices": sm.n_slices,
        "slice_spacing_mm": sm.slice_spacing_mm,
        "n_holes": sm.n_holes,
        "n_gaps": sm.n_gaps,
        "outside_dose_grid_fraction": stats["outside_fraction"],
        "volumes_cm3": {
            "planimetric_slab": dm.planimetric_volume_cm3(contours, "slab"),
            "planimetric_eclipse": dm.planimetric_volume_cm3(contours, "eclipse"),
            "raster_slab": tv if volume_model == "slab" else tv_alt,
            "raster_eclipse": tv if volume_model == "eclipse" else tv_alt,
        },
        "components": {
            "tv_cm3": tv, "piv_cm3": piv_v, "piv_global_cm3": lv100.volume_cm3,
            "tv_piv_cm3": tv_piv, "piv50_cm3": piv50_v,
            "piv50_global_cm3": lv50.volume_cm3,
        },
        "indices": indices,
        "dvh_stats": {
            "d2_gy": stats["d2"], "d50_gy": stats["d50"], "d95_gy": stats["d95"],
            "d98_gy": stats["d98"], "dmin_gy": stats["dmin"], "dmax_gy": stats["dmax"],
            "dmean_gy": stats["dmean"], "v95_pct": stats["v95_pct"],
            "v100_pct": stats["v100_pct"], "v100_cm3": v100_cm3,
            "n_samples": stats["n_samples"],
        },
        "volume_models": {
            volume_model: {"tv_cm3": tv, "tv_piv_cm3": tv_piv,
                           "ci_paddick": indices["ci_paddick"],
                           "hi_icru83": indices["hi_icru83"]},
            alt_model: {"tv_cm3": tv_alt, "tv_piv_cm3": tv_piv_alt,
                        "ci_paddick": indices_alt["ci_paddick"],
                        "hi_icru83": indices_alt["hi_icru83"]},
        },
        "piv": info,
        "helper_rois": {
            "intersection": {"volume_cm3": dm.mask_volume_cm3(inter, grid)},
            "underdosed": {"volume_cm3": dm.mask_volume_cm3(under, grid)},
            "spill": {"volume_cm3": dm.mask_volume_cm3(spill, grid)},
        },
        "warnings": warnings,
    }
    return TargetMasks(
        name=name, roi_number=int(roi_number), color=color, rt_type=rt_type,
        contours=contours, structure=sm, structure_alt=sm_alt, piv=piv, piv50=piv50,
        inter=inter, under=under, spill=spill, dose_samples=samples,
        sample_weights=sm.slab_w[k_idx], piv_info=info, result=result,
    )


def evaluate_global(targets: dict, grid: dm.FineGrid, dose_fine: np.ndarray,
                    rx_gy: float, levels: dict) -> dict:
    """Union aller Ziele gegen die globalen Isodosen (nur bei > 1 Ziel)."""
    w = None
    for tm in targets.values():
        wt = tm.structure.weighted()
        w = wt if w is None else np.maximum(w, wt)
    union_mask = w > 0
    lv100, lv50 = levels["100"], levels["50"]
    tv = dm.mask_volume_cm3(w, grid)
    piv = lv100.volume_cm3
    tv_piv = dm.intersection_volume_cm3(w, lv100.mask, grid)
    piv50 = lv50.volume_cm3
    stats = dm.weighted_dose_statistics(dose_fine[union_mask], w[union_mask], rx_gy)
    return {
        "targets": list(targets.keys()),
        "components": {"tv_cm3": tv, "piv_cm3": piv, "tv_piv_cm3": tv_piv, "piv50_cm3": piv50},
        "indices": {
            "ci_paddick": paddick_ci(tv, piv, tv_piv),
            "coverage": coverage(tv, tv_piv),
            "selectivity": selectivity(piv, tv_piv),
            "ci_rtog": rtog_ci(tv, piv),
            "dice": dice_index(tv, piv, tv_piv),
            "gi": gradient_index(piv50, piv),
            "gm_cm": gradient_measure_cm(piv50, piv),
            "hi_icru83": homogeneity_index(stats["d2"], stats["d98"], stats["d50"]),
        },
        "dvh_stats": {
            "d2_gy": stats["d2"], "d50_gy": stats["d50"], "d95_gy": stats["d95"],
            "d98_gy": stats["d98"], "dmin_gy": stats["dmin"], "dmax_gy": stats["dmax"],
            "dmean_gy": stats["dmean"], "v95_pct": stats["v95_pct"],
            "v100_pct": stats["v100_pct"],
        },
    }


def evaluate_on_grid(grid: dm.FineGrid, dose_fine: np.ndarray, dose: dm.DoseGrid,
                     rs_ds: Optional[pydicom.Dataset], target_specs: list, rx_gy: float,
                     rx_source: str, level_specs: list, volume_model: str,
                     piv_scope: str, settings: Optional[dict] = None) -> DoseIndexArtifacts:
    """
    Numerischer Kern ohne I/O.  ``target_specs`` = Liste von Dicts
    ``{name, roi_number, contours, color, rt_type}``; ``level_specs`` aus
    ``parse_isodose_levels``.  Der Self-Test injiziert hier ein analytisches Feld.
    """
    ctx = _runtime.current()
    warnings = []
    levels = {}
    for i, lv in enumerate(level_specs):
        ctx.check_cancel()
        ctx.progress(i, len(level_specs), f"Isodose {lv['label']}")
        mask = dm.isodose_mask(dose_fine, lv["gy"])
        labels, n = dm.label_components(mask)
        levels[lv["key"]] = IsodoseLevel(
            key=lv["key"], label=lv["label"], pct=lv["pct"], gy=lv["gy"], mask=mask,
            labels=labels, n_components=n, volume_cm3=dm.mask_volume_cm3(mask, grid),
            color=_level_color(lv["pct"]),
        )
        if not mask.any():
            warnings.append(f"Isodose {lv['label']} ({lv['gy']:.2f} Gy) ist im Gitter leer.")
    lowest = min(levels.values(), key=lambda l: l.gy).mask
    if lowest.any() and (lowest[0].any() or lowest[-1].any() or lowest[:, 0, :].any()
                         or lowest[:, -1, :].any() or lowest[:, :, 0].any()
                         or lowest[:, :, -1].any()):
        warnings.append("Isodose beruehrt den Gitterrand (BBox oder Dosisgitter zu klein).")

    targets = {}
    for i, spec in enumerate(target_specs):
        ctx.check_cancel()
        ctx.progress(i, len(target_specs), f"Ziel {spec['name']}")
        tm = evaluate_target(
            spec["name"], spec["roi_number"], spec["contours"], spec.get("color", (255, 0, 0)),
            spec.get("rt_type", ""), grid, dose_fine, rx_gy, levels, volume_model, piv_scope,
        )
        targets[spec["name"]] = tm
    global_result = (evaluate_global(targets, grid, dose_fine, rx_gy, levels)
                     if len(targets) > 1 else None)
    return DoseIndexArtifacts(
        grid=grid, dose_fine=dose_fine, dose=dose, rs_ds=rs_ds, rx_gy=rx_gy,
        rx_source=rx_source, levels=levels, targets=targets, global_result=global_result,
        settings=dict(settings or {}), warnings=warnings,
    )


def plan_fine_grid(rs_ds: pydicom.Dataset, dose: dm.DoseGrid, target_list: list,
                   level_specs: list, grid_mm: float = 0.25,
                   restrict_z_to: Optional[np.ndarray] = None,
                   align: Optional[tuple] = None,
                   max_voxels: float = dm.MAX_FINE_VOXELS) -> tuple:
    """
    Zielkonturen und Feingitter wie ``compute_dose_indices``, ohne die Dosis
    abzutasten: ``(target_specs, grid)``.  Fuer Vorschau und Speicherschaetzung
    (``grid.n_voxels``); ``ValueError`` wie in der Rechnung (z.B. mehr als
    ``max_voxels``; ``float('inf')`` misst auch ein zu grosses Gitter).
    """
    types = ana.get_structure_type(rs_ds)
    specs, contour_sets = [], []
    for num, name in target_list:
        contours = dm.closed_planar_contours(rs_ds, num)
        if not contours:
            raise ValueError(f"Zielvolumen {name!r} (ROI {num}) hat keine CLOSED_PLANAR-Konturen.")
        specs.append({"name": name, "roi_number": int(num), "contours": contours,
                      "color": _roi_color(rs_ds, num), "rt_type": types.get(num, "")})
        contour_sets.append(contours)

    lo, hi = dm.contours_bbox(contour_sets)
    min_gy = min(lv["gy"] for lv in level_specs)
    lvl_bbox = dm.native_level_bbox(dose, min_gy, margin_voxels=2)
    if lvl_bbox is not None:
        lo, hi = np.minimum(lo, lvl_bbox[0]), np.maximum(hi, lvl_bbox[1])
    contour_z = np.array(sorted({round(float(c[0, 2]), 3) for cs in contour_sets for c in cs}))
    grid = dm.build_fine_grid(dose, lo, hi, grid_mm, contour_z=contour_z,
                              restrict_z_to=restrict_z_to, align=align, max_voxels=max_voxels)
    return specs, grid


def eclipse_compat_settings(mode: str, ct_index: Optional[dict]) -> dict:
    """
    Effektive Einstellungen von ``--eclipse-compat`` (``high`` | ``default``):
    Raster und Ursprung aus dem CT-Pixelraster, Volumenmodell eclipse, PIV
    global, linear, Feld-Isolinien ohne Vereinfachung, dazu der Hinweistext.
    ``ValueError`` ohne CT oder bei nicht quadratischen Pixeln.
    """
    if mode not in ("high", "default"):
        raise ValueError("--eclipse-compat muss 'high' oder 'default' sein.")
    if ct_index is None:
        raise ValueError("--eclipse-compat braucht den CT-Ordner (CT-Pixelraster).")
    psp = ct_index["pixel_spacing"]
    if abs(psp[0] - psp[1]) > 1e-6:
        raise ValueError("--eclipse-compat: CT-Pixel sind nicht quadratisch.")
    px = float(psp[0])
    x0, y0 = ct_index["ipp_xy"]
    if mode == "high":
        grid_mm, align = px, (x0, y0)
    else:
        grid_mm, align = 2.0 * px, (x0 + px / 2.0, y0 + px / 2.0)
    note = (f"Eclipse-kompatibel ({mode}): Raster {grid_mm:.5f} mm auf dem CT-Pixelgitter "
            f"(Ursprung x={align[0]:.4f}, y={align[1]:.4f}), Volumenmodell eclipse, PIV global, "
            "Interpolation linear, Isodosen als Feld-Isolinien ohne Vereinfachung.")
    return {"grid_mm": grid_mm, "align": align, "volume_model": "eclipse", "piv_scope": "global",
            "dose_interp": "linear", "iso_contours": "field", "simplify_mm": 0.0, "note": note}


def compute_dose_indices(rs_ds: pydicom.Dataset, dose: dm.DoseGrid, target_list: list,
                         rx_gy: float, rx_source: str, level_specs: list,
                         grid_mm: float = 0.25, dose_interp: str = "linear",
                         volume_model: str = "slab", piv_scope: str = "component",
                         restrict_z_to: Optional[np.ndarray] = None,
                         align: Optional[tuple] = None,
                         extra_settings: Optional[dict] = None) -> DoseIndexArtifacts:
    """BBox -> Feingitter (``plan_fine_grid``) -> Dosis-Sampling -> ``evaluate_on_grid``."""
    order = INTERP_ORDER[dose_interp]
    specs, grid = plan_fine_grid(rs_ds, dose, target_list, level_specs, grid_mm,
                                 restrict_z_to=restrict_z_to, align=align)
    dose_fine = dm.sample_dose_on_grid(dose, grid, order)
    settings = {
        "rx_gy": rx_gy, "rx_source": rx_source, "grid_mm": float(grid_mm),
        "dz_mm": float(grid.dz), "dose_interp": dose_interp, "volume_model": volume_model,
        "piv_scope": piv_scope,
        "grid_align": [float(v) for v in align] if align else None,
        "isodose_levels": [{"key": lv["key"], "pct": lv["pct"], "gy": lv["gy"]} for lv in level_specs],
    }
    settings.update(extra_settings or {})
    return evaluate_on_grid(grid, dose_fine, dose, rs_ds, specs, rx_gy, rx_source,
                            level_specs, volume_model, piv_scope, settings)


# ---------------------------------------------------------------------------
# 4. Report (JSON-Struktur, Konsole/TXT, CSV)
# ---------------------------------------------------------------------------

def _r(v, nd: int):
    """Rundung mit NaN/None-Durchreichung."""
    if v is None:
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return v
    if math.isnan(f) or math.isinf(f):
        return None
    return round(f, nd)


def _round_block(d: dict) -> dict:
    out = {}
    for k, v in d.items():
        if isinstance(v, dict):
            out[k] = _round_block(v)
        elif isinstance(v, (list, tuple)):
            out[k] = [_round_block(x) if isinstance(x, dict) else x for x in v]
        elif isinstance(v, (float, np.floating)):
            if k.endswith("_cm3"):
                out[k] = _r(v, 3)
            elif k.endswith("_gy") or k == "slice_spacing_mm" or k.endswith("_mm"):
                out[k] = _r(v, 2 if k.endswith("_gy") else 3)
            elif k.endswith("_pct") or k.endswith("_fraction"):
                out[k] = _r(v, 1 if k.endswith("_pct") else 4)
            else:
                out[k] = _r(v, 3)
        elif isinstance(v, (np.integer,)):
            out[k] = int(v)
        else:
            out[k] = v
    return out


def build_report(art: DoseIndexArtifacts, meta: dict, outputs: dict) -> dict:
    """Vollstaendige JSON-Struktur (gerundet) aus den Artefakten."""
    dose = art.dose
    grid = art.grid
    lo, hi = grid.bbox
    report = {
        "meta": {
            "tool": TOOL_NAME, "tool_version": TOOL_VERSION,
            "timestamp": meta.get("timestamp"),
            "case_id": meta.get("case_id"),
            "files": {"rs": meta.get("rs_path"), "rd": meta.get("rd_path"), "rp": meta.get("rp_path")},
            "rtstruct": {"label": meta.get("rs_label"), "sop_instance_uid": meta.get("rs_sop_uid"),
                         "n_rois": meta.get("rs_n_rois")},
            "frame_of_reference_uid": dose.frame_of_reference_uid,
            "rtplan": meta.get("rtplan"),
            "notes": list(meta.get("notes", [])),
            "eclipse_reference": meta.get("eclipse_reference"),
            "settings": dict(art.settings),
            "fine_grid": {"shape": list(grid.shape), "res_xy_mm": grid.res_xy, "dz_mm": grid.dz,
                          "bbox_mm": [[round(float(v), 3) for v in lo], [round(float(v), 3) for v in hi]],
                          "n_voxels": grid.n_voxels},
        },
        "dose_grid": {
            "shape": list(dose.shape), "spacing_mm": [round(float(s), 4) for s in dose.spacing],
            "origin_mm": [round(float(v), 3) for v in dose.origin], "dmax_gy": _r(dose.dmax, 3),
            "units": dose.units, "dose_type": dose.dose_type, "summation_type": dose.summation_type,
            "gfov_mode": dose.gfov_mode, "sop_instance_uid": dose.sop_instance_uid,
        },
        "isodoses": {
            lv.key: {"label": lv.label, "pct": _r(lv.pct, 3), "level_gy": _r(lv.gy, 2),
                     "volume_cm3": _r(lv.volume_cm3, 3), "r_eq_cm": _r(equivalent_radius_cm(lv.volume_cm3), 3),
                     "n_components": lv.n_components, "roi_name": lv.roi_name}
            for lv in art.levels.values()
        },
        "targets": {},
        "global": None,
        "outputs": dict(outputs),
        "warnings": list(art.warnings),
    }
    for name, tm in art.targets.items():
        block = _round_block(tm.result)
        block["helper_rois"] = {
            kind: {"roi_name": tm.helper_names.get(kind),
                   "volume_cm3": block["helper_rois"][kind]["volume_cm3"]}
            for kind in ("intersection", "underdosed", "spill")
        }
        block["eclipse"] = _round_block(art.eclipse[name]) if name in art.eclipse else None
        report["targets"][name] = block
    if art.global_result is not None:
        report["global"] = _round_block(art.global_result)
    return report


def _fmt(v, nd: int = 3, width: int = 0) -> str:
    if v is None or (isinstance(v, float) and (math.isnan(v) or math.isinf(v))):
        s = "n/a"
    else:
        s = f"{float(v):.{nd}f}"
    return s.rjust(width) if width else s


def format_report(report: dict) -> list:
    """Zeilen fuer Konsole und ``indices.txt`` (ASCII)."""
    m = report["meta"]
    s = m["settings"]
    dg = report["dose_grid"]
    L = []
    L.append("=" * 60)
    L.append(f"Dosisindex-Berechnung: {m.get('case_id')}")
    L.append("=" * 60)
    rs = m["rtstruct"]
    L.append(f"  RTSTRUCT : {Path(m['files']['rs']).name if m['files']['rs'] else '-'}  "
             f"({rs.get('n_rois')} ROIs, Label '{rs.get('label')}')")
    L.append(f"  RTDOSE   : {Path(m['files']['rd']).name if m['files']['rd'] else '-'}  "
             f"({dg['summation_type']}, {dg['units']}, "
             f"{dg['shape'][0]} x {dg['shape'][1]} x {dg['shape'][2]} Voxel @ "
             f"{dg['spacing_mm'][0]:.2f}/{dg['spacing_mm'][1]:.2f}/{dg['spacing_mm'][2]:.2f} mm, "
             f"Dmax {_fmt(dg['dmax_gy'], 2)} Gy)")
    rp = m.get("rtplan") or {}
    if m["files"].get("rp"):
        L.append(f"  RTPLAN   : {Path(m['files']['rp']).name}  (Label '{rp.get('label')}', "
                 f"{rp.get('fractions') if rp.get('fractions') is not None else '?'} Fraktion(en))")
    else:
        L.append("  RTPLAN   : -")
    L.append(f"  Rx       : {_fmt(s['rx_gy'], 2)} Gy  (Quelle: {s['rx_source']}"
             f"{', ' + m['rx_detail'] if m.get('rx_detail') else ''})")
    fg = m["fine_grid"]
    L.append(f"  Einstellungen: Gitter {s['grid_mm']:.4g} mm (z {s['dz_mm']:.2f} mm), "
             f"Interpolation {s['dose_interp']}, Volumenmodell {s['volume_model']}, "
             f"PIV-Scope {s['piv_scope']}, Isodosen-Konturen {s.get('iso_contours', 'mask')}"
             f"{' (Eclipse-kompatibel: ' + s['eclipse_compat'] + ')' if s.get('eclipse_compat') else ''}")
    L.append(f"  Feingitter: {fg['shape'][2]} x {fg['shape'][1]} x {fg['shape'][0]} Voxel "
             f"({fg['n_voxels'] / 1e6:.2f} M), BBox x[{fg['bbox_mm'][0][0]:.1f},{fg['bbox_mm'][1][0]:.1f}] "
             f"y[{fg['bbox_mm'][0][1]:.1f},{fg['bbox_mm'][1][1]:.1f}] "
             f"z[{fg['bbox_mm'][0][2]:.1f},{fg['bbox_mm'][1][2]:.1f}]")
    for n in m.get("notes", []):
        L.append(f"  Hinweis: {n}")

    L.append("")
    L.append("ISODOSEN (global)")
    L.append("-" * 60)
    L.append(f"  {'Level':<9}{'Gy':>8}{'Volumen cm3':>14}{'r_eq cm':>10}{'Komponenten':>13}   ROI")
    for key, lv in report["isodoses"].items():
        L.append(f"  {lv['label']:<9}{_fmt(lv['level_gy'], 2, 8)}{_fmt(lv['volume_cm3'], 3, 14)}"
                 f"{_fmt(lv['r_eq_cm'], 3, 10)}{lv['n_components']:>13}   {lv.get('roi_name') or '-'}")

    for name, t in report["targets"].items():
        c, ix, dv = t["components"], t["indices"], t["dvh_stats"]
        L.append("")
        L.append(f"ZIELVOLUMEN: {name}  (ROI {t['roi_number']}, {t['n_contours']} Konturen auf "
                 f"{t['n_slices']} Schichten, dz {_fmt(t.get('slice_spacing_mm'), 2)} mm, "
                 f"{t['n_holes']} Loecher)")
        L.append("-" * 60)
        L.append(f"  Komponenten [cm3]  TV {_fmt(c['tv_cm3'])} | PIV {_fmt(c['piv_cm3'])} | "
                 f"TV&PIV {_fmt(c['tv_piv_cm3'])} | PIV50 {_fmt(c['piv50_cm3'])}")
        L.append(f"  Konformitaet       CI Paddick {_fmt(ix['ci_paddick'])} "
                 f"(Coverage {_fmt(ix['coverage'])} x Selektivitaet {_fmt(ix['selectivity'])}) | "
                 f"RTOG {_fmt(ix['ci_rtog'])} | Dice {_fmt(ix['dice'])}")
        L.append(f"  Gradient           GI {_fmt(ix['gi'], 2)} | GM {_fmt(ix['gm_cm'], 2)} cm")
        L.append(f"  Homogenitaet       HI ICRU83 {_fmt(ix['hi_icru83'])} | D2 {_fmt(dv['d2_gy'], 2)} | "
                 f"D50 {_fmt(dv['d50_gy'], 2)} | D98 {_fmt(dv['d98_gy'], 2)} Gy")
        L.append(f"  Ziel-DVH           Dmin {_fmt(dv['dmin_gy'], 2)} | Dmax {_fmt(dv['dmax_gy'], 2)} | "
                 f"Dmean {_fmt(dv['dmean_gy'], 2)} | D95 {_fmt(dv['d95_gy'], 2)} Gy | "
                 f"V95 {_fmt(dv['v95_pct'], 1)} % | V100 {_fmt(dv['v100_pct'], 1)} %")
        cur = s["volume_model"]
        alt = "eclipse" if cur == "slab" else "slab"
        va = t["volume_models"].get(alt, {})
        L.append(f"  Volumenmodell      {alt}: TV {_fmt(va.get('tv_cm3'))} | "
                 f"TV&PIV {_fmt(va.get('tv_piv_cm3'))} | CI {_fmt(va.get('ci_paddick'))} | "
                 f"HI {_fmt(va.get('hi_icru83'))}")
        p = t["piv"]
        L.append(f"  PIV-Scope          {p.get('scope')}: Komponente(n) {p.get('components_used')} "
                 f"von {p.get('n_components_total')}"
                 f"{'  [Fallback: naechste Komponente]' if p.get('fallback') else ''}"
                 f"  (global PIV {_fmt(c.get('piv_global_cm3'))} cm3)")
        h = t["helper_rois"]
        L.append(f"  Hilfs-ROIs [cm3]   Schnitt {_fmt(h['intersection']['volume_cm3'])}"
                 f" ({h['intersection'].get('roi_name') or '-'}) | "
                 f"Unterdos. {_fmt(h['underdosed']['volume_cm3'])}"
                 f" ({h['underdosed'].get('roi_name') or '-'}) | "
                 f"Spill {_fmt(h['spill']['volume_cm3'])} ({h['spill'].get('roi_name') or '-'})")
        L.extend(_format_eclipse_block(t.get("eclipse"), m.get("eclipse_reference") or {}))
        if t.get("warnings"):
            for w in t["warnings"]:
                L.append(f"  ! WARNUNG !        {w}")
        else:
            L.append("  Warnungen          -")

    g = report.get("global")
    if g:
        c, ix, dv = g["components"], g["indices"], g["dvh_stats"]
        L.append("")
        L.append(f"GESAMT (Union aller Ziele: {', '.join(g['targets'])})")
        L.append("-" * 60)
        L.append(f"  Komponenten [cm3]  TV {_fmt(c['tv_cm3'])} | PIV {_fmt(c['piv_cm3'])} | "
                 f"TV&PIV {_fmt(c['tv_piv_cm3'])} | PIV50 {_fmt(c['piv50_cm3'])}")
        L.append(f"  Konformitaet       CI Paddick {_fmt(ix['ci_paddick'])} | RTOG {_fmt(ix['ci_rtog'])} | "
                 f"Dice {_fmt(ix['dice'])} | GI {_fmt(ix['gi'], 2)} | GM {_fmt(ix['gm_cm'], 2)} cm | "
                 f"HI {_fmt(ix['hi_icru83'])}")

    if report.get("warnings"):
        L.append("")
        for w in report["warnings"]:
            L.append(f"! WARNUNG ! {w}")

    o = report.get("outputs", {})
    L.append("")
    L.append("AUSGABEN")
    L.append("-" * 60)
    L.append(f"  RTSTRUCT : {o.get('rs_path') or '-'}")
    L.append(f"  JSON     : {o.get('json_path') or '-'}")
    L.append(f"  TXT      : {o.get('txt_path') or '-'}")
    L.append(f"  CSV      : {o.get('csv_path') or '-'}"
             f"{'  (+ ' + o['append_csv_path'] + ')' if o.get('append_csv_path') else ''}")
    L.append(f"  HTML     : {o.get('viz_html_path') or '-'}")
    L.append(f"  PNG      : {o.get('viz_png_path') or '-'}")
    return L


def _format_eclipse_block(ec: Optional[dict], er: dict) -> list:
    """Zeilen des Blocks 'Abgleich Eclipse' eines Ziels (ASCII)."""
    if not ec:
        return ["  Abgleich Eclipse   - (keine Referenz: kein DVH dieses Ziels in der RTDOSE, "
                "kein eclipse_ref.json, keine --eclipse-values)"]
    src = []
    if "dvh" in ec["sources"]:
        rois = ", ".join(str(r) for r in er.get("dvh_rois", []) or [])
        src.append(f"dvh = RTDOSE-DVHSequence (ROI {rois or '?'})")
    if "json" in ec["sources"]:
        src.append(f"json = {Path(er['json_path']).name if er.get('json_path') else 'eclipse_ref.json'}")
    if "cli" in ec["sources"]:
        src.append("cli = --eclipse-values")
    if "derived" in ec["sources"]:
        src.append("derived = aus Eclipse-Werten abgeleitet")
    L = [f"  Abgleich Eclipse   Toleranz {ec['tol_pct']:g} % | Quellen: {'; '.join(src) or '-'}",
         f"    {'Komponente':<18}{'Tool':>9}{'Eclipse':>9}{'Diff':>9}{'Diff %':>9}  Quelle"]
    for row in ec["rows"]:
        nd = 2 if row["key"].endswith("_gy") else 3
        mark = "  !" if row.get("within_tol") is False else ""
        note = f"  ({row['note']})" if row.get("note") else ""
        L.append(f"    {row['label']:<18}{_fmt(row['tool'], nd, 9)}{_fmt(row['eclipse'], nd, 9)}"
                 f"{_fmt(row['diff_abs'], nd, 9)}{_fmt(row['diff_pct'], 1, 9)}  "
                 f"{(row['source'] or '-'):<8}{mark}{note}")
    if ec["n_flagged"]:
        names = ", ".join(ECLIPSE_LABELS[k].split(" [")[0] for k in ec["flagged"])
        L.append(f"    {ec['n_flagged']} von {ec['n_compared']} Werten ausserhalb der Toleranz: {names}")
    else:
        L.append(f"    {ec['n_compared']} Werte innerhalb der Toleranz")
    if ec.get("piv_scope_note"):
        L.append(f"    Hinweis: {ec['piv_scope_note']}")
    return L


def csv_rows(report: dict) -> list:
    """Eine Zeile je Ziel (Dicts mit ``CSV_COLUMNS``)."""
    m, s = report["meta"], report["meta"]["settings"]
    rows = []
    for name, t in report["targets"].items():
        c, ix, dv = t["components"], t["indices"], t["dvh_stats"]
        rows.append({
            "case_id": m.get("case_id"), "target": name, "rx_gy": s["rx_gy"],
            "rx_source": s["rx_source"], "grid_mm": s["grid_mm"], "dose_interp": s["dose_interp"],
            "volume_model": s["volume_model"], "piv_scope": s["piv_scope"],
            "tv_cm3": c["tv_cm3"], "piv_cm3": c["piv_cm3"], "tv_piv_cm3": c["tv_piv_cm3"],
            "piv50_cm3": c["piv50_cm3"], "ci_paddick": ix["ci_paddick"], "ci_rtog": ix["ci_rtog"],
            "coverage": ix["coverage"], "dice": ix["dice"], "gi": ix["gi"], "gm_cm": ix["gm_cm"],
            "hi_icru83": ix["hi_icru83"], "d2_gy": dv["d2_gy"], "d50_gy": dv["d50_gy"],
            "d98_gy": dv["d98_gy"], "dmin_gy": dv["dmin_gy"], "dmax_gy": dv["dmax_gy"],
            "dmean_gy": dv["dmean_gy"], "d95_gy": dv["d95_gy"], "v95_pct": dv["v95_pct"],
            "v100_pct": dv["v100_pct"],
            "rs_file": Path(report["outputs"]["rs_path"]).name if report["outputs"].get("rs_path") else "",
            "rd_sop_uid": report["dose_grid"]["sop_instance_uid"],
            "run_timestamp": m.get("timestamp"),
        })
        ec = t.get("eclipse") or {}
        by = {row["key"]: row for row in ec.get("rows", [])}
        e = lambda k: by.get(k, {}).get("eclipse")     # noqa: E731
        d = lambda k: by.get(k, {}).get("diff_pct")    # noqa: E731
        rows[-1].update({
            "ecl_source": "+".join(ec.get("sources", [])),
            "ecl_tv_cm3": e("tv_cm3"), "ecl_tv_piv_cm3": e("tv_piv_cm3"), "ecl_piv_cm3": e("piv_cm3"),
            "ecl_piv50_cm3": e("piv50_cm3"), "ecl_ci_paddick": e("ci_paddick"), "ecl_gi": e("gi"),
            "ecl_hi_icru83": e("hi_icru83"), "ecl_d98_gy": e("d98_gy"), "ecl_d50_gy": e("d50_gy"),
            "ecl_d2_gy": e("d2_gy"),
            "d_tv_pct": d("tv_cm3"), "d_tv_piv_pct": d("tv_piv_cm3"), "d_piv_pct": d("piv_cm3"),
            "d_ci_paddick_pct": d("ci_paddick"), "d_gi_pct": d("gi"), "d_hi_icru83_pct": d("hi_icru83"),
            "d_d98_pct": d("d98_gy"),
            "ecl_tol_pct": ec.get("tol_pct"), "ecl_n_flagged": ec.get("n_flagged"),
        })
    return rows


def _csv_header_mismatch(found: list, expected: list) -> Optional[str]:
    """Meldung, wenn der vorhandene Spaltenkopf nicht ``expected`` ist, sonst None."""
    if list(found) == list(expected):
        return None
    return (f"Spaltenkopf weicht ab ({len(found)} Spalten, erwartet {len(expected)}; "
            "vermutlich aeltere Tool-Version).")


def _check_csv_header(path) -> None:
    """``ValueError`` bei Append auf eine Sammel-CSV mit anderem Spaltenkopf."""
    if path is None:
        return
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return
    with open(p, newline="", encoding="utf-8-sig") as fh:
        found = next(csv.reader(fh), [])
    msg = _csv_header_mismatch(found, CSV_COLUMNS)
    if msg:
        raise ValueError(f"Sammel-CSV {p}: {msg} Neue Datei angeben (z.B. {p.stem}_v2{p.suffix}) "
                         "oder die alte umbenennen.")


def write_csv(rows: list, path: Path, append: bool = False) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if append:
        _check_csv_header(path)
    need_header = not append or not path.exists() or path.stat().st_size == 0
    with open(path, "a" if append else "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_COLUMNS)
        if need_header:
            w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in CSV_COLUMNS})


def _json_safe(obj):
    """NaN/Inf -> None, numpy -> Python (ergaenzt analyzer._results_to_jsonable)."""
    obj = ana._results_to_jsonable(obj)
    if isinstance(obj, float) and (math.isnan(obj) or math.isinf(obj)):
        return None
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_json_safe(v) for v in obj]
    return obj


def write_json(report: dict, path: Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(_json_safe(report), fh, indent=2, ensure_ascii=False)


def write_txt(lines: list, path: Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# 5. Orchestrierung
# ---------------------------------------------------------------------------

def _say(quiet: bool, msg: str = "") -> None:
    if not quiet:
        print(msg)


def _remove_quietly(path: Path) -> None:
    try:
        Path(path).unlink(missing_ok=True)
    except OSError:
        pass                         # z.B. von einem Virenscanner gesperrt


@dataclass
class DoseRunPlan:
    """
    Ergebnis von ``prepare_dose_run``: geladene Objekte, Auswahl und die
    effektiven Einstellungen (nach ``--eclipse-compat``), noch ohne Rechnung.
    ``options`` enthaelt die uebrigen Laufoptionen, wie uebergeben.
    """
    files: dict
    case_id: str
    rs_ds: pydicom.Dataset
    rd_ds: pydicom.Dataset
    rp_ds: Optional[pydicom.Dataset]
    dose: dm.DoseGrid
    rp_refs: list
    targets: list
    rx_gy: float
    rx_source: str
    rx_detail: str
    level_specs: list
    ecl_ref: EclipseReference
    dvh_map: dict
    ct_index: Optional[dict]
    align: Optional[tuple]
    grid_mm: float
    dose_interp: str
    volume_model: str
    piv_scope: str
    iso_contours: str
    simplify_mm: float
    write_rs: bool
    warnings: list                    # echte Warnungen (! WARNUNG !)
    notes: list                       # informative Hinweise
    options: dict


def run_dose_indices(case_dir: Optional[str] = None, *, rs: Optional[str] = None,
                     rd: Optional[str] = None, rp: Optional[str] = None,
                     target: Optional[str] = None, rx: Optional[float] = None,
                     rx_pct_of_max: Optional[float] = None, isodose: str = DEFAULT_ISODOSE,
                     grid_mm: float = 0.25, dose_interp: str = "linear",
                     volume_model: str = "slab", piv_scope: str = "component",
                     output: str = "output", label: str = "_IDX", write_rs: bool = True,
                     include_target: bool = False, simplify_mm: float = 0.1,
                     transfer_syntax: str = "explicit", max_name_len: int = 64,
                     append_csv: Optional[str] = None, eclipse_compat: Optional[str] = None,
                     iso_contours: str = "mask", quiet: bool = False,
                     eclipse_ref: Optional[str] = None, eclipse_values: Optional[str] = None,
                     eclipse_tol_pct: float = 5.0, no_eclipse_dvh: bool = False,
                     no_viz: bool = False, viz_ct: bool = True) -> dict:
    """
    Kompletter Lauf: Discovery -> Berechnung -> Eclipse-Abgleich -> RS-Export
    -> Validierungsansicht -> Report/JSON/TXT/CSV.

    ``eclipse_compat`` ("high" | "default") setzt alle Parameter auf die Eclipse-
    Konventionen: Feingitter exakt auf dem CT-Pixelraster (High = 1 Pixel auf den
    Pixelzentren, Default = 2 Pixel mit Halbpixel-Versatz), Volumenmodell
    ``eclipse``, PIV global, lineare Interpolation, Isodosen-Konturen als
    Feld-Isolinien (Vertices auf den Gitterlinien) ohne Vereinfachung.

    Ablauf in drei Stufen: ``prepare_dose_run`` -> ``compute_dose_run`` ->
    ``write_dose_outputs``; ``run_dose_indices_ex`` liefert zusaetzlich die
    Artefakte.
    """
    report, _art = run_dose_indices_ex(**locals())      # alle Parameter unveraendert weiter
    return report


def run_dose_indices_ex(case_dir: Optional[str] = None, **options) -> tuple:
    """Wie ``run_dose_indices``, liefert aber ``(report, DoseIndexArtifacts)``."""
    plan = prepare_dose_run(case_dir, **options)
    art = compute_dose_run(plan)
    return write_dose_outputs(plan, art), art


@dataclass
class DoseInputs:
    """
    Ergebnis von ``load_dose_inputs``: gelesene und gegeneinander gepruefte
    Dateien, noch ohne Auswahl.  ``plan_dose_run`` veraendert es nicht; eine
    Oberflaeche kann damit beliebig oft neu planen.
    """
    files: dict
    case_id: str
    rs_ds: pydicom.Dataset
    rd_ds: pydicom.Dataset
    rp_ds: Optional[pydicom.Dataset]
    dose: dm.DoseGrid
    rp_refs: list
    ct_index: Optional[dict]
    ct_error: Optional[Exception]     # CT vorhanden, aber nicht verwendbar
    warnings: list                    # Dosisgitter und Objektverweise
    ct_warnings: list                 # fehlende CT-Schichten (validate_index_against_rs)


def load_dose_inputs(files: dict, quiet: bool = False) -> DoseInputs:
    """
    RS, RD und RP lesen, Dosisgitter und Verweise pruefen, CT-Schichtindex
    bauen; druckt den Kopf des Laufs.  Ein unbrauchbares CT wird nur vermerkt
    (``ct_error``): ob das ein Fehler ist, entscheidet ``plan_dose_run``.
    """
    def say(msg=""):
        _say(quiet, msg)

    case_id = files["case_id"]
    say(f"\nDosisindex-Berechnung fuer Case '{case_id}'")
    say(f"  RS: {files['rs'].name}\n  RD: {files['rd'].name}\n  RP: {files['rp'].name if files['rp'] else '-'}")

    rs_ds = ana.load_rtstruct(str(files["rs"]))
    rd_ds = dm.load_rtdose(str(files["rd"]))
    dose = dm.dose_grid_from_dataset(rd_ds, str(files["rd"]))
    rp_ds = dm.load_rtplan(str(files["rp"])) if files["rp"] else None
    warns = list(dose.warnings)                      # echte Warnungen (! WARNUNG !)
    warns += dm.validate_dose_against_rtstruct(dose, rs_ds, rp_ds)
    say(f"  Dosisgitter: {dose.shape[2]} x {dose.shape[1]} x {dose.shape[0]} Voxel @ "
        f"{dose.spacing[2]:.2f}/{dose.spacing[1]:.2f}/{dose.spacing[0]:.2f} mm, "
        f"Dmax {dose.dmax:.2f} Gy, FoR OK")
    rp_refs = dm.prescription_references(rp_ds) if rp_ds is not None else []

    # CT-Schichtindex (fuer RS-Export, z-Beschraenkung und Eclipse-Raster)
    ct_index, ct_error, ct_warns = None, None, []
    ct_source = files.get("ct_files") or files.get("ct_dir")
    if ct_source is not None:
        try:
            ct_index = rw.build_ct_slice_index(ct_source if files.get("ct_files") else str(ct_source))
            ct_warns = rw.validate_index_against_rs(ct_index, rs_ds)
        except _runtime.JobCancelled:
            raise
        except Exception as e:  # noqa: BLE001 - Entscheidung in plan_dose_run
            ct_index, ct_error, ct_warns = None, e, []
    return DoseInputs(
        files=files, case_id=case_id, rs_ds=rs_ds, rd_ds=rd_ds, rp_ds=rp_ds, dose=dose,
        rp_refs=rp_refs, ct_index=ct_index, ct_error=ct_error, warnings=warns,
        ct_warnings=ct_warns,
    )


def plan_dose_run(inputs: DoseInputs, *, target: Optional[str] = None, rx: Optional[float] = None,
                  rx_pct_of_max: Optional[float] = None, isodose: str = DEFAULT_ISODOSE,
                  grid_mm: float = 0.25, dose_interp: str = "linear",
                  volume_model: str = "slab", piv_scope: str = "component",
                  write_rs: bool = True, simplify_mm: float = 0.1,
                  eclipse_compat: Optional[str] = None, iso_contours: str = "mask",
                  eclipse_values: Optional[str] = None, no_eclipse_dvh: bool = False,
                  quiet: bool = False, options: Optional[dict] = None) -> DoseRunPlan:
    """
    Ziel-, Rx- und Level-Auswahl, Eclipse-Referenz, Entscheidung ueber das CT
    und ``--eclipse-compat`` auf gelesenen Eingaben; druckt Ziel, Rx und
    Eclipse-Quellen.  Eingabefehler -> ``ValueError`` (CLI-Exit 2).
    """
    def say(msg=""):
        _say(quiet, msg)

    files, rs_ds, rd_ds, dose, rp_refs = (inputs.files, inputs.rs_ds, inputs.rd_ds,
                                          inputs.dose, inputs.rp_refs)
    warns = list(inputs.warnings)
    notes = []                                       # informative Hinweise
    targets, t_notes = select_targets(rs_ds, target, rp_refs)
    notes += t_notes
    rx_gy, rx_source, rx_detail = resolve_prescription(
        rx, rx_pct_of_max, rp_refs, dose, [n for _, n in targets])
    level_specs, l_notes = parse_isodose_levels(isodose, rx_gy)
    notes += l_notes
    say(f"  Ziel(e): {', '.join(repr(n) for _, n in targets)}")
    say(f"  Rx: {rx_gy:.2f} Gy ({rx_source}: {rx_detail})")

    # Eclipse-Referenz: DVHSequence (auto) + eclipse_ref.json + --eclipse-values
    ecl_ref, dvh_map, e_notes = build_eclipse_reference(
        rd_ds, rs_ds, targets, rx_gy, json_path=files.get("eclipse_ref"),
        cli_string=eclipse_values, use_dvh=not no_eclipse_dvh, dose=dose)
    notes += e_notes
    src_bits = []
    if ecl_ref.meta.get("dvh_rois"):
        src_bits.append(f"DVHSequence ROI {', '.join(str(r) for r in ecl_ref.meta['dvh_rois'])}")
    if files.get("eclipse_ref"):
        src_bits.append(str(files["eclipse_ref"].name))
    if eclipse_values:
        src_bits.append("--eclipse-values")
    say(f"  Eclipse-Referenz: {'; '.join(src_bits) if src_bits else 'keine'}")

    # Ohne RS-Export und ohne --eclipse-compat wird das CT nicht zwingend
    # gebraucht: ein unbrauchbarer CT-Ordner ist dann eine Warnung statt eines Abbruchs.
    ct_index = inputs.ct_index
    if inputs.ct_error is not None:
        e = inputs.ct_error
        if write_rs or eclipse_compat:
            raise e
        warns.append(f"CT-Ordner nicht verwendbar ({type(e).__name__}: {e}); ohne CT-Bezug "
                     "gerechnet, Validierungsansicht ohne CT-Hintergrund.")
    elif ct_index is not None:
        warns += inputs.ct_warnings
    elif write_rs:
        notes.append("Kein CT-Ordner gefunden; RS-Export uebersprungen (nur Indizes).")
        write_rs = False

    align = None
    if eclipse_compat:
        ec = eclipse_compat_settings(eclipse_compat, ct_index)
        grid_mm, align = ec["grid_mm"], ec["align"]
        volume_model, piv_scope, dose_interp = ec["volume_model"], ec["piv_scope"], ec["dose_interp"]
        iso_contours, simplify_mm = ec["iso_contours"], ec["simplify_mm"]
        notes.append(ec["note"])
    return DoseRunPlan(
        files=files, case_id=inputs.case_id, rs_ds=rs_ds, rd_ds=rd_ds, rp_ds=inputs.rp_ds,
        dose=dose, rp_refs=rp_refs, targets=targets, rx_gy=rx_gy, rx_source=rx_source,
        rx_detail=rx_detail, level_specs=level_specs, ecl_ref=ecl_ref, dvh_map=dvh_map,
        ct_index=ct_index, align=align, grid_mm=grid_mm, dose_interp=dose_interp,
        volume_model=volume_model, piv_scope=piv_scope, iso_contours=iso_contours,
        simplify_mm=simplify_mm, write_rs=write_rs, warnings=warns, notes=notes,
        options=dict(options or {}),
    )


def prepare_dose_run(case_dir: Optional[str] = None, *, rs: Optional[str] = None,
                     rd: Optional[str] = None, rp: Optional[str] = None,
                     target: Optional[str] = None, rx: Optional[float] = None,
                     rx_pct_of_max: Optional[float] = None, isodose: str = DEFAULT_ISODOSE,
                     grid_mm: float = 0.25, dose_interp: str = "linear",
                     volume_model: str = "slab", piv_scope: str = "component",
                     output: str = "output", label: str = "_IDX", write_rs: bool = True,
                     include_target: bool = False, simplify_mm: float = 0.1,
                     transfer_syntax: str = "explicit", max_name_len: int = 64,
                     append_csv: Optional[str] = None, eclipse_compat: Optional[str] = None,
                     iso_contours: str = "mask", quiet: bool = False,
                     eclipse_ref: Optional[str] = None, eclipse_values: Optional[str] = None,
                     eclipse_tol_pct: float = 5.0, no_eclipse_dvh: bool = False,
                     no_viz: bool = False, viz_ct: bool = True,
                     files: Optional[dict] = None, out_dir: Optional[str] = None) -> DoseRunPlan:
    """
    Stufe 1: Discovery, ``load_dose_inputs`` (Laden und Pruefen der Dateien,
    CT-Schichtindex) und ``plan_dose_run`` (Ziel-, Rx- und Level-Auswahl,
    Eclipse-Referenz, ``--eclipse-compat``).  Druckt den Kopf des Laufs.
    Eingabefehler -> ``ValueError`` bzw. ``FileNotFoundError`` (CLI-Exit 2).

    ``files`` (aus ``dose_files``) ersetzt die Discovery; ``out_dir`` ist der
    genaue Ergebnisordner statt ``<output>/<case_id><label>``.
    """
    ctx = _runtime.current()
    ctx.stage("prepare", "Daten lesen und pruefen")
    options = {"output": output, "label": label, "include_target": include_target,
               "transfer_syntax": transfer_syntax, "max_name_len": max_name_len,
               "append_csv": append_csv, "eclipse_compat": eclipse_compat, "quiet": quiet,
               "eclipse_values": eclipse_values, "eclipse_tol_pct": eclipse_tol_pct,
               "no_viz": no_viz, "viz_ct": viz_ct, "out_dir": out_dir}
    if files is None:
        files = discover_dose_case(case_dir, rs, rd, rp, eclipse_ref)
    _check_csv_header(append_csv)                    # fail fast, bevor gerechnet wird
    inputs = load_dose_inputs(files, quiet)
    return plan_dose_run(
        inputs, target=target, rx=rx, rx_pct_of_max=rx_pct_of_max, isodose=isodose,
        grid_mm=grid_mm, dose_interp=dose_interp, volume_model=volume_model,
        piv_scope=piv_scope, write_rs=write_rs, simplify_mm=simplify_mm,
        eclipse_compat=eclipse_compat, iso_contours=iso_contours,
        eclipse_values=eclipse_values, no_eclipse_dvh=no_eclipse_dvh, quiet=quiet,
        options=options,
    )


def compute_dose_run(plan: DoseRunPlan) -> DoseIndexArtifacts:
    """Stufe 2: Feingitter, Dosis-Sampling, Indizes und Eclipse-Abgleich (ohne Dateien)."""
    ctx = _runtime.current()
    ctx.check_cancel()
    ctx.stage("compute", "Feingitter, Dosis und Indizes berechnen")
    o = plan.options
    restrict_z = plan.ct_index["z_values"] if plan.ct_index is not None else None
    extra = {"eclipse_compat": o["eclipse_compat"], "iso_contours": plan.iso_contours,
             "simplify_mm": float(plan.simplify_mm)}
    art = compute_dose_indices(plan.rs_ds, plan.dose, plan.targets, plan.rx_gy, plan.rx_source,
                               plan.level_specs, grid_mm=plan.grid_mm, dose_interp=plan.dose_interp,
                               volume_model=plan.volume_model, piv_scope=plan.piv_scope,
                               restrict_z_to=restrict_z, align=plan.align, extra_settings=extra)
    _say(o["quiet"], f"  Feingitter: {art.grid.shape[2]} x {art.grid.shape[1]} x {art.grid.shape[0]} Voxel "
         f"({art.grid.n_voxels / 1e6:.2f} M) @ {plan.grid_mm:g} mm, z {art.grid.dz:.2f} mm")
    art.eclipse_dvh = plan.dvh_map
    art.eclipse = compare_with_eclipse(art, plan.ecl_ref, o["eclipse_tol_pct"])
    return art


def write_dose_outputs(plan: DoseRunPlan, art: DoseIndexArtifacts) -> dict:
    """
    Stufe 3: Isodosen-RTSTRUCT, Validierungsansicht, Report und JSON/TXT/CSV;
    liefert den Report (auch in ``art.results``).  Fehler beim RS-Export und in
    der Ansicht werden Warnungen; das RS entsteht als ``.tmp`` und wird erst nach
    der Pruefung unter seinem Namen abgelegt.
    """
    ctx = _runtime.current()
    o = plan.options
    quiet, label = o["quiet"], o["label"]

    def say(msg=""):
        _say(quiet, msg)

    files, case_id, rs_ds, ct_index = plan.files, plan.case_id, plan.rs_ds, plan.ct_index
    warns = plan.warnings
    out_dir = Path(o["out_dir"]) if o.get("out_dir") else Path(o["output"]) / f"{case_id}{label}"
    out_dir.mkdir(parents=True, exist_ok=True)
    outputs = {"rs_path": None, "json_path": str(out_dir / f"{case_id}_indices.json"),
               "txt_path": str(out_dir / "indices.txt"), "csv_path": str(out_dir / "indices.csv"),
               "append_csv_path": str(o["append_csv"]) if o["append_csv"] else None,
               "viz_html_path": None, "viz_png_path": None}

    # RS-Export (Fehler duerfen die Indizes nicht verwerfen)
    specs = None
    if plan.write_rs:
        ctx.check_cancel()
        ctx.stage("rs_export", "Isodosen-RTSTRUCT schreiben")
        rs_out = out_dir / f"RS_{case_id}{label}.dcm"
        rs_tmp = rs_out.with_name(rs_out.name + ".tmp")
        try:
            desc = rw.summary_description(art)
            specs = rw.build_roi_specs(art, include_target=o["include_target"],
                                       max_name_len=o["max_name_len"], simplify_mm=plan.simplify_mm,
                                       iso_contours=plan.iso_contours)
            rw.write_isodose_rtstruct(rs_ds, ct_index, specs, rs_tmp, label=label,
                                      description=desc, transfer_syntax=o["transfer_syntax"],
                                      max_name_len=o["max_name_len"])
            problems = rw.verify_rtstruct(rs_tmp, ct_index)
            os.replace(rs_tmp, rs_out)
            if problems:
                warns.append("RS-Pruefung meldet Probleme: " + "; ".join(problems))
                say("  RS-Pruefung: PROBLEME (siehe Warnungen)")
            else:
                say(f"  RS-Pruefung: OK ({rs_out.name}, {len(specs)} ROIs)")
            outputs["rs_path"] = str(rs_out)
        except _runtime.JobCancelled:
            raise
        except Exception as e:  # noqa: BLE001 - Export darf den Lauf nicht abbrechen
            warns.append(f"RS-Export fehlgeschlagen ({type(e).__name__}: {e}).")
            say(f"  ! WARNUNG ! RS-Export fehlgeschlagen: {e}")
        finally:
            _remove_quietly(rs_tmp)

    # Validierungsansicht (validation.html + dose_overview.png); Fehler hier
    # duerfen die Indexdateien nicht entwerten -> defensiv abgefangen.
    if not o["no_viz"]:
        ctx.check_cancel()
        ctx.stage("viz", "Validierungsansicht erstellen")
        try:
            from . import dose_viz          # erst hier: laedt plotly/matplotlib
            if specs is None:          # --no-rs oder RS-Export fehlgeschlagen
                specs = rw.build_roi_specs(art, include_target=o["include_target"],
                                           max_name_len=o["max_name_len"],
                                           simplify_mm=plan.simplify_mm,
                                           iso_contours=plan.iso_contours)
            outputs.update(dose_viz.run_dose_visualization(
                art, specs, ct_index, out_dir, ct_background=o["viz_ct"],
                case_id=case_id, label=label, verbose=not quiet))
        except _runtime.JobCancelled:
            raise
        except Exception as e:  # noqa: BLE001
            warns.append(f"Visualisierung fehlgeschlagen ({type(e).__name__}: {e}); "
                         "Indexdateien bleiben gueltig.")
            say(f"  ! WARNUNG ! Visualisierung fehlgeschlagen: {e}")

    ctx.check_cancel()
    ctx.stage("reports", "Berichte schreiben")
    art.warnings = warns + art.warnings
    ecl_ref, rp_ds = plan.ecl_ref, plan.rp_ds
    meta = {
        "timestamp": _dt.datetime.now().isoformat(timespec="seconds"),
        "case_id": case_id,
        "rs_path": str(files["rs"]), "rd_path": str(files["rd"]),
        "rp_path": str(files["rp"]) if files["rp"] else None,
        "rs_label": str(rs_ds.get("StructureSetLabel", "")),
        "rs_sop_uid": str(rs_ds.get("SOPInstanceUID", "")),
        "rs_n_rois": len(rs_ds.get("StructureSetROISequence", [])),
        "rx_detail": plan.rx_detail,
        "eclipse_reference": {
            "json_path": ecl_ref.meta.get("json_path"), "cli": ecl_ref.meta.get("cli"),
            "dvh_rois": list(ecl_ref.meta.get("dvh_rois", []) or []),
            "body_roi": ecl_ref.meta.get("body_roi"), "tol_pct": float(o["eclipse_tol_pct"]),
            "sources": ecl_ref.sources(), "source": ecl_ref.meta.get("source"),
            "date": ecl_ref.meta.get("date"),
        },
        "rtplan": ({"label": str(rp_ds.get("RTPlanLabel", "")),
                    "sop_instance_uid": str(rp_ds.get("SOPInstanceUID", "")),
                    "fractions": dm.fractions_planned(rp_ds),
                    "dose_references": plan.rp_refs} if rp_ds is not None else None),
        "notes": plan.notes,
    }
    report = build_report(art, meta, outputs)
    art.results = report
    lines = format_report(report)
    say("")
    for ln in lines:
        say(ln)
    write_json(report, outputs["json_path"])
    write_txt(lines, outputs["txt_path"])
    rows = csv_rows(report)
    write_csv(rows, outputs["csv_path"], append=False)
    if o["append_csv"]:
        write_csv(rows, Path(o["append_csv"]), append=True)
    say("\nFertig.")
    return report


def list_rois(case_dir: Optional[str], rs: Optional[str], rd: Optional[str],
              rp: Optional[str]) -> int:
    """``--list``: ROI-Tabelle und Verschreibungen, dann Ende (ohne RTDOSE)."""
    files = discover_dose_case(case_dir, rs, rd, rp, need_rd=False)
    rs_ds = ana.load_rtstruct(str(files["rs"]))
    print(f"\nROIs in {files['rs'].name}:")
    print(f"  {'Nr':>4}  {'Name':<32}{'Typ':<12}{'Kategorie':<14}{'Konturen':>9}{'Vol cm3':>10}")
    for num, name, rt, cat in sorted(roi_table(rs_ds), key=lambda r: r[0]):
        cs = dm.closed_planar_contours(rs_ds, num)
        vol = ana.compute_volume(cs) if cs else 0.0
        print(f"  {num:>4}  {name[:31]:<32}{rt[:11]:<12}{cat:<14}{len(cs):>9}{vol:>10.3f}")
    if files["rp"]:
        rp_ds = dm.load_rtplan(str(files["rp"]))
        print(f"\nVerschreibungen in {files['rp'].name} (Label '{rp_ds.get('RTPlanLabel', '')}', "
              f"{dm.fractions_planned(rp_ds)} Fraktion(en)):")
        for r in dm.prescription_references(rp_ds):
            print(f"  #{r['number']} {r['description']!r} {r['structure_type']}/{r['reference_type']}: "
                  f"{r['target_prescription_dose_gy']} Gy")
    return 0


# ---------------------------------------------------------------------------
# 6. Self-Test (analytisches Phantom, keine Dateien)
# ---------------------------------------------------------------------------

def _lens_area(a: float, b: float, d: float) -> float:
    """Schnittflaeche zweier Kreise (Radien a, b, Mittelpunktsabstand d)."""
    if a <= 0 or b <= 0 or d >= a + b:
        return 0.0
    if d <= abs(a - b):
        return math.pi * min(a, b) ** 2
    t1 = a * a * math.acos((d * d + a * a - b * b) / (2 * d * a))
    t2 = b * b * math.acos((d * d + b * b - a * a) / (2 * d * b))
    t3 = 0.5 * math.sqrt((-d + a + b) * (d + a - b) * (d - a + b) * (d + a + b))
    return t1 + t2 - t3


class _Phantom:
    """
    Kugel-Ziel (Radius R bei 0) + glattes radiales Dosisfeld um CD mit
    geschlossener Form fuer alle Erwartungswerte:
      D(r) = Dmax / (1 + (r / r0)^p)   (Hill-Profil, C-unendlich),
      r(L) = r0 * (Dmax / L - 1)^(1/p);  r(Rx) = R100 = 10.5 mm, GI ~ 2.4.
    Ein Profil mit Knick am Rx-Level (z.B. quadratisch/exponentiell) wuerde
    die lineare Interpolation systematisch verzerren und den Sampler-Test
    unbrauchbar machen.
    """
    R = 10.0            # Zielradius
    DMAX, RX = 25.0, 20.0
    R100, P = 10.5, 6.0
    R0 = 10.5 / (25.0 / 20.0 - 1.0) ** (1.0 / 6.0)
    CD = np.array([2.0, 0.0, 0.0])   # Dosiszentrum (Ziel bei 0)
    DZ = 1.0

    def __init__(self):
        self.contour_z = np.arange(-9.5, 9.5 + 1e-9, self.DZ)

    def contours(self, center=(0.0, 0.0, 0.0), n: int = 360) -> list:
        cx, cy, cz = center
        return [ana._synthetic_circle(cz + z, math.sqrt(self.R ** 2 - z * z), cx, cy, n)
                for z in self.contour_z]

    def dose_at(self, pts: np.ndarray, center=None) -> np.ndarray:
        c = self.CD if center is None else np.asarray(center, float)
        r = np.linalg.norm(np.asarray(pts, float) - c, axis=1)
        return self.DMAX / (1.0 + (r / self.R0) ** self.P)

    def level_radius(self, level: float) -> float:
        if level >= self.DMAX:
            return 0.0
        return self.R0 * (self.DMAX / level - 1.0) ** (1.0 / self.P)

    def native_grid(self, extent: float = 32.0, res_xy: float = 0.5) -> dm.DoseGrid:
        zs = np.arange(-extent + 0.5, extent, self.DZ)
        xs = np.arange(-extent + res_xy / 2, extent, res_xy)
        Z, Y, X = np.meshgrid(zs, xs, xs, indexing="ij")
        pts = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
        arr = self.dose_at(pts).reshape(Z.shape).astype(np.float32)
        A = np.zeros((4, 4))
        A[3, 3] = 1.0
        A[2, 0] = self.DZ          # k -> z
        A[1, 1] = res_xy           # j -> y
        A[0, 2] = res_xy           # i -> x
        A[:3, 3] = [xs[0], xs[0], zs[0]]
        return dm.DoseGrid(array=arr, affine=A, spacing=(self.DZ, res_xy, res_xy),
                           origin=A[:3, 3].copy(), units="GY", dose_type="PHYSICAL",
                           summation_type="PLAN", dmax=float(arr.max()),
                           frame_of_reference_uid="1.2.3", sop_instance_uid="1.2.3.4")

    # --- Erwartungswerte als Scheibenstapel (exakt fuer das Slab-Modell) ---
    def _disc_radius(self, z: float) -> float:
        return math.sqrt(max(self.R ** 2 - z * z, 0.0))

    def _weights(self, model: str) -> np.ndarray:
        w = np.ones(len(self.contour_z))
        if model == "eclipse":
            w[0] = w[-1] = 0.5
        return w

    @staticmethod
    def _poly_f(n_polygon: int) -> float:
        # Polygonflaeche eines regelmaessigen n-Ecks = pi r^2 * (n/(2 pi)) sin(2 pi / n)
        return n_polygon / (2 * math.pi) * math.sin(2 * math.pi / n_polygon)

    def target_volume(self, model: str = "slab", n_polygon: int = 360) -> float:
        poly_f, w = self._poly_f(n_polygon), self._weights(model)
        return sum(wk * math.pi * self._disc_radius(z) ** 2 * poly_f
                   for wk, z in zip(w, self.contour_z)) * self.DZ / 1000.0

    def isodose_volume(self, level: float) -> float:
        rl = self.level_radius(level)
        zs = np.arange(-31.5, 32.0, self.DZ)
        return sum(math.pi * max(rl * rl - z * z, 0.0) for z in zs) * self.DZ / 1000.0

    def cumulative_volume(self, level: float, model: str = "slab", n_polygon: int = 360) -> float:
        """V(D >= level) des Ziels in cm3 (Scheibenstapel aus Kreis-Kreis-Linsen)."""
        if level <= 0:
            return self.target_volume(model, n_polygon)
        poly_f, w = self._poly_f(n_polygon), self._weights(model)
        d = float(np.linalg.norm(self.CD))
        rl = self.level_radius(level)
        tot = 0.0
        for wk, z in zip(w, self.contour_z):
            a = self._disc_radius(z) * math.sqrt(poly_f)
            b = math.sqrt(max(rl * rl - z * z, 0.0))
            tot += wk * _lens_area(a, b, d)
        return tot * self.DZ / 1000.0

    def cumulative_dvh(self, levels, model: str = "slab") -> np.ndarray:
        """Kumulatives DVH V(D >= L) fuer alle ``levels`` (cm3)."""
        return np.array([self.cumulative_volume(float(L), model) for L in levels])

    def expected(self, model: str = "slab", n_polygon: int = 360) -> dict:
        TV = self.target_volume(model, n_polygon)
        PIV, PIV50 = self.isodose_volume(self.RX), self.isodose_volume(0.5 * self.RX)
        I = self.cumulative_volume(self.RX, model, n_polygon)
        inter = lambda level: self.cumulative_volume(level, model, n_polygon)  # noqa: E731

        def d_at(frac):   # Dosis, die frac des Volumens erhaelt (bisection auf inter(L)/TV)
            lo, hi = 0.01, self.DMAX
            for _ in range(80):
                mid = 0.5 * (lo + hi)
                if inter(mid) / TV >= frac:
                    lo = mid
                else:
                    hi = mid
            return 0.5 * (lo + hi)

        levels = np.linspace(0.0, self.DMAX, 801)
        vcum = np.array([inter(L) if L > 0 else TV for L in levels])
        dmean = float(trapezoid(vcum, levels) / TV)
        d2, d50, d98, d95 = d_at(0.02), d_at(0.50), d_at(0.98), d_at(0.95)
        return {
            "tv": TV, "piv": PIV, "tv_piv": I, "piv50": PIV50,
            "ci_paddick": I * I / (TV * PIV), "coverage": I / TV, "selectivity": I / PIV,
            "ci_rtog": PIV / TV, "dice": 2 * I / (TV + PIV), "gi": PIV50 / PIV,
            "gm_cm": gradient_measure_cm(PIV50, PIV),
            "d2": d2, "d50": d50, "d98": d98, "d95": d95, "dmean": dmean,
            "hi": (d2 - d98) / d50,
            "v100_pct": 100.0 * I / TV, "v95_pct": 100.0 * inter(0.95 * self.RX) / TV,
        }


def _synthetic_rtdose_dataset(dose: dm.DoseGrid, absolute_gfov: bool = False) -> pydicom.Dataset:
    """Minimaler In-Memory-RTDOSE (uint32) aus einem DoseGrid (fuer den Loader-Test)."""
    from pydicom.dataset import Dataset, FileMetaDataset
    from pydicom.uid import ExplicitVRLittleEndian

    # DS hat hoechstens 16 Zeichen; mit dem gekuerzten Wert rechnen, damit Pixel
    # und DoseGridScaling zusammenpassen
    scaling = float(format_number_as_ds(float(dose.dmax) / 4.0e9))
    px = np.round(dose.array.astype(np.float64) / scaling).astype(np.uint32)
    ds = Dataset()
    ds.file_meta = FileMetaDataset()
    ds.file_meta.TransferSyntaxUID = ExplicitVRLittleEndian
    if PYDICOM_MAJOR < 3:   # pydicom 2.x braucht die Flags fuer pixel_array im Speicher
        ds.is_little_endian, ds.is_implicit_VR = True, False
    ds.Modality = "RTDOSE"
    ds.SOPInstanceUID = "1.2.3.4"
    ds.FrameOfReferenceUID = "1.2.3"
    ds.DoseUnits, ds.DoseType, ds.DoseSummationType = "GY", "PHYSICAL", "PLAN"
    ds.ImageOrientationPatient = [1, 0, 0, 0, 1, 0]
    ds.ImagePositionPatient = [float(v) for v in dose.origin]
    ds.PixelSpacing = [float(dose.spacing[1]), float(dose.spacing[2])]
    nk = dose.shape[0]
    offs = np.arange(nk) * dose.spacing[0]
    ds.GridFrameOffsetVector = [float(v) for v in (offs + dose.origin[2] if absolute_gfov else offs)]
    ds.NumberOfFrames = nk
    ds.Rows, ds.Columns = dose.shape[1], dose.shape[2]
    ds.SamplesPerPixel, ds.PhotometricInterpretation = 1, "MONOCHROME2"
    ds.BitsAllocated, ds.BitsStored, ds.HighBit, ds.PixelRepresentation = 32, 32, 31, 0
    ds.DoseGridScaling = format_number_as_ds(scaling)
    ds.PixelData = px.tobytes()
    return ds


def _run_self_test() -> int:
    """Analytischer Phantomtest; 0 = alle PASS, 1 = mindestens ein FAIL."""
    results = []

    def check(name, ok, detail=""):
        results.append(bool(ok))
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))

    def rel(a, b):
        return abs(a - b) / max(abs(b), 1e-12)

    print("\nDosisindex Self-Test (analytisches Kugelphantom)")
    print("-" * 70)
    ph = _Phantom()
    native = ph.native_grid()

    # 1) Loader mit synthetischem RTDOSE (relativer + absoluter GFOV)
    for absolute in (False, True):
        ds = _synthetic_rtdose_dataset(native, absolute_gfov=absolute)
        g = dm.dose_grid_from_dataset(ds)
        ok = (np.allclose(g.affine, native.affine, atol=1e-6) and rel(g.dmax, native.dmax) < 1e-6
              and g.gfov_mode == ("absolute" if absolute else "relative")
              and np.allclose(g.patient_to_index([[0.25, 0.25, 0.5]]),
                              native.patient_to_index([[0.25, 0.25, 0.5]])))
        check(f"Loader: synthetische RTDOSE, GFOV {'absolut' if absolute else 'relativ'}", ok,
              f"dmax={g.dmax:.4f}, mode={g.gfov_mode}")

    # 2) Sampler-Exaktheit auf linearem Feld (Achsen-Transposition wuerde auffallen)
    lin = dm.DoseGrid(array=native.array.copy(), affine=native.affine.copy(),
                      spacing=native.spacing, origin=native.origin.copy(), units="GY",
                      dose_type="PHYSICAL", summation_type="PLAN", dmax=0.0,
                      frame_of_reference_uid="1", sop_instance_uid="1")
    nk, nj, ni = lin.shape
    K, J, I = np.meshgrid(np.arange(nk), np.arange(nj), np.arange(ni), indexing="ij")
    P = lin.index_to_patient(np.column_stack([K.ravel(), J.ravel(), I.ravel()]))
    field_fn = lambda p: 10.0 + 0.5 * p[:, 0] + 0.25 * p[:, 1] - 0.1 * p[:, 2]
    lin.array = field_fn(P).reshape(lin.shape).astype(np.float32)
    grid_l = dm.build_fine_grid(lin, [-8, -8, -8], [8, 8, 8], 0.25)
    for name, order in (("linear", 1), ("kubisch", 3)):
        df = dm.sample_dose_on_grid(lin, grid_l, order)
        exp = np.stack([field_fn(grid_l.plane_points(k)).reshape(len(grid_l.gy), len(grid_l.gx))
                        for k in range(len(grid_l.gz))])
        err = float(np.nanmax(np.abs(df - exp)))
        check(f"Sampler {name}: lineares Feld exakt reproduziert", err < 2e-3, f"max|err|={err:.2e} Gy")
    outside = dm.sample_dose_at_points(lin, np.array([[100.0, 0.0, 0.0]]))
    check("Sampler: Punkt ausserhalb -> NaN", bool(np.isnan(outside[0])))

    # 3) Kern-Mathematik mit analytischem Feld auf dem Feingitter
    levels, _ = parse_isodose_levels("100,50", ph.RX)
    contours = ph.contours()
    cz = ph.contour_z
    for res, tol_v, tol_ci, tol_gi, tol_hi, tol_d in ((0.25, 0.005, 0.005, 0.01, 0.01, 0.1),
                                                    (1.0, 0.02, 0.02, 0.03, 0.03, 0.3)):
        lo, hi = dm.contours_bbox([contours])
        lb = dm.native_level_bbox(native, 0.5 * ph.RX)
        grid = dm.build_fine_grid(native, np.minimum(lo, lb[0]), np.maximum(hi, lb[1]), res, contour_z=cz)
        analytic = np.stack([ph.dose_at(grid.plane_points(k)).reshape(len(grid.gy), len(grid.gx))
                             for k in range(len(grid.gz))]).astype(np.float32)
        for model in ("slab", "eclipse"):
            exp = ph.expected(model)
            art = evaluate_on_grid(grid, analytic, native, None,
                                   [{"name": "Kugel", "roi_number": 1, "contours": contours}],
                                   ph.RX, "cli", levels, model, "global")
            r = art.targets["Kugel"].result
            c, ix, dv = r["components"], r["indices"], r["dvh_stats"]
            ok_v = (rel(c["tv_cm3"], exp["tv"]) < tol_v and rel(c["piv_cm3"], exp["piv"]) < tol_v
                    and rel(c["tv_piv_cm3"], exp["tv_piv"]) < tol_v and rel(c["piv50_cm3"], exp["piv50"]) < tol_v)
            check(f"Phantom {res} mm/{model}: Volumina TV/PIV/TV&PIV/PIV50 innerhalb {tol_v*100:g} %", ok_v,
                  f"TV {c['tv_cm3']:.4f}/{exp['tv']:.4f}, PIV {c['piv_cm3']:.4f}/{exp['piv']:.4f}, "
                  f"I {c['tv_piv_cm3']:.4f}/{exp['tv_piv']:.4f}, PIV50 {c['piv50_cm3']:.3f}/{exp['piv50']:.3f}")
            ok_i = (abs(ix["ci_paddick"] - exp["ci_paddick"]) < tol_ci and abs(ix["coverage"] - exp["coverage"]) < tol_ci
                    and abs(ix["dice"] - exp["dice"]) < tol_ci and rel(ix["gi"], exp["gi"]) < tol_gi
                    and abs(ix["gm_cm"] - exp["gm_cm"]) < 0.01 and abs(ix["hi_icru83"] - exp["hi"]) < tol_hi)
            check(f"Phantom {res} mm/{model}: CI/Coverage/Dice/GI/GM/HI", ok_i,
                  f"CI {ix['ci_paddick']:.4f}/{exp['ci_paddick']:.4f}, GI {ix['gi']:.3f}/{exp['gi']:.3f}, "
                  f"GM {ix['gm_cm']:.3f}/{exp['gm_cm']:.3f}, HI {ix['hi_icru83']:.4f}/{exp['hi']:.4f}")
            ok_d = (abs(dv["d2_gy"] - exp["d2"]) < tol_d and abs(dv["d50_gy"] - exp["d50"]) < tol_d
                    and abs(dv["d98_gy"] - exp["d98"]) < tol_d and abs(dv["d95_gy"] - exp["d95"]) < tol_d
                    and abs(dv["dmean_gy"] - exp["dmean"]) < tol_d)
            check(f"Phantom {res} mm/{model}: D2/D50/D95/D98/Dmean innerhalb {tol_d:g} Gy", ok_d,
                  f"D2 {dv['d2_gy']:.2f}/{exp['d2']:.2f}, D50 {dv['d50_gy']:.2f}/{exp['d50']:.2f}, "
                  f"D98 {dv['d98_gy']:.2f}/{exp['d98']:.2f}, Dmean {dv['dmean_gy']:.2f}/{exp['dmean']:.2f}")
            if res == 0.25 and model == "slab":
                art_ref = art
        # 4) volle Pipeline (Sampling vom nativen 0.5/1.0-mm-Gitter), nur slab, 0.25 mm
        if res == 0.25:
            exp = ph.expected("slab")
            art_s = compute_dose_indices.__wrapped__(grid, native, contours, levels) \
                if hasattr(compute_dose_indices, "__wrapped__") else None
            df = dm.sample_dose_on_grid(native, grid, 1)
            art_s = evaluate_on_grid(grid, df, native, None,
                                     [{"name": "Kugel", "roi_number": 1, "contours": contours}],
                                     ph.RX, "cli", levels, "slab", "global")
            r = art_s.targets["Kugel"].result
            c, ix = r["components"], r["indices"]
            check("Pipeline 0.25 mm (Sampling linear vom nativen Gitter): Volumina innerhalb 1.5 %, CI +-0.015",
                  rel(c["tv_cm3"], exp["tv"]) < 0.015 and rel(c["piv_cm3"], exp["piv"]) < 0.015
                  and abs(ix["ci_paddick"] - exp["ci_paddick"]) < 0.015 and rel(ix["gi"], exp["gi"]) < 0.03,
                  f"PIV {c['piv_cm3']:.4f}/{exp['piv']:.4f}, CI {ix['ci_paddick']:.4f}/{exp['ci_paddick']:.4f}")

    # 5) Volumenmodelle: TV_slab - TV_eclipse = 0.5*(A_erste + A_letzte)*dz
    e_s, e_e = ph.expected("slab"), ph.expected("eclipse")
    r = art_ref.targets["Kugel"].result
    diff = r["volume_models"]["slab"]["tv_cm3"] - r["volume_models"]["eclipse"]["tv_cm3"]
    check("Volumenmodell: TV_slab - TV_eclipse = halbe Endschichten (2 %)", rel(diff, e_s["tv"] - e_e["tv"]) < 0.02,
          f"{diff:.5f} vs {e_s['tv'] - e_e['tv']:.5f} cm3")
    pv_s = dm.planimetric_volume_cm3(contours, "slab")
    pv_e = dm.planimetric_volume_cm3(contours, "eclipse")
    check("Planimetrie slab/eclipse = Erwartung (1e-6)", rel(pv_s, e_s["tv"]) < 1e-6 and rel(pv_e, e_e["tv"]) < 1e-6,
          f"{pv_s:.5f}/{e_s['tv']:.5f}, {pv_e:.5f}/{e_e['tv']:.5f}")

    # 6) Ring-Ziel (Loch): XOR-Raster = Planimetrie, n_holes = 20
    ring = [c for z in cz for c in (ana._synthetic_circle(z, 8.0), ana._synthetic_circle(z, 4.0))]
    grid_r = dm.build_fine_grid(native, [-10, -10, -10], [10, 10, 10], 0.25, contour_z=cz)
    smr = dm.rasterize_structure(ring, grid_r, "slab", "Ring", 2)
    check("Ring-Ziel: Rastervolumen = Even-Odd-Planimetrie (0.5 %), 20 Loecher",
          rel(smr.volume_cm3(grid_r), dm.planimetric_volume_cm3(ring)) < 0.005 and smr.n_holes == 20,
          f"{smr.volume_cm3(grid_r):.4f} vs {dm.planimetric_volume_cm3(ring):.4f}, holes={smr.n_holes}")

    # 7) Komponenten-Scope: zweiter Hotspot bei (32,0,0) mit eigenem Ziel bei (30,0,0)
    #    (groesseres natives Gitter, damit beide Isodosen vollstaendig enthalten sind)
    wide = ph.native_grid(extent=48.0, res_xy=1.0)
    nk2, nj2, ni2 = wide.shape
    K, J, I = np.meshgrid(np.arange(nk2), np.arange(nj2), np.arange(ni2), indexing="ij")
    P = wide.index_to_patient(np.column_stack([K.ravel(), J.ravel(), I.ravel()]))
    two = np.maximum(wide.array, ph.dose_at(P, center=[32.0, 0.0, 0.0]).reshape(wide.shape).astype(np.float32))
    dose2 = dm.DoseGrid(array=two, affine=wide.affine.copy(), spacing=wide.spacing,
                        origin=wide.origin.copy(), units="GY", dose_type="PHYSICAL",
                        summation_type="PLAN", dmax=float(two.max()), frame_of_reference_uid="1",
                        sop_instance_uid="2")
    c1, c2 = ph.contours(), ph.contours(center=(30.0, 0.0, 0.0))
    lo, hi = dm.contours_bbox([c1, c2])
    lb = dm.native_level_bbox(dose2, 0.5 * ph.RX)
    grid2 = dm.build_fine_grid(dose2, np.minimum(lo, lb[0]), np.maximum(hi, lb[1]), 0.5, contour_z=cz)
    df2 = dm.sample_dose_on_grid(dose2, grid2, 1)
    specs2 = [{"name": "A", "roi_number": 1, "contours": c1}, {"name": "B", "roi_number": 2, "contours": c2}]
    art_g = evaluate_on_grid(grid2, df2, dose2, None, specs2, ph.RX, "cli", levels, "slab", "global")
    art_c = evaluate_on_grid(grid2, df2, dose2, None, specs2, ph.RX, "cli", levels, "slab", "component")
    piv_g = art_g.targets["A"].result["components"]["piv_cm3"]
    piv_c = art_c.targets["A"].result["components"]["piv_cm3"]
    exp = ph.expected("slab")
    check("Komponenten-Scope: global PIV = 2 Hotspots, component PIV = 1 Hotspot",
          rel(piv_g, 2 * exp["piv"]) < 0.03 and rel(piv_c, exp["piv"]) < 0.03
          and art_c.levels["100"].n_components == 2 and art_g.global_result is not None,
          f"global {piv_g:.3f}, component {piv_c:.3f}, Erwartung {exp['piv']:.3f}, "
          f"CI_A component {art_c.targets['A'].result['indices']['ci_paddick']:.3f}")

    # 8) Roundtrip Maske -> Konturen -> Maske
    d_iso = dm.contours_roundtrip_dice(art_ref.levels["100"].mask, art_ref.grid)
    d_tv = dm.contours_roundtrip_dice(art_ref.targets["Kugel"].structure.mask, art_ref.grid)
    check("Roundtrip mask_to_contours -> rasterize: Dice > 0.99 (Isodose) / > 0.995 (Ziel)",
          d_iso > 0.99 and d_tv > 0.995, f"Dice iso={d_iso:.4f}, Ziel={d_tv:.4f}")

    # 9) Eclipse-Abgleich: DVHSequence-Leser, Referenzmodell, Vergleich, CSV-Kopf
    from pydicom.dataset import Dataset as _DS
    from pydicom.sequence import Sequence as _Seq
    exp_s = ph.expected("slab")
    ds_dvh = _synthetic_rtdose_dataset(native)
    width, scaling = 0.1, 0.5                 # DVHData-Breite 0.2 * DVHDoseScaling 0.5 = 0.1 Gy
    n_bins = int(round(ph.DMAX / width)) + 2
    edges = np.arange(n_bins) * width
    vols = ph.cumulative_dvh(edges, "slab")

    def _dvh_item(roi, dvh_type="CUMULATIVE", units="GY"):
        it = _DS()
        it.DVHType, it.DoseUnits, it.DoseType = dvh_type, units, "PHYSICAL"
        it.DVHDoseScaling, it.DVHVolumeUnits, it.DVHNumberOfBins = scaling, "CM3", n_bins
        it.DVHData = [format_number_as_ds(float(v))
                      for pair in zip([width / scaling] * n_bins, vols.tolist()) for v in pair]
        r_ = _DS()
        r_.ReferencedROINumber, r_.DVHROIContributionType = roi, "INCLUDED"
        it.DVHReferencedROISequence = _Seq([r_])
        return it

    ds_dvh.DVHSequence = _Seq([_dvh_item(1), _dvh_item(2, dvh_type="DIFFERENTIAL"),
                               _dvh_item(3, units="RELATIVE")])
    dvh_map, dvh_notes = dm.read_dvh_sequence(ds_dvh)
    check("Eclipse-DVH: Leser nimmt das kumulative GY/CM3-DVH, ueberspringt DIFFERENTIAL/RELATIVE mit Hinweis",
          set(dvh_map) == {1} and len(dvh_notes) == 2, f"ROIs {sorted(dvh_map)}, {len(dvh_notes)} Hinweise")
    dvh1 = dvh_map[1]
    st_e = dm.dvh_statistics(dvh1, ph.RX)
    check("Eclipse-DVH: Gesamtvolumen = TV (1e-6), V(Rx) = TV&PIV (0.5 %), Bin 0.1 Gy (Skalierung angewendet)",
          rel(dvh1.total_volume_cm3, exp_s["tv"]) < 1e-6 and rel(st_e["v_rx_cm3"], exp_s["tv_piv"]) < 0.005
          and abs(dvh1.bin_width_gy - width) < 1e-9,
          f"V {dvh1.total_volume_cm3:.5f}/{exp_s['tv']:.5f}, V(Rx) {st_e['v_rx_cm3']:.5f}/{exp_s['tv_piv']:.5f}, "
          f"dbin {dvh1.bin_width_gy:.3f}")
    check("Eclipse-DVH: D2/D50/D98/Dmean aus der Kurve innerhalb 0.1 Gy, V(0) = TV, V(30 Gy) = 0",
          abs(st_e["d2_gy"] - exp_s["d2"]) < 0.1 and abs(st_e["d50_gy"] - exp_s["d50"]) < 0.1
          and abs(st_e["d98_gy"] - exp_s["d98"]) < 0.1 and abs(st_e["dmean_gy"] - exp_s["dmean"]) < 0.1
          and rel(dvh1.v_at(0.0), exp_s["tv"]) < 1e-9 and dvh1.v_at(30.0) == 0.0,
          f"D2 {st_e['d2_gy']:.2f}/{exp_s['d2']:.2f}, D50 {st_e['d50_gy']:.2f}/{exp_s['d50']:.2f}, "
          f"D98 {st_e['d98_gy']:.2f}/{exp_s['d98']:.2f}, Dmean {st_e['dmean_gy']:.2f}/{exp_s['dmean']:.2f}")

    r_ref = art_ref.targets["Kugel"].result
    c0, ix0, dv0 = r_ref["components"], r_ref["indices"], r_ref["dvh_stats"]
    good = {"*": {"TV": c0["tv_cm3"], "V20Gy": c0["tv_piv_cm3"], "PIV": c0["piv_global_cm3"],
                  "V10Gy": c0["piv50_global_cm3"], "CI": ix0["ci_paddick"], "GI": ix0["gi"],
                  "HI": ix0["hi_icru83"], "D98": dv0["d98_gy"], "D50": dv0["d50_gy"], "D2": dv0["d2_gy"],
                  "Dmean": dv0["dmean_gy"], "Dmin": dv0["dmin_gy"], "Dmax": dv0["dmax_gy"]},
            "_meta": {"source": "Self-Test"}}
    ref_ok = eclipse_reference_from_dict(good, ["Kugel"], ph.RX, "json")
    n_warn0 = len(r_ref["warnings"])
    cmp_ok = compare_with_eclipse(art_ref, ref_ok, 5.0)["Kugel"]
    check("Eclipse-Abgleich: identische Referenz -> 13 Werte verglichen, Diff 0, nichts markiert, keine Warnung",
          cmp_ok["n_compared"] == 13 and max(abs(row["diff_abs"]) for row in cmp_ok["rows"]) < 1e-9
          and cmp_ok["n_flagged"] == 0 and len(r_ref["warnings"]) == n_warn0
          and ref_ok.meta.get("source") == "Self-Test",
          f"n={cmp_ok['n_compared']}, markiert={cmp_ok['n_flagged']}")
    bad = {"*": {k: v for k, v in good["*"].items() if k not in ("CI", "GI")}}
    bad["*"]["PIV"] = good["*"]["PIV"] * 1.3
    ref_bad = eclipse_reference_from_dict(bad, ["Kugel"], ph.RX, "json")
    derive_eclipse_values(ref_bad, "Kugel")
    cmp_bad = compare_with_eclipse(art_ref, ref_bad, 5.0)["Kugel"]
    by = {row["key"]: row for row in cmp_bad["rows"]}
    check("Eclipse-Abgleich: PIV +30 % -> PIV/CI/GI markiert (CI, GI abgeleitet), TV innerhalb, Warnung im Zielblock",
          set(cmp_bad["flagged"]) == {"piv_cm3", "ci_paddick", "gi"} and by["ci_paddick"]["source"] == "derived"
          and by["gi"]["source"] == "derived" and by["tv_cm3"]["within_tol"] is True
          and any(w.startswith("Abgleich Eclipse: PIV") for w in r_ref["warnings"]),
          f"markiert={cmp_bad['flagged']}, CI {by['ci_paddick']['diff_pct']:+.1f} %")
    ref_ci = eclipse_reference_from_dict(
        {"*": {"tv_cm3": c0["tv_cm3"], "vrx": c0["tv_piv_cm3"], "ci": ix0["ci_paddick"]}}, ["Kugel"], ph.RX, "json")
    derive_eclipse_values(ref_ci, "Kugel")
    e_piv = ref_ci.get("Kugel", "piv_cm3")
    check("Eclipse-Abgleich: PIV aus CI abgeleitet = (TV&PIV)^2 / (TV*CI) = PIV des Tools (1e-9)",
          e_piv is not None and e_piv["source"] == "derived" and rel(e_piv["value"], c0["piv_global_cm3"]) < 1e-9,
          f"{e_piv['value'] if e_piv else 'n/a'} vs {c0['piv_global_cm3']:.5f}")
    ref_cli = parse_eclipse_values("TV=1.20, v20gy=1.17,PIV=1.45,V50=6.1", ["Kugel"], 20.0)
    n_err = 0
    for spec_bad, names in (("Foo:CI=1", ["Kugel"]), ("V15Gy=1", ["Kugel"]), ("CI=1", ["A", "B"]),
                            ("TV", ["Kugel"])):
        try:
            parse_eclipse_values(spec_bad, names, 20.0)
        except ValueError:
            n_err += 1
    check("Eclipse-Werte (CLI): Aliase TV/V20Gy/PIV/V50 -> kanonisch, Quelle cli; 4 ungueltige Angaben -> ValueError",
          set(ref_cli.values.get("Kugel", {})) == {"tv_cm3", "tv_piv_cm3", "piv_cm3", "piv50_cm3"}
          and all(v["source"] == "cli" for v in ref_cli.values["Kugel"].values()) and n_err == 4,
          f"{sorted(ref_cli.values.get('Kugel', {}))}, {n_err} Fehler")
    check("Sammel-CSV: abweichender Spaltenkopf wird erkannt, identischer nicht",
          _csv_header_mismatch(CSV_COLUMNS[:31], CSV_COLUMNS) is not None
          and _csv_header_mismatch(list(CSV_COLUMNS), CSV_COLUMNS) is None)

    # 10) Visualisierung: DVH-Kurve gegen weighted_dose_statistics + Smoke-Test ohne CT
    import tempfile
    from types import SimpleNamespace
    from . import dose_viz
    tm_ref = art_ref.targets["Kugel"]
    curve = dose_viz.build_dvh_curve(tm_ref.dose_samples, tm_ref.sample_weights, art_ref.grid.voxel_volume_mm3)
    st_w = dm.weighted_dose_statistics(tm_ref.dose_samples, tm_ref.sample_weights, ph.RX)
    check("DVH-Kurve: D98/D50/D2 identisch mit weighted_dose_statistics (1e-9), Gesamtvolumen = TV (1e-6)",
          abs(curve.d98_gy - st_w["d98"]) < 1e-9 and abs(curve.d50_gy - st_w["d50"]) < 1e-9
          and abs(curve.d2_gy - st_w["d2"]) < 1e-9 and rel(curve.total_cm3, c0["tv_cm3"]) < 1e-6
          and len(curve.dose_gy) <= dose_viz.VIZ_DVH_MAX_PTS + 1,
          f"D98 {curve.d98_gy:.4f}, V {curve.total_cm3:.5f}, {len(curve.dose_gy)} Punkte")
    art_ref.eclipse_dvh = {1: SimpleNamespace(dose_gy=edges, volume_cm3=vols)}
    art_ref.eclipse = compare_with_eclipse(art_ref, ref_ok, 5.0)
    specs_v = rw.build_roi_specs(art_ref)
    with tempfile.TemporaryDirectory(prefix="dose_viz_selftest_") as tmp:
        paths = dose_viz.run_dose_visualization(art_ref, specs_v, None, Path(tmp), case_id="Phantom",
                                                label="_IDX", verbose=False)
        h = Path(paths["viz_html_path"]) if paths.get("viz_html_path") else None
        p = Path(paths["viz_png_path"]) if paths.get("viz_png_path") else None
        ok_v = (h is not None and h.is_file() and h.stat().st_size > 100_000
                and p is not None and p.is_file() and p.stat().st_size > 10_000)
        detail = (f"html {h.stat().st_size / 1e6:.1f} MB, png {p.stat().st_size / 1e3:.0f} kB" if ok_v
                  else f"html={paths.get('viz_html_path')}, png={paths.get('viz_png_path')}")
    check("Visualisierung: validation.html (> 100 kB) und dose_overview.png (> 10 kB) geschrieben", ok_v, detail)

    n_fail = results.count(False)
    print("-" * 70)
    print(f"Gesamt: {'PASS' if n_fail == 0 else 'FAIL'}  "
          f"({len(results) - n_fail}/{len(results)} Pruefungen bestanden)")
    return 0 if n_fail == 0 else 1


# ---------------------------------------------------------------------------
# 7. CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="python -m dicom_file_modifier.dose_indices",
        description="Dosisindizes (Paddick-CI, GI, HI ICRU83, ...) aus RTSTRUCT + RTDOSE (+ RTPLAN) "
                    "und separate RTSTRUCT mit Isodosen- und Schnittkonturen.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Beispiele:
  %(prog)s data/<case-id>
  %(prog)s data/<case-id> --target PTV_1 --grid 0.5 --volume-model eclipse
  %(prog)s --rs RS.dcm --rd RD.dcm --rx 20 --isodose 100,80,50,12Gy
  %(prog)s data/<case-id> --list
  %(prog)s --self-test
""",
    )
    p.add_argument("case_dir", nargs="?", default=None,
                   help="Case-Ordner mit RS*.dcm, RD*.dcm, optional RP*.dcm und CT/ (fuer den RS-Export)")
    p.add_argument("--rs", dest="rs", default=None, help="RTSTRUCT-Datei (Override)")
    p.add_argument("--rd", dest="rd", default=None, help="RTDOSE-Datei (Override)")
    p.add_argument("--rp", dest="rp", default=None, help="RTPLAN-Datei (Override, fuer die Verschreibung)")
    p.add_argument("--list", action="store_true", help="ROI-Tabelle und Verschreibungen anzeigen, dann Ende")
    p.add_argument("--target", default=None,
                   help="Zielvolumen NAME[,NAME...] (exakt oder eindeutiger Teilstring); Default: PTV automatisch")
    p.add_argument("--rx", type=float, default=None, help="Verschreibungsdosis in Gy (Default: RTPLAN)")
    p.add_argument("--rx-pct-of-max", type=float, default=None,
                   help="Verschreibung als Prozent der Maximaldosis (SRS-Konvention), z.B. 80")
    p.add_argument("--isodose", default=DEFAULT_ISODOSE,
                   help=f"Isodosen-Level in %% von Rx oder absolut mit Gy (z.B. 100,80,50,12Gy); Default {DEFAULT_ISODOSE}")
    p.add_argument("--grid", type=float, default=0.25, choices=GRID_CHOICES,
                   help="In-Plane-Aufloesung des Feingitters in mm (z bleibt auf den Dosisebenen); Default 0.25")
    p.add_argument("--dose-interp", choices=("linear", "cubic"), default="linear",
                   help="Dosis-Interpolation: linear (Default) oder cubic (B-Spline)")
    p.add_argument("--volume-model", choices=("slab", "eclipse"), default="slab",
                   help="slab = jede Konturschicht volle Dicke (Default); eclipse = Endschichten halb")
    p.add_argument("--piv-scope", choices=("global", "component"), default="component",
                   help="PIV global oder nur die Isodosen-Komponente(n), die das Ziel ueberlappen (Default)")
    p.add_argument("--output", "-o", default="output", help="Basis-Ausgabeordner (Default: output)")
    p.add_argument("--label", default="_IDX", help="Suffix fuer Ausgabeordner/RS-Datei/StructureSetLabel (Default _IDX)")
    p.add_argument("--no-rs", action="store_true", help="Keine Isodosen-RTSTRUCT schreiben")
    p.add_argument("--include-target", action="store_true", help="Zielkopie(n) mit in die Isodosen-RTSTRUCT schreiben")
    p.add_argument("--simplify-mm", type=float, default=0.1,
                   help="Douglas-Peucker-Toleranz der Isodosen-Konturen in mm (< Raster/2; Default 0.1)")
    p.add_argument("--transfer-syntax", choices=("explicit", "implicit"), default="explicit",
                   help="Transfer-Syntax der RS-Datei (Default explicit VR little endian)")
    p.add_argument("--max-name-len", type=int, default=64, help="Maximale ROI-Namenslaenge (Default 64)")
    p.add_argument("--iso-contours", choices=("mask", "field"), default="mask",
                   help="Isodosen-Konturen aus der Maske (Kantenmitten, Default) oder als Feld-Isolinie "
                        "(Vertices auf den Gitterlinien wie in Eclipse)")
    p.add_argument("--eclipse-compat", choices=("high", "default"), default=None,
                   help="Alles auf Eclipse-Konventionen setzen: Feingitter auf dem CT-Pixelraster "
                        "(high = 1 Pixel auf den Pixelzentren, default = 2 Pixel mit Halbpixel-Versatz), "
                        "Volumenmodell eclipse, PIV global, linear, Feld-Isolinien ohne Vereinfachung")
    p.add_argument("--append-csv", default=None, help="Ergebniszeilen zusaetzlich an diese Sammel-CSV anhaengen")
    p.add_argument("--eclipse-ref", default=None,
                   help="JSON mit Eclipse-Referenzwerten je Ziel (Default: <case>/eclipse_ref.json, falls vorhanden)")
    p.add_argument("--eclipse-values", default=None,
                   help="Eclipse-Werte direkt: TV=1.20,VRX=1.17,PIV=1.45,V50=6.1 (auch V20Gy/V10Gy, CI, GI, "
                        "D98, ...; Dezimalpunkt; mehrere Ziele: NAME:KEY=WERT)")
    p.add_argument("--eclipse-tol-pct", type=float, default=5.0,
                   help="Toleranz des Eclipse-Abgleichs in %% (Default 5); Abweichungen darueber werden markiert")
    p.add_argument("--no-eclipse-dvh", action="store_true",
                   help="DVHSequence der RTDOSE nicht als Eclipse-Referenz benutzen")
    p.add_argument("--no-viz", action="store_true",
                   help="Keine Validierungsansicht (validation.html, dose_overview.png) schreiben")
    p.add_argument("--no-viz-ct", action="store_true",
                   help="Kein CT-Hintergrund in der Validierungsansicht (CT-Pixel werden nicht geladen)")
    p.add_argument("--self-test", action="store_true", help="Analytischer Phantomtest; Exit 0 = pass")
    return p


def main(argv: Optional[list] = None) -> int:
    args = _build_parser().parse_args(argv)
    try:
        if args.self_test:
            return _run_self_test()
        if args.list:
            return list_rois(args.case_dir, args.rs, args.rd, args.rp)
        if args.case_dir is None and (args.rs is None or args.rd is None):
            print("Fehler: Case-Ordner oder --rs und --rd angeben (siehe --help).", file=sys.stderr)
            return 2
        run_dose_indices(
            args.case_dir, rs=args.rs, rd=args.rd, rp=args.rp, target=args.target, rx=args.rx,
            rx_pct_of_max=args.rx_pct_of_max, isodose=args.isodose, grid_mm=args.grid,
            dose_interp=args.dose_interp, volume_model=args.volume_model, piv_scope=args.piv_scope,
            output=args.output, label=args.label, write_rs=not args.no_rs,
            include_target=args.include_target, simplify_mm=args.simplify_mm,
            transfer_syntax=args.transfer_syntax, max_name_len=args.max_name_len,
            append_csv=args.append_csv, eclipse_compat=args.eclipse_compat,
            iso_contours=args.iso_contours, eclipse_ref=args.eclipse_ref,
            eclipse_values=args.eclipse_values, eclipse_tol_pct=args.eclipse_tol_pct,
            no_eclipse_dvh=args.no_eclipse_dvh, no_viz=args.no_viz, viz_ct=not args.no_viz_ct,
        )
    except (FileNotFoundError, ValueError, KeyError) as e:
        print(f"\nFehler: {e}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
