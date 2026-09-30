"""
Workflow Dosisindizes (``dose_indices``): Paddick-CI, GI, HI, Eclipse-Abgleich,
Isodosen-RTSTRUCT und Validierungsansicht.

``inspect`` liest RS, RD, RP und den CT-Schichtindex einmal
(``dose_indices.load_dose_inputs``); ``preview`` plant darauf mit denselben
Funktionen wie der Lauf (``plan_dose_run``, ``plan_fine_grid``) und meldet
Ziele, Rx, Isodosen mit ROI-Namen, das Feingitter mit Speicherschaetzung und
alle Befunde.  CLI-Entsprechung: ``dfm dose-indices``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from .. import analyzer as ana
from .. import dose as dm
from .. import dose_indices as di
from .. import rtstruct_writer as rw
from .fields import SettingsBase, setting
from .issues import Issue, has_errors, issue_from_exception
from .outputs import folder_name_problem
from .results import JobResult, blocked, run_guarded
from .selection import CaseSelection
from .sysinfo import available_memory_bytes, format_bytes

WORKFLOW = "dose"
TOOL = "dose-indices"
NON_SETTINGS = {TOOL: {"case_dir", "rs", "rd", "rp", "list", "output", "eclipse_ref", "self_test",
                       "help"}}
# Von --eclipse-compat festgelegte Felder
ECLIPSE_LOCKED = ("grid_mm", "dose_interp", "volume_model", "piv_scope", "iso_contours",
                  "simplify_mm")


def _locked_by_mode(s, info) -> Optional[str]:
    return f"durch Eclipse-kompatibel ({s.eclipse_compat}) festgelegt" if s.eclipse_compat else None


def _needs_rs(s, info) -> Optional[str]:
    return None if s.write_rs else "nur mit Isodosen-RTSTRUCT"


def _needs_viz(s, info) -> Optional[str]:
    return None if s.viz else "nur mit Validierungsansicht"


@dataclass
class DoseIndexSettings(SettingsBase):
    WORKFLOW = WORKFLOW

    # -- Basis ---------------------------------------------------------------
    target: Optional[list[str]] = setting(
        None, label="Zielvolumen", kind="names", cli_flag="--target",
        help="ROI-Namen (exakt oder eindeutiger Teil); leer = alle PTVs, eine passende "
             "RTPLAN-Verschreibung grenzt auf ein PTV ein.")
    rx: Optional[float] = setting(
        None, label="Verschreibung", unit="Gy", min=0.01, cli_flag="--rx",
        help="Verschreibungsdosis; leer = aus dem RTPLAN.")
    rx_pct_of_max: Optional[float] = setting(
        None, label="Verschreibung in % von Dmax", unit="%", min=0.01, max=100.0,
        cli_flag="--rx-pct-of-max", help="SRS-Konvention, z.B. 80 (statt einer Dosis in Gy).")
    isodose: str = setting(
        di.DEFAULT_ISODOSE, label="Isodosen", kind="spec", cli_flag="--isodose",
        help="Level in % von Rx oder absolut mit Gy, z.B. 100,80,50,12Gy; 100 und 50 werden "
             "immer ausgewertet.")
    eclipse_compat: Optional[str] = setting(
        None, label="Eclipse-kompatibel", choices=("high", "default"), cli_flag="--eclipse-compat",
        help="Alle Parameter auf Eclipse-Konventionen: Raster auf dem CT-Pixelgitter (high = "
             "1 Pixel, default = 2 Pixel), Volumenmodell eclipse, PIV global, linear, "
             "Feld-Isolinien ohne Vereinfachung.")
    write_rs: bool = setting(
        True, label="Isodosen-RTSTRUCT schreiben", cli_flag="--no-rs", cli_invert=True,
        help="Separate RTSTRUCT mit Isodosen sowie Schnitt-, Unterdosierungs- und Spill-Konturen "
             "(braucht das CT).")
    viz: bool = setting(
        True, label="Validierungsansicht", cli_flag="--no-viz", cli_invert=True,
        help="validation.html und dose_overview.png schreiben.")
    label: str = setting(
        "_IDX", label="Kennung", cli_flag="--label",
        help="Suffix fuer Ergebnisordner, RS-Datei und StructureSetLabel.")
    # -- Erweitert -----------------------------------------------------------
    grid_mm: float = setting(
        0.25, label="Feingitter", unit="mm", level="advanced", choices=di.GRID_CHOICES,
        cli_flag="--grid", enabled_if=_locked_by_mode,
        help="In-Plane-Aufloesung; z bleibt auf den Dosisebenen.")
    dose_interp: str = setting(
        "linear", label="Dosis-Interpolation", level="advanced", choices=("linear", "cubic"),
        cli_flag="--dose-interp", enabled_if=_locked_by_mode,
        help="linear oder cubic (B-Spline).")
    volume_model: str = setting(
        "slab", label="Volumenmodell", level="advanced", choices=("slab", "eclipse"),
        cli_flag="--volume-model", enabled_if=_locked_by_mode,
        help="slab = jede Konturschicht volle Dicke; eclipse = Endschichten halb.")
    piv_scope: str = setting(
        "component", label="PIV-Bereich", level="advanced", choices=("global", "component"),
        cli_flag="--piv-scope", enabled_if=_locked_by_mode,
        help="component = nur Isodosen-Anteile, die das Ziel ueberlappen; global = ganze Isodose.")
    iso_contours: str = setting(
        "mask", label="Isodosen-Konturen", level="advanced", choices=("mask", "field"),
        cli_flag="--iso-contours", enabled_if=_locked_by_mode,
        help="mask = aus der Maske (Kantenmitten); field = Isolinie des Dosisfelds wie Eclipse.")
    eclipse_values: Optional[str] = setting(
        None, label="Eclipse-Werte", level="advanced", kind="spec", cli_flag="--eclipse-values",
        help="Referenzwerte, z.B. TV=1.20,PIV=1.45,V50=6.1 (Dezimalpunkt; mehrere Ziele: "
             "NAME:KEY=WERT).")
    eclipse_tol_pct: float = setting(
        5.0, label="Toleranz Eclipse-Abgleich", unit="%", level="advanced", min=0.0,
        cli_flag="--eclipse-tol-pct", help="Abweichungen darueber werden markiert.")
    eclipse_dvh: bool = setting(
        True, label="Eclipse-DVH der RTDOSE verwenden", level="advanced",
        cli_flag="--no-eclipse-dvh", cli_invert=True,
        help="DVHSequence der RTDOSE automatisch als Eclipse-Referenz nutzen.")
    append_csv: Optional[str] = setting(
        None, label="Sammel-CSV", level="advanced", kind="path", cli_flag="--append-csv",
        help="Ergebniszeilen zusaetzlich an diese CSV anhaengen (Kopfzeile wird vorab geprueft).")
    include_target: bool = setting(
        False, label="Zielkopie ins RTSTRUCT", level="advanced", cli_flag="--include-target",
        enabled_if=_needs_rs, help="Originalkonturen der Ziele mit in die Isodosen-RTSTRUCT.")
    viz_ct: bool = setting(
        True, label="CT-Hintergrund", level="advanced", cli_flag="--no-viz-ct", cli_invert=True,
        enabled_if=_needs_viz, help="CT-Schichten in der Validierungsansicht zeigen.")
    # -- Experte -------------------------------------------------------------
    simplify_mm: float = setting(
        0.1, label="Konturvereinfachung", unit="mm", level="expert", min=0.0,
        cli_flag="--simplify-mm", enabled_if=_locked_by_mode,
        help="Douglas-Peucker-Toleranz der Isodosen-Konturen (kleiner als das halbe Raster).")
    transfer_syntax: str = setting(
        "explicit", label="Transfer-Syntax", level="expert", choices=("explicit", "implicit"),
        cli_flag="--transfer-syntax", enabled_if=_needs_rs,
        help="Explicit oder Implicit VR Little Endian fuer die RS-Datei.")
    max_name_len: int = setting(
        64, label="Max. ROI-Namenslaenge", level="expert", min=8, max=64,
        cli_flag="--max-name-len", enabled_if=_needs_rs, help="DICOM erlaubt bis zu 64 Zeichen.")

    def _check(self, info) -> list:
        issues = []
        if self.rx is not None and self.rx_pct_of_max is not None:
            issues.append(Issue("error", "DOSE.RX_CONFLICT",
                                "Verschreibung entweder in Gy oder in % von Dmax angeben, nicht beides.",
                                field="rx"))
        problem = folder_name_problem(f"RS_x{self.label}") if self.label else None
        if problem:
            issues.append(Issue("error", "DOSE.LABEL_INVALID", f"Kennung {self.label!r}: {problem}",
                                field="label"))
        if self.target and any("," in t for t in self.target):
            issues.append(Issue("error", "DOSE.TARGET_COMMA",
                                "ROI-Namen mit Komma sind als Zielvolumen nicht waehlbar.",
                                field="target"))
        return issues


Settings = DoseIndexSettings


def plan_kwargs(s: DoseIndexSettings) -> dict:
    """Argumente fuer ``dose_indices.plan_dose_run``."""
    return {"target": ",".join(s.target) if s.target else None, "rx": s.rx,
            "rx_pct_of_max": s.rx_pct_of_max, "isodose": s.isodose, "grid_mm": s.grid_mm,
            "dose_interp": s.dose_interp, "volume_model": s.volume_model,
            "piv_scope": s.piv_scope, "write_rs": s.write_rs, "simplify_mm": s.simplify_mm,
            "eclipse_compat": s.eclipse_compat, "iso_contours": s.iso_contours,
            "eclipse_values": s.eclipse_values, "no_eclipse_dvh": not s.eclipse_dvh}


def run_kwargs(s: DoseIndexSettings) -> dict:
    """Argumente fuer ``run_dose_indices`` (wie ``dose_indices.main`` sie aus der CLI baut)."""
    kw = plan_kwargs(s)
    kw.update({"label": s.label, "include_target": s.include_target,
               "transfer_syntax": s.transfer_syntax, "max_name_len": s.max_name_len,
               "append_csv": s.append_csv, "eclipse_tol_pct": s.eclipse_tol_pct,
               "no_viz": not s.viz, "viz_ct": s.viz_ct})
    return kw


# ---------------------------------------------------------------------------
# Inspektion
# ---------------------------------------------------------------------------

@dataclass
class DoseCaseInfo:
    selection: CaseSelection
    inputs: Optional[di.DoseInputs] = field(default=None, repr=False)
    rois: list = field(default_factory=list)        # [{number, name, rt_type, category}]
    targets: dict = field(default_factory=dict)     # {default, reason, ptvs, all, notes}
    prescriptions: dict = field(default_factory=dict)   # {refs, default, reason, dmax_gy}
    dose: dict = field(default_factory=dict)
    labels: dict = field(default_factory=dict)      # {rs: StructureSetLabel, rp: RTPlanLabel}
    dvh: dict = field(default_factory=dict)         # {available, trusted, rois, body, notes}
    ct: Optional[dict] = None
    rs_export: dict = field(default_factory=dict)   # {possible, reason}
    eclipse_compat: dict = field(default_factory=dict)  # {high|default: {grid_mm, align} | {reason}}
    issues: list = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return self.inputs is not None and not has_errors(self.issues)


def _files_for(sel: CaseSelection) -> dict:
    return di.dose_files(sel.rs, sel.rd, sel.rp, ct=sel.ct_files or None,
                         eclipse_ref=sel.eclipse_ref, case_id=sel.case_id or None)


def inspect(selection: CaseSelection) -> DoseCaseInfo:
    """Liest RS, RD, RP und CT-Header; Kandidaten fuer Ziel und Rx, DVH, CT; still."""
    info = DoseCaseInfo(selection=selection)
    try:
        inputs = di.load_dose_inputs(_files_for(selection), quiet=True)
    except Exception as exc:  # noqa: BLE001 - Befund statt Absturz
        info.issues.append(issue_from_exception(exc))
        return info
    info.inputs = inputs
    rs_ds, dose = inputs.rs_ds, inputs.dose
    info.issues += [Issue("warning", "DOSE.INPUT_WARNING", w) for w in inputs.warnings]

    table = di.roi_table(rs_ds)
    info.rois = [{"number": n, "name": nm, "rt_type": rt, "category": cat} for n, nm, rt, cat in table]
    cand = di.target_candidates(rs_ds, inputs.rp_refs)
    info.targets = {"default": [nm for _, nm in cand["default"]], "reason": cand["reason"],
                    "ptvs": [r[1] for r in cand["ptvs"]], "all": [r[1] for r in cand["targets"]],
                    "notes": list(cand["notes"])}
    if not cand["ptvs"]:
        info.issues.append(Issue("warning", "DOSE.NO_PTV", "Kein PTV erkannt.",
                                 hint_de="Zielvolumen von Hand waehlen.", field="target"))
    rxc = di.rx_candidates(inputs.rp_refs, dose, info.targets["default"])
    info.prescriptions = dict(rxc)
    if rxc["default"] is None:
        info.issues.append(Issue("warning", "DOSE.RX_MANUAL",
                                 f"Verschreibung nicht automatisch bestimmbar ({rxc['reason']}).",
                                 hint_de="Dosis in Gy oder % von Dmax angeben.", field="rx"))
    info.dose = {"file": Path(inputs.files["rd"]).name, "shape": list(dose.shape),
                 "spacing_mm": [float(v) for v in dose.spacing], "dmax_gy": float(dose.dmax),
                 "units": dose.units, "dose_type": dose.dose_type,
                 "summation_type": dose.summation_type}
    rp_ds = inputs.rp_ds
    info.labels = {"rs": str(rs_ds.get("StructureSetLabel", "") or ""),
                   "rp": str(rp_ds.get("RTPlanLabel", "") or "") if rp_ds is not None else ""}

    dvh_map, dvh_notes = dm.read_dvh_sequence(inputs.rd_ds)
    rs_uid = str(rs_ds.get("SOPInstanceUID", ""))
    trusted = not (dose.referenced_rtstruct_uid and rs_uid and dose.referenced_rtstruct_uid != rs_uid)
    names = {n: nm for n, nm, _, _ in table}
    body = next((nm for n, nm, _, cat in table if cat == ana.CAT_EXTERNAL and n in dvh_map), None)
    info.dvh = {"available": bool(dvh_map), "trusted": trusted,
                "rois": [names.get(n, f"ROI {n}") for n in sorted(dvh_map)], "body": body,
                "notes": list(dvh_notes)}

    if inputs.ct_index is not None:
        ci = inputs.ct_index
        info.ct = {"n_slices": ci["n_slices"], "pixel_spacing": list(ci["pixel_spacing"]),
                   "ipp_xy": list(ci["ipp_xy"]), "z_min": float(ci["z_values"][0]),
                   "z_max": float(ci["z_values"][-1]), "series_uid": ci["series_uid"]}
        info.issues += [Issue("warning", "DOSE.CT_SLICES", w) for w in inputs.ct_warnings]
        info.rs_export = {"possible": True, "reason": ""}
    elif inputs.ct_error is not None:
        e = inputs.ct_error
        reason = f"CT nicht verwendbar ({type(e).__name__}: {e})"
        info.rs_export = {"possible": False, "reason": reason}
        info.issues.append(Issue("warning", "DOSE.CT_UNUSABLE", reason + ".",
                                 hint_de="Ohne CT nur Indizes: RS-Export und Eclipse-kompatibel aus.",
                                 field="write_rs"))
    else:
        info.rs_export = {"possible": False, "reason": "kein CT gewaehlt"}
    for mode in ("high", "default"):
        try:
            ec = di.eclipse_compat_settings(mode, inputs.ct_index)
            info.eclipse_compat[mode] = {"grid_mm": ec["grid_mm"], "align": list(ec["align"])}
        except ValueError as exc:
            info.eclipse_compat[mode] = {"reason": str(exc)}
    return info


# ---------------------------------------------------------------------------
# Vorschau
# ---------------------------------------------------------------------------

def estimate_memory_bytes(n_voxels: int, n_levels: int, n_targets: int, dose: dm.DoseGrid) -> int:
    """
    Grobe Spitze der Rechnung: Dosis float32 und Arbeitskopien (~12 B), je Level
    Maske + Komponenten-Labels (5 B), je Ziel Masken und Gewichte (~7 B) pro
    Feingitter-Voxel, dazu das native Dosisgitter.
    """
    per_voxel = 12 + 5 * n_levels + 7 * n_targets
    return int(n_voxels * per_voxel + 4 * dose.array.size)


def resolve_effective(settings: DoseIndexSettings, info: DoseCaseInfo) -> tuple:
    """
    ``(werte, gesperrt, issues)``: effektive Werte der vom Modus betroffenen
    Felder, die gesperrten Felder mit Grund und Befunde (Eclipse-kompatibel
    ohne verwendbares CT ist ein Fehler).
    """
    values = {k: getattr(settings, k) for k in ECLIPSE_LOCKED}
    values["align"] = None
    locked, issues = {}, []
    mode = settings.eclipse_compat
    if mode:
        ec = info.eclipse_compat.get(mode, {"reason": "Fall nicht gelesen"})
        if "reason" in ec:
            issues.append(Issue("error", "DOSE.ECLIPSE_COMPAT_UNAVAILABLE",
                                f"Eclipse-kompatibel ({mode}) nicht moeglich: {ec['reason']}",
                                field="eclipse_compat"))
        else:
            values.update(grid_mm=ec["grid_mm"], align=tuple(ec["align"]), volume_model="eclipse",
                          piv_scope="global", dose_interp="linear", iso_contours="field",
                          simplify_mm=0.0)
        locked = {k: _locked_by_mode(settings, info) for k in ECLIPSE_LOCKED}
    return values, locked, issues


@dataclass
class DosePreview:
    ok: bool
    issues: list = field(default_factory=list)
    effective: dict = field(default_factory=dict)   # grid_mm, align, dose_interp, ..., write_rs
    locked: dict = field(default_factory=dict)      # Feld -> Grund
    targets: list = field(default_factory=list)     # Zielnamen
    rx: dict = field(default_factory=dict)          # {gy, source, detail}
    levels: list = field(default_factory=list)      # [{key, label, pct, gy}]
    roi_names: list = field(default_factory=list)   # [(kind, name, bezug)] der Isodosen-RTSTRUCT
    grid: dict = field(default_factory=dict)        # {shape, n_voxels, res_xy, dz, limit, memory_bytes}
    eclipse: dict = field(default_factory=dict)     # {sources, keys: {ziel: [schluessel]}}


_FAILED = object()


def _step(issues: list, fld: str, fn, *args):
    """``fn(*args)``; jede Ausnahme -> Befund am Feld ``fld``, Rueckgabe ``_FAILED``."""
    try:
        return fn(*args)
    except Exception as exc:  # noqa: BLE001 - die Vorschau wirft nie
        issues.append(issue_from_exception(exc, field=fld))
        return _FAILED


def preview(info: DoseCaseInfo, settings: DoseIndexSettings) -> DosePreview:
    """Pruefung vor dem Start mit den Funktionen des Laufs (nichts wird geschrieben)."""
    issues = list(info.issues) + settings.validate(info)
    values, locked, mode_issues = resolve_effective(settings, info)
    issues += mode_issues
    pv = DosePreview(ok=False, issues=issues, effective=values, locked=locked)
    if not info.ok or has_errors(issues):
        return pv
    inputs = info.inputs
    kw = plan_kwargs(settings)
    # Einzelschritte zuerst, damit ein Fehler am richtigen Feld steht
    picked = _step(issues, "target", di.select_targets, inputs.rs_ds, kw["target"], inputs.rp_refs)
    if picked is _FAILED:
        return pv
    names = [nm for _, nm in picked[0]]
    rx = _step(issues, "rx", di.resolve_prescription, settings.rx, settings.rx_pct_of_max,
               inputs.rp_refs, inputs.dose, names)
    if rx is _FAILED:
        return pv
    _step(issues, "isodose", di.parse_isodose_levels, settings.isodose, rx[0])
    if settings.eclipse_values:
        _step(issues, "eclipse_values", di.parse_eclipse_values, settings.eclipse_values, names, rx[0])
    if settings.append_csv:
        _step(issues, "append_csv", di._check_csv_header, settings.append_csv)
    if has_errors(issues):
        return pv
    try:
        plan = di.plan_dose_run(inputs, **kw, quiet=True)
    except Exception as exc:  # noqa: BLE001 - z.B. CT unbrauchbar bei RS-Export
        issues.append(issue_from_exception(exc, field="write_rs" if inputs.ct_error else ""))
        return pv

    issues += [Issue("warning", "DOSE.WARNING", w) for w in plan.warnings[len(inputs.warnings):]]
    issues += [Issue("info", "DOSE.NOTE", n) for n in plan.notes]
    pv.effective = {"grid_mm": plan.grid_mm, "align": plan.align, "dose_interp": plan.dose_interp,
                    "volume_model": plan.volume_model, "piv_scope": plan.piv_scope,
                    "iso_contours": plan.iso_contours, "simplify_mm": plan.simplify_mm,
                    "write_rs": plan.write_rs}
    pv.targets = [nm for _, nm in plan.targets]
    pv.rx = {"gy": plan.rx_gy, "source": plan.rx_source, "detail": plan.rx_detail}
    pv.levels = [{k: lv[k] for k in ("key", "label", "pct", "gy")} for lv in plan.level_specs]
    if plan.write_rs:
        pv.roi_names = rw.planned_roi_names(plan.level_specs, pv.targets, settings.include_target,
                                            settings.max_name_len)
    pv.eclipse = {"sources": plan.ecl_ref.sources(),
                  "keys": {t: sorted(plan.ecl_ref.values.get(t, {})) for t in pv.targets}}
    restrict = plan.ct_index["z_values"] if plan.ct_index is not None else None
    try:
        _, grid = di.plan_fine_grid(plan.rs_ds, plan.dose, plan.targets, plan.level_specs,
                                    plan.grid_mm, restrict_z_to=restrict, align=plan.align,
                                    max_voxels=float("inf"))
    except ValueError as exc:
        issues.append(issue_from_exception(exc, field="grid_mm"))
        return pv
    mem = estimate_memory_bytes(grid.n_voxels, len(plan.level_specs), len(plan.targets), plan.dose)
    pv.grid = {"shape": list(grid.shape), "n_voxels": grid.n_voxels, "res_xy": grid.res_xy,
               "dz": grid.dz, "limit": dm.MAX_FINE_VOXELS, "memory_bytes": mem}
    if grid.n_voxels > dm.MAX_FINE_VOXELS:
        issues.append(Issue("error", "DOSE.GRID_TOO_LARGE",
                            f"Feingitter zu gross ({grid.n_voxels / 1e6:.1f} M Voxel, erlaubt "
                            f"{dm.MAX_FINE_VOXELS / 1e6:.0f} M).",
                            hint_de="Groeberes Raster waehlen (0,5 oder 1,0 mm).", field="grid_mm"))
    avail = available_memory_bytes()
    if avail is not None and mem > avail:
        issues.append(Issue("warning", "SYSTEM.MEMORY_LOW",
                            f"Geschaetzter Speicherbedarf {format_bytes(mem)}, frei "
                            f"{format_bytes(avail)}.",
                            hint_de="Andere Programme schliessen oder ein groeberes Raster waehlen.",
                            field="grid_mm"))
    pv.ok = not has_errors(issues)
    return pv


# ---------------------------------------------------------------------------
# Befehl und Lauf
# ---------------------------------------------------------------------------

def default_folder(settings: DoseIndexSettings, selection: CaseSelection) -> str:
    return f"{selection.case_id}{settings.label}"


def command(settings: DoseIndexSettings, selection: CaseSelection, out_dir) -> list:
    """
    ``dfm dose-indices`` mit denselben Eingaben: der Fallordner, wenn die
    Auswahl der CLI-Discovery entspricht, sonst ``--rs/--rd/--rp/--eclipse-ref``.
    Die CLI legt ``<case_id><label>`` unter ``--output`` an: Heisst ``out_dir``
    so, ist ``--output`` sein Elternordner (exakt); sonst ``out_dir`` selbst
    (Ergebnis eine Ebene tiefer, ueberschreibt aber nichts).  Was die CLI nicht
    ausdruecken kann, nennt ``command_gaps``.
    """
    argv = ["dfm", TOOL]
    if _case_dir_matches(selection):
        argv.append(str(selection.case_dir))
    else:
        argv += ["--rs", str(selection.rs), "--rd", str(selection.rd)]
        if selection.rp:
            argv += ["--rp", str(selection.rp)]
        if selection.eclipse_ref:
            argv += ["--eclipse-ref", str(selection.eclipse_ref)]
    argv += settings.to_argv(TOOL)
    out_dir = Path(out_dir)
    exact = out_dir.name == default_folder(settings, selection)
    argv += ["--output", str(out_dir.parent if exact else out_dir)]
    return [argv]


def _key(sel: CaseSelection) -> tuple:
    def r(p):
        return str(Path(p).resolve()) if p else None
    return (r(sel.rs), r(sel.rd), r(sel.rp), r(sel.eclipse_ref), tuple(r(f) for f in sel.ct_files))


def _case_dir_matches(sel: CaseSelection) -> bool:
    """Trifft die CLI-Discovery im Fallordner genau diese Auswahl?"""
    from .selection import for_dose
    if not sel.case_dir:
        return False
    try:
        return _key(for_dose(sel.case_dir)) == _key(sel)
    except (OSError, ValueError):
        return False


def _cli_key(sel: CaseSelection) -> Optional[tuple]:
    """Auswahl, die der Befehl von ``command`` tatsaechlich verwenden wuerde."""
    from .selection import for_dose
    if _case_dir_matches(sel):
        return _key(sel)
    try:
        return _key(for_dose(None, sel.rs, sel.rd, sel.rp, sel.eclipse_ref))
    except (OSError, ValueError):
        return None


def command_gaps(selection: CaseSelection) -> list:
    """Unterschiede zwischen der Auswahl und dem, was der CLI-Befehl verwenden wuerde."""
    got = _cli_key(selection)
    want = _key(selection)
    if got is None or got == want:
        return []
    gaps = []
    if got[3] != want[3]:
        gaps.append("Eclipse-Referenz (die CLI nimmt eclipse_ref.json neben der RTDOSE automatisch)")
    if got[4] != want[4]:
        gaps.append("CT-Schichten (die CLI sucht das CT nur im Unterordner CT neben der RTDOSE)")
    return gaps


def _issues_from_report(report: dict) -> list:
    issues = [Issue("warning", "DOSE.WARNING", w) for w in report.get("warnings", [])]
    for name, block in report.get("targets", {}).items():
        for w in block.get("warnings", []):
            code = "DOSE.ECLIPSE_TOLERANCE" if w.startswith("Abgleich Eclipse") else "DOSE.TARGET"
            issues.append(Issue("warning", code, f"{name}: {w}"))
    issues += [Issue("info", "DOSE.NOTE", n) for n in report.get("meta", {}).get("notes", [])]
    return issues


def _summary(report: dict) -> dict:
    out = {"rx_gy": None, "n_voxels": report["meta"]["fine_grid"]["n_voxels"], "targets": {}}
    s = report["meta"]["settings"]
    out["rx_gy"] = s.get("rx_gy")
    for name, b in report["targets"].items():
        c, ix, dv = b["components"], b["indices"], b["dvh_stats"]
        ec = b.get("eclipse") or {}
        out["targets"][name] = {
            "tv_cm3": c["tv_cm3"], "piv_cm3": c["piv_cm3"], "tv_piv_cm3": c["tv_piv_cm3"],
            "piv50_cm3": c["piv50_cm3"], "ci_paddick": ix["ci_paddick"], "gi": ix["gi"],
            "hi_icru83": ix["hi_icru83"], "d98_gy": dv["d98_gy"], "d50_gy": dv["d50_gy"],
            "d2_gy": dv["d2_gy"], "eclipse_n_flagged": ec.get("n_flagged"),
            "eclipse_n_compared": ec.get("n_compared"),
        }
    return out


def planned_stages(settings: DoseIndexSettings, selection: CaseSelection) -> list:
    """Stufen-Schluessel des Laufs in Reihenfolge (fuer ``stage i/n`` im Worker)."""
    return (["prepare", "compute"] + (["rs_export"] if settings.write_rs else [])
            + (["viz"] if settings.viz else []) + ["reports"])


def job_hooks(settings: DoseIndexSettings) -> tuple:
    """
    Fuer Jobs: ``(settings fuer den Lauf, Befunde vorab, Nachtrag nach dem Commit)``.
    Die Sammel-CSV wird erst geschrieben, wenn der Ergebnisordner steht; ihre
    Kopfzeile wird vorher geprueft.  Die angehaengten Zeilen sind die der
    ``indices.csv`` im Ergebnisordner (dieselben, die die CLI anhaengt).
    """
    if not settings.append_csv:
        return settings, [], None
    issues = []
    try:
        di._check_csv_header(settings.append_csv)
    except Exception as exc:  # noqa: BLE001 - z.B. gesperrt oder alte Kopfzeile
        issues.append(issue_from_exception(exc, field="append_csv"))
    target = Path(settings.append_csv)

    def after_commit(result_dir: Path) -> dict:
        import csv
        with open(Path(result_dir) / "indices.csv", newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        di.write_csv(rows, target, append=True)
        return {"append_csv": str(target.resolve())}

    return settings.replace(append_csv=None), issues, after_commit


def run(settings: DoseIndexSettings, selection: CaseSelection, out_dir) -> JobResult:
    """Kompletter Lauf nach ``out_dir`` (Konsolentext wie die CLI; im Worker: Protokoll)."""
    out_dir = Path(out_dir)
    cmd = command(settings, selection, out_dir)
    pre = settings.validate()
    if has_errors(pre):
        return blocked(WORKFLOW, out_dir, pre, cmd)

    def body():
        files = _files_for(selection)
        report, _art = di.run_dose_indices_ex(None, files=files, out_dir=str(out_dir),
                                              output=str(out_dir.parent), **run_kwargs(settings))
        roles = {"json_path": "json", "txt_path": "txt", "csv_path": "csv", "rs_path": "rs",
                 "viz_html_path": "viz_html", "viz_png_path": "viz_png",
                 "append_csv_path": "append_csv"}
        outputs = {roles.get(k, k): v for k, v in report["outputs"].items() if v}
        return outputs, _issues_from_report(report), _summary(report)

    gaps = command_gaps(selection)
    if gaps:
        pre.append(Issue("info", "DOSE.COMMAND_APPROX",
                         "Der CLI-Befehl gibt die Dateiauswahl nicht genau wieder: " + "; ".join(gaps)))
    return run_guarded(WORKFLOW, out_dir, body, pre, cmd)
