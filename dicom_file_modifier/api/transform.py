"""
Workflow Transformation: starre Bewegung von CT + RTSTRUCT (``case_modifier``)
oder nur des CT (``modifier``), je nachdem ob die Auswahl ein RTSTRUCT enthaelt.

``inspect`` prueft nur Header (Geometrie, Orientierung, FoR, Marker);
``preview`` rechnet die Matrix, den Drehpunkt, das transformierte RTSTRUCT
(Clipping) und die Speicherschaetzung, ohne etwas zu schreiben.
CLI-Entsprechung: ``dfm case-transform`` bzw. ``dfm ct-transform``.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Optional

import numpy as np

from .. import case_modifier as cm
from .. import modifier as mod
from .fields import SettingsBase, setting
from .issues import Issue, UserInputError, has_errors, issue_from_exception
from .results import JobResult, blocked, run_guarded
from .selection import CaseSelection
from .sysinfo import available_memory_bytes, format_bytes

WORKFLOW = "transform"
CASE_TOOL, CT_TOOL = "case-transform", "ct-transform"
TOOLS = (CASE_TOOL, CT_TOOL)
NON_SETTINGS = {
    CASE_TOOL: {"case_dir", "output", "rs_override", "list_markers", "non_interactive", "dry_run",
                "self_test", "help"},
    CT_TOOL: {"ct_dir", "output", "save_viz", "help"},
}
VIZ_HTML = "visualization_3d.html"       # CT-only-Ansicht im Ergebnisordner


def _needs_rs(s, info) -> Optional[str]:
    return None if info is None or info.has_rs else "nur mit RTSTRUCT"


def _needs_resample(s, info) -> Optional[str]:
    return None if s.method == "resample" else "nur bei Neuabtastung (resample)"


def _needs_case_viz(s, info) -> Optional[str]:
    return _needs_rs(s, info) or (None if s.viz else "nur mit Vorher/Nachher-Ansicht")


@dataclass
class TransformSettings(SettingsBase):
    WORKFLOW = WORKFLOW

    tx: float = setting(0.0, label="Verschiebung X", unit="mm", cli_flag="--tx",
                        help="Positiv = nach links (Patient).")
    ty: float = setting(0.0, label="Verschiebung Y", unit="mm", cli_flag="--ty",
                        help="Positiv = nach posterior.")
    tz: float = setting(0.0, label="Verschiebung Z", unit="mm", cli_flag="--tz",
                        help="Positiv = nach superior (kranial).")
    rx: float = setting(0.0, label="Rotation um die Links-Rechts-Achse", unit="Grad", min=-180.0,
                        max=180.0, cli_flag="--rx",
                        help="Um die X-Achse; positiv nach der Rechte-Hand-Regel: +Y (posterior) dreht "
                             "nach +Z (superior). In Rueckenlage Pitch.")
    ry: float = setting(0.0, label="Rotation um die anterior-posteriore Achse", unit="Grad",
                        min=-180.0, max=180.0, cli_flag="--ry",
                        help="Um die Y-Achse; positiv: +Z (superior) dreht nach +X (links). "
                             "In Rueckenlage Yaw.")
    rz: float = setting(0.0, label="Rotation um die Kopf-Fuss-Achse", unit="Grad", min=-180.0,
                        max=180.0, cli_flag="--rz",
                        help="Um die Z-Achse; positiv: +X (links) dreht nach +Y (posterior). "
                             "In Rueckenlage Roll.")
    center: str = setting(
        "volume", label="Drehpunkt", kind="spec", cli_flag="--center", cli_default=None,
        tools=(CASE_TOOL,), enabled_if=_needs_rs,
        help="volume (Volumenmitte), marker:NAME (POINT-Marker) oder x,y,z (LPS, mm).")
    method: str = setting(
        "resample", label="Methode", choices=("resample", "metadata"), cli_flag="--method",
        help="resample = Neuabtastung, achsiale Schichten; metadata = Pixel unveraendert, nur "
             "Lage-Tags (schraege Schichten, nicht jedes TPS nimmt sie an).")
    label: str = setting(
        "_RB", label="Kennung", cli_flag="--label", tools=(CASE_TOOL,),
        help="Suffix des Ergebnisordners; mit RTSTRUCT auch fuer RS-Datei, StructureSetLabel und "
             "SeriesDescription.")
    new_frame_of_reference: bool = setting(
        False, label="Neue FrameOfReferenceUID", cli_flag="--new-frame-of-reference",
        tools=(CASE_TOOL,), enabled_if=_needs_rs,
        help="Neue FoR fuer CT und RS: das TPS legt dann keine alten Plaene/Dosen darueber. "
             "Ohne: gleiche FoR wie das Original (Warnung).")
    order: int = setting(
        1, label="Interpolation", choices=(0, 1, 3), level="advanced", cli_flag="--order",
        enabled_if=_needs_resample,
        help="0 = naechster Nachbar (HU unveraendert), 1 = linear, 3 = kubisch.")
    verify: bool = setting(
        False, label="Centroide pruefen", level="advanced", cli_flag="--verify",
        tools=(CASE_TOOL,), enabled_if=_needs_rs,
        help="Nach dem Schreiben je ROI den Schwerpunkt gegen T pruefen (Schwelle 1e-3 mm).")
    viz: bool = setting(
        True, label="Vorher/Nachher-Ansicht", level="advanced", cli_flag="--no-viz",
        cli_invert=True, help="Plots und 3D-Ansicht der Transformation schreiben.")
    viz_ct_surface: bool = setting(
        False, label="CT-Koerperoberflaeche in 3D", level="expert", cli_flag="--viz-ct-surface",
        tools=(CASE_TOOL,), enabled_if=_needs_case_viz,
        help="Zusaetzlich die Koerperoberflaeche aus dem CT extrahieren (rechenintensiv).")

    def _check(self, info) -> list:
        issues = []
        if info is not None and not info.has_rs and self.center.strip().lower() != "volume":
            issues.append(Issue("error", "TRANSFORM.CENTER_NEEDS_RS",
                                "Ohne RTSTRUCT dreht das CT immer um die Volumenmitte.",
                                hint_de="Drehpunkt 'volume' waehlen oder ein RTSTRUCT hinzunehmen.",
                                field="center"))
        if not any((self.tx, self.ty, self.tz, self.rx, self.ry, self.rz)):
            issues.append(Issue("warning", "TRANSFORM.IDENTITY", "Keine Bewegung eingestellt.",
                                hint_de="Das Ergebnis ist eine Kopie mit neuen UIDs.", field="tx"))
        return issues


Settings = TransformSettings


def describe_motion(s: TransformSettings) -> str:
    """Klartext, z.B. ``10 mm nach links · 15° um die Kopf-Fuss-Achse`` (LPS, patientenbezogen)."""
    parts = []
    for v, pos, neg in ((s.tx, "links", "rechts"), (s.ty, "posterior", "anterior"),
                        (s.tz, "superior", "inferior")):
        if v:
            parts.append(f"{abs(v):g} mm nach {pos if v > 0 else neg}")
    for v, axis in ((s.rx, "die Links-Rechts-Achse"), (s.ry, "die anterior-posteriore Achse"),
                    (s.rz, "die Kopf-Fuss-Achse")):
        if v:
            parts.append(f"{v:+g}° um {axis}")
    return " · ".join(parts) or "keine Bewegung"


def estimate_memory_bytes(shape, method: str, order: int, ct_surface: bool = False) -> int:
    """
    Grobe Spitze: resample ~20 B/Voxel (Pixeldaten, Cache, HU float32, float64-Kopie,
    Ergebnis) bzw. 28 B bei kubisch, dazu die Koordinaten eines 20-Schichten-Blocks;
    metadata ~6 B/Voxel (Pixeldaten und Kopien beim Schreiben).
    """
    nz, ny, nx = (int(v) for v in shape)
    n = nz * ny * nx
    if method == "metadata":
        est = 6 * n + (4 * n if ct_surface else 0)
    else:
        est = (20 if order <= 1 else 28) * n + 96 * min(20, nz) * ny * nx
    return int(est + (8 * n if ct_surface else 0))


@dataclass
class TransformCaseInfo:
    selection: CaseSelection
    has_rs: bool = False
    case_id: str = ""
    preflight: Optional[cm.CasePreflight] = field(default=None, repr=False)   # mit RTSTRUCT
    ct_headers: Optional[list] = field(default=None, repr=False)              # nur CT
    geom: Optional[dict] = field(default=None, repr=False)
    ct: dict = field(default_factory=dict)
    rtstruct: dict = field(default_factory=dict)
    markers: list = field(default_factory=list)    # [{name, position_mm}]
    issues: list = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return self.geom is not None and not has_errors(self.issues)

    @property
    def volume_center(self) -> Optional[np.ndarray]:
        return mod.volume_center(self.geom) if self.geom is not None else None


def inspect(selection: CaseSelection) -> TransformCaseInfo:
    """CT-Header, Geometrie, Orientierung und (mit RTSTRUCT) FoR und Marker; still."""
    info = TransformCaseInfo(selection=selection, has_rs=bool(selection.rs),
                             case_id=selection.case_id)
    if not selection.ct_files:
        info.issues.append(Issue("error", "TRANSFORM.NO_CT", "Keine CT-Serie gewaehlt.",
                                 field="ct_files"))
        return info
    try:
        if selection.rs:
            pre = cm.preflight_case(selection.case_dir, rs_override=selection.rs,
                                    ct_files=selection.ct_files, case_id=selection.case_id or None,
                                    siblings=selection.related, quiet=True)
            info.preflight, info.ct_headers, info.geom = pre, pre.ct_headers, pre.geom
            info.case_id = pre.case_id
            info.issues += pre.issues
            info.markers = [{"name": n, "position_mm": [float(v) for v in p]} for n, p in pre.markers]
            rs = pre.rs_ds
            info.rtstruct = {"file": pre.rs_path.name, "label": str(rs.get("StructureSetLabel", "")),
                             "n_rois": len(rs.get("StructureSetROISequence", []))}
        else:
            headers = mod.load_ct_headers(selection.ct_files)
            mod.validate_ct_geometry(headers)
            info.ct_headers, info.geom = headers, mod.extract_geometry(headers)
    except Exception as exc:  # noqa: BLE001 - Befund statt Absturz
        info.issues.append(issue_from_exception(exc))
        return info
    g, h0 = info.geom, info.ct_headers[0]
    c = mod.volume_center(g)
    info.ct = {"n_slices": len(info.ct_headers), "shape": list(g["shape"]),
               "spacing_mm": {"dz": g["dz"], "dr": g["dr"], "dc": g["dc"]},
               "volume_center_mm": [float(v) for v in c],
               "series_description": str(h0.get("SeriesDescription", "")),
               "series_uid": str(h0.get("SeriesInstanceUID", "")),
               "frame_of_reference_uid": str(h0.get("FrameOfReferenceUID", ""))}
    return info


@dataclass
class TransformPreview:
    ok: bool
    tool: str
    issues: list = field(default_factory=list)
    T: list = field(default_factory=list)             # 4x4, Patient -> Patient
    center_mm: list = field(default_factory=list)
    center_label: str = ""
    drehpunkt_mm: Optional[list] = None               # Marker-ROI im Ergebnis (nur mit RS)
    description: str = ""
    planned: dict = field(default_factory=dict)       # Rolle -> Pfad relativ zum Ergebnisordner
    clipping: list = field(default_factory=list)      # [{roi, n_outside, n_total, fraction}]
    memory_bytes: int = 0


def preview(info: TransformCaseInfo, settings: TransformSettings) -> TransformPreview:
    """Matrix, Drehpunkt, Clipping und Speicher, ohne zu schreiben."""
    tool = CASE_TOOL if info.has_rs else CT_TOOL
    issues = list(info.issues) + settings.validate(info)
    pv = TransformPreview(ok=False, tool=tool, issues=issues, description=describe_motion(settings))
    if not info.ok or has_errors(issues):
        return pv
    s = settings
    vol_c = info.volume_center
    pv.memory_bytes = estimate_memory_bytes(info.geom["shape"], s.method, s.order,
                                            s.viz and s.viz_ct_surface and info.has_rs)
    try:
        cm.validate_label(s.label, info.case_id)                # auch der Ordnername ohne RTSTRUCT
    except UserInputError as exc:
        issues.append(exc.issue)
        return pv
    if info.has_rs:
        pre = info.preflight
        try:
            center, center_label = cm.resolve_center(s.center, pre.rs_ds, vol_c, interactive=False)
            plan = cm.plan_transform(pre, s.tx, s.ty, s.tz, s.rx, s.ry, s.rz, out_dir=".",
                                     method=s.method, order=s.order, label=s.label, center=center,
                                     center_label=center_label,
                                     new_frame_of_reference=s.new_frame_of_reference, quiet=True)
        except UserInputError as exc:
            issues.append(exc.issue)
            return pv
        except Exception as exc:  # noqa: BLE001 - z.B. unbekannter Marker
            issue = issue_from_exception(exc, field="center")
            if not s.center.strip().lower().startswith(("marker:", "volume")):
                issue = replace(issue, hint_de="Koordinate als x,y,z in mm (LPS) angeben, z. B. 12.5,-3,0.")
            issues.append(issue)
            return pv
        issues += plan.issues
        pv.T, pv.center_mm, pv.center_label = plan.T.tolist(), plan.center.tolist(), plan.center_label
        pv.drehpunkt_mm = plan.drehpunkt_pos.tolist()
        pv.clipping = cm._clipping_dicts(plan.clipping)
        pv.planned = {"ct_dir": "CT", "rs": f"RS{s.label}.dcm"}
    else:
        T = mod.build_rigid_transform(s.rx, s.ry, s.rz, s.tx, s.ty, s.tz, vol_c)
        pv.T, pv.center_mm, pv.center_label = T.tolist(), vol_c.tolist(), "Volumenmitte"
        pv.planned = {"ct_dir": "CT"}
        if s.viz:
            pv.planned["viz_html"] = VIZ_HTML
    avail = available_memory_bytes()
    if avail is not None and pv.memory_bytes > avail:
        issues.append(Issue("warning", "SYSTEM.MEMORY_LOW",
                            f"Geschaetzter Speicherbedarf {format_bytes(pv.memory_bytes)}, frei "
                            f"{format_bytes(avail)}.",
                            hint_de="Andere Programme schliessen oder die Methode 'metadata' waehlen.",
                            field="method"))
    pv.ok = not has_errors(issues)
    return pv


def default_folder(settings: TransformSettings, selection: CaseSelection) -> str:
    return f"{selection.case_id}{settings.label}"


def command(settings: TransformSettings, selection: CaseSelection, out_dir) -> list:
    """
    ``dfm case-transform CASE --rs RS ...`` (die CLI legt ``<case_id><label>``
    unter ``--output`` an; wie bei den Dosisindizes) bzw. ``dfm ct-transform
    CT_DIR ... --output <out_dir>/CT --save-viz <out_dir>/visualization_3d.html``.
    """
    out_dir = Path(out_dir)
    if selection.rs:
        case_dir = selection.case_dir or str(Path(selection.rs).parent)
        exact = out_dir.name == default_folder(settings, selection)
        return [["dfm", CASE_TOOL, str(case_dir), "--rs", str(selection.rs),
                 *settings.to_argv(CASE_TOOL), "--output", str(out_dir.parent if exact else out_dir)]]
    argv = ["dfm", CT_TOOL, str(selection.ct_dir or ""), *settings.to_argv(CT_TOOL),
            "--output", str(out_dir / "CT")]
    if settings.viz:
        argv += ["--save-viz", str(out_dir / VIZ_HTML)]
    return [argv]


def command_gaps(selection: CaseSelection) -> list:
    """Unterschiede zwischen der Auswahl und dem, was der CLI-Befehl verwenden wuerde."""
    def r(files):
        return [str(Path(f).resolve()) for f in files]

    want = r(selection.ct_files)
    if selection.rs:
        case_dir = Path(selection.case_dir or Path(selection.rs).parent)
        got = r(mod.ct_dir_files(case_dir / "CT")) if (case_dir / "CT").is_dir() else []
        where = "im Unterordner CT des Fallordners"
    else:
        got = r(mod.ct_dir_files(selection.ct_dir)) if selection.ct_dir else []
        where = "als alle .dcm-Dateien eines Ordners"
    return [] if got == want else [f"CT-Schichten (die CLI erwartet sie {where})"]


def planned_stages(settings: TransformSettings, selection: CaseSelection) -> list:
    """Stufen-Schluessel des Laufs in Reihenfolge (fuer ``stage i/n`` im Worker)."""
    if selection.rs:
        return (["preflight", "plan", "ct_transform", "rs_write"]
                + (["verify"] if settings.verify else []) + (["viz"] if settings.viz else []))
    return ["load", "transform"] + (["viz"] if settings.viz else [])


def job_checks(settings: TransformSettings, selection: CaseSelection) -> list:
    """Vor dem Lauf im Worker: reicht der freie Arbeitsspeicher (nur Header gelesen)?"""
    avail = available_memory_bytes()
    if avail is None or not selection.ct_files:
        return []
    try:
        shape = mod.extract_geometry(mod.load_ct_headers(selection.ct_files))["shape"]
    except Exception:  # noqa: BLE001 - der Lauf meldet den Fehler selbst
        return []
    need = estimate_memory_bytes(shape, settings.method, settings.order,
                                 settings.viz and settings.viz_ct_surface and bool(selection.rs))
    if need <= avail:
        return []
    return [Issue("warning", "SYSTEM.MEMORY_LOW",
                  f"Geschaetzter Speicherbedarf {format_bytes(need)}, frei {format_bytes(avail)}.",
                  hint_de="Andere Programme schliessen oder die Methode 'metadata' waehlen.",
                  field="method")]


def run(settings: TransformSettings, selection: CaseSelection, out_dir) -> JobResult:
    """Transformation nach ``out_dir`` (``CT/`` und ggf. ``RS<label>.dcm``, Ansichten)."""
    out_dir = Path(out_dir)
    cmd = command(settings, selection, out_dir)
    pre_issues = settings.validate(inspect_light(selection))
    if not selection.ct_files:
        pre_issues.append(Issue("error", "TRANSFORM.NO_CT", "Keine CT-Serie gewaehlt.",
                                field="ct_files"))
    if has_errors(pre_issues):
        return blocked(WORKFLOW, out_dir, pre_issues, cmd)
    gaps = command_gaps(selection)
    if gaps:
        pre_issues.append(Issue("info", "TRANSFORM.COMMAND_APPROX",
                                "Der CLI-Befehl gibt die Dateiauswahl nicht genau wieder: "
                                + "; ".join(gaps)))
    s = settings

    def body_case():
        pre = cm.preflight_case(selection.case_dir, rs_override=selection.rs, label=s.label,
                                ct_files=selection.ct_files, case_id=selection.case_id or None,
                                siblings=selection.related)
        center, center_label = cm.resolve_center(s.center, pre.rs_ds, pre.volume_center,
                                                 interactive=False)
        plan = cm.plan_transform(pre, s.tx, s.ty, s.tz, s.rx, s.ry, s.rz, out_dir=str(out_dir),
                                 method=s.method, order=s.order, label=s.label, center=center,
                                 center_label=center_label,
                                 new_frame_of_reference=s.new_frame_of_reference)
        res = cm.execute_transform(pre, plan, verify=s.verify, no_viz=not s.viz,
                                   viz_ct_surface=s.viz_ct_surface)
        outputs = {"ct_dir": res["ct_output_dir"], "rs": res["rs_output_path"]}
        viz = res.get("viz") or {"written": {}, "skipped": {}}
        outputs.update({f"viz:{k}": v for k, v in viz["written"].items()})
        issues = [Issue(**d) for d in res["issues"]]
        issues += [Issue("warning", "TRANSFORM.VIZ_SKIPPED", f"Ansicht {k} uebersprungen: {v}")
                   for k, v in viz["skipped"].items()]
        summary = {"case_id": res["case_id"], "tool": CASE_TOOL, "method": res["method"],
                   "for_strategy": res["for_strategy"], "n_slices": len(res["sop_map"]),
                   "rotation_center_mm": res["rotation_center"],
                   "rotation_center_label": res["rotation_center_label"],
                   "drehpunkt_mm": res["drehpunkt_pos"], "n_clipped_rois": len(res["clipping"]),
                   "clipping": res["clipping"], "verify": res["verify"], "T": plan.T.tolist(),
                   "description": describe_motion(s)}
        return outputs, issues, summary

    def body_ct():
        res = mod.run_ct_transform(list(selection.ct_files), out_dir / "CT", s.tx, s.ty, s.tz,
                                   s.rx, s.ry, s.rz, method=s.method, order=s.order, viz=s.viz,
                                   viz_html=str(out_dir / VIZ_HTML))
        outputs = {"ct_dir": res["output_dir"], "viz_html": res["viz_html_path"]}
        summary = {"case_id": selection.case_id, "tool": CT_TOOL, "method": res["method"],
                   "n_slices": res["n_slices"], "rotation_center_mm": res["rotation_center"],
                   "rotation_center_label": "Volumenmitte", "T": res["T"], "description": describe_motion(s)}
        return outputs, [], summary

    return run_guarded(WORKFLOW, out_dir, body_case if selection.rs else body_ct, pre_issues, cmd)


class _LightInfo:
    """Minimal-``info`` fuer die Einstellungspruefung ohne Header-Scan."""

    def __init__(self, has_rs: bool):
        self.has_rs = has_rs


def inspect_light(selection: CaseSelection) -> _LightInfo:
    return _LightInfo(bool(selection.rs))
