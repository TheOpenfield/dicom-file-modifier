"""
Workflow Strukturanalyse: ``analyzer`` (JSON) + ``visualizer`` (Plots,
statistics.txt) auf einem RTSTRUCT.

CLI-Entsprechung: ``dfm analyze RS --output DIR`` und ``dfm visualize RS
--output DIR`` mit denselben ``--targets``/``--oars``; ``run`` analysiert nur
einmal (die Zahlen sind je Lauf deterministisch) und schreibt beides.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from .. import _runtime
from .. import analyzer as ana
from .fields import SettingsBase, setting
from .issues import Issue, has_errors, issue_from_exception
from .results import JobResult, blocked, run_guarded
from .selection import CaseSelection

WORKFLOW = "structures"
TOOLS = ("analyze", "visualize")
# CLI-Optionen ohne Settings-Feld: Auswahl, Ausgabe, Betriebsarten
NON_SETTINGS = {"analyze": {"file", "list", "output", "self_test", "help"},
                "visualize": {"file", "output", "help"}}


@dataclass
class StructureSettings(SettingsBase):
    WORKFLOW = WORKFLOW

    targets: Optional[list[str]] = setting(
        None, label="Zielvolumen", kind="names", cli_flag="--targets",
        help="Namen oder Namensteile (Gross-/Kleinschreibung egal); leer = automatisch nach "
             "DICOM-Typ und Name.")
    oars: Optional[list[str]] = setting(
        None, label="Risikoorgane", kind="names", cli_flag="--oars",
        help="Namen oder Namensteile; leer = automatisch (serielle und parallele Organe).")
    write_json: bool = setting(
        True, label="Analyse als JSON", tools=("analyze",),
        help="<RS>_analysis.json mit allen Kennzahlen schreiben.")
    plots: bool = setting(
        True, label="Plots und statistics.txt", tools=("visualize",),
        help="Volumen-, Form-, Abstands- und SRS-Plots sowie die Zahlenzusammenfassung.")

    def _check(self, info) -> list:
        if not (self.write_json or self.plots):
            return [Issue("error", "STRUCT.NOTHING_TO_DO", "Weder JSON noch Plots gewaehlt.",
                          field="plots")]
        return []


Settings = StructureSettings


@dataclass
class StructureCaseInfo:
    selection: CaseSelection
    rtstruct: dict = field(default_factory=dict)     # analyzer.inspect_rtstruct
    auto_targets: list = field(default_factory=list)
    auto_oars: list = field(default_factory=list)
    issues: list = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not has_errors(self.issues)


def inspect(selection: CaseSelection) -> StructureCaseInfo:
    """ROI-Tabelle mit Kategorien, Marker, FoRs und die automatische Auswahl; still."""
    info = StructureCaseInfo(selection=selection)
    if not selection.rs:
        info.issues.append(Issue("error", "STRUCT.NO_RS", "Kein RTSTRUCT gewaehlt.", field="rs"))
        return info
    try:
        info.rtstruct = ana.inspect_rtstruct(selection.rs, volumes=True)
    except Exception as exc:  # noqa: BLE001 - Befund statt Absturz
        info.issues.append(issue_from_exception(exc, field="rs"))
        return info
    rois = info.rtstruct["rois"]
    info.auto_targets = [r["name"] for r in rois if r["category"] == ana.CAT_TARGET]
    info.auto_oars = [r["name"] for r in rois
                      if r["category"] in (ana.CAT_OAR_SERIAL, ana.CAT_OAR_PARALLEL)]
    if not info.auto_targets:
        info.issues.append(Issue("warning", "STRUCT.NO_TARGETS", "Keine Zielvolumen erkannt.",
                                 hint_de="Zielvolumen per Name waehlen.", field="targets"))
    return info


def preview(info: StructureCaseInfo, settings: StructureSettings) -> dict:
    """Welche ROIs als Ziel, Risikoorgan und Hilfsstruktur ausgewertet wuerden."""
    issues = list(info.issues) + settings.validate(info)
    if not info.rtstruct:
        return {"targets": [], "oars": [], "helpers": [], "issues": issues, "ok": False}
    rois = info.rtstruct["rois"]
    names = {r["number"]: r["name"] for r in rois}
    cats = {r["number"]: r["category"] for r in rois}
    t, o, h = ana.select_rois(names, cats, settings.targets, settings.oars)
    pick = [names[n] for n in sorted(t)], [names[n] for n in sorted(o)], [names[n] for n in sorted(h)]
    if settings.targets and not t:
        issues.append(Issue("warning", "STRUCT.TARGETS_UNMATCHED",
                            "Kein ROI-Name passt zu den Zielvolumen-Angaben.", field="targets"))
    if settings.oars and not o:
        issues.append(Issue("warning", "STRUCT.OARS_UNMATCHED",
                            "Kein ROI-Name passt zu den Risikoorgan-Angaben.", field="oars"))
    return {"targets": pick[0], "oars": pick[1], "helpers": pick[2], "issues": issues,
            "ok": not has_errors(issues)}


def default_folder(settings: StructureSettings, selection: CaseSelection) -> str:
    return f"{selection.case_id or Path(selection.rs).parent.name}_STRUCT"


def command(settings: StructureSettings, selection: CaseSelection, out_dir) -> list:
    cmds = []
    if settings.write_json:
        cmds.append(["dfm", "analyze", str(selection.rs), *settings.to_argv("analyze"),
                     "--output", str(out_dir)])
    if settings.plots:
        cmds.append(["dfm", "visualize", str(selection.rs), *settings.to_argv("visualize"),
                     "--output", str(out_dir)])
    return cmds


def planned_stages(settings: StructureSettings, selection: CaseSelection) -> list:
    """Stufen-Schluessel des Laufs in Reihenfolge (fuer ``stage i/n`` im Worker)."""
    return ["analysis"] + (["plots"] if settings.plots else [])


def run(settings: StructureSettings, selection: CaseSelection, out_dir) -> JobResult:
    """Analyse (und Plots) nach ``out_dir``; Konsolentext wie die CLI (im Worker: Protokoll)."""
    out_dir = Path(out_dir)
    cmd = command(settings, selection, out_dir)
    pre = settings.validate() + ([] if selection.rs else
                                 [Issue("error", "STRUCT.NO_RS", "Kein RTSTRUCT gewaehlt.", field="rs")])
    if has_errors(pre):
        return blocked(WORKFLOW, out_dir, pre, cmd)

    def body():
        ctx = _runtime.current()
        ctx.stage("analysis", "Strukturen analysieren")
        results, info = ana.run_analysis(str(selection.rs), settings.targets, settings.oars,
                                         return_info=True)
        outputs, issues = {}, list(info.get("issues", []))
        if settings.write_json:
            outputs["analysis_json"] = str(ana.write_analysis_json(
                results, out_dir, Path(selection.rs).stem))
        if settings.plots:
            ctx.check_cancel()
            ctx.stage("plots", "Plots erstellen")
            from .. import visualizer as viz          # erst hier: laedt matplotlib
            report = viz.run_visualization(results, out_dir)
            for name, path in report["written"].items():
                outputs["statistics" if name == "statistics" else f"plot:{name}"] = path
            for name, reason in report["skipped"].items():
                no_data = reason == "keine passenden Daten"
                issues.append(Issue("info" if no_data else "warning", "STRUCT.PLOT_SKIPPED",
                                    f"Plot {name} uebersprungen: {reason}"))
        summary = {
            "n_targets": len(results["targets"]), "n_oars": len(results["oars"]),
            "n_helpers": len(results["helpers"]), "n_distances": len(results["distances"]),
            "closest": [{"a": d["structure_a"], "b": d["structure_b"],
                         "min_distance_mm": d["min_distance_mm"]}
                        for d in results["distances"][:5]],
        }
        return outputs, issues, summary

    return run_guarded(WORKFLOW, out_dir, body, pre, cmd)
