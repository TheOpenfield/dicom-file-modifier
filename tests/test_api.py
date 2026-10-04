"""
api/ (Plan P0.8b): Settings <-> argparse, API = CLI, stille Inspektion,
Vorschau = Lauf, Ergebnisse und Fehlerklassen.  Alles auf dem Demo-Fall.
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import importlib
import json
import shutil
import subprocess
import sys
from pathlib import Path
from unittest import mock

import pydicom
import pytest

from dicom_file_modifier import _runtime
from dicom_file_modifier import case_modifier as cm
from dicom_file_modifier import dose_indices as di
from dicom_file_modifier.api import dose, selection, structures, transform
from dicom_file_modifier.api.fields import command_string
from dicom_file_modifier.api.issues import exit_code_for, issue_from_exception
from dicom_file_modifier.api.outputs import OutputSpec
from dicom_file_modifier.api.results import JobResult
from dicom_file_modifier.demo import DemoSpec, make_demo_case
from dicom_file_modifier.issues import Issue, UserInputError


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def _hashes(root) -> dict:
    return {p: _sha(p) for p in Path(root).rglob("*") if p.is_file()}


@pytest.fixture(scope="module")
def demo(tmp_path_factory):
    case = make_demo_case(tmp_path_factory.mktemp("demo") / "case", DemoSpec())
    case.hashes = _hashes(case.root)          # letzter Test: Eingaben unveraendert
    return case


# -- Settings <-> argparse ----------------------------------------------------

class _Captured(Exception):
    pass


def _parser(module: str) -> argparse.ArgumentParser:
    """Den Parser abfangen, den ``main()`` des Moduls baut."""
    def fake(self, args=None, namespace=None):
        raise _Captured(self)

    mod = importlib.import_module(f"dicom_file_modifier.{module}")
    with mock.patch.object(argparse.ArgumentParser, "parse_args", fake):
        try:
            mod.main([])
        except _Captured as c:
            return c.args[0]
    raise AssertionError(f"{module}.main() baut keinen Parser")


PARITY = [  # (Settings, dfm-Befehl, Modul, Positionsargumente)
    (structures.StructureSettings, "analyze", "analyzer", ["RS.dcm"]),
    (structures.StructureSettings, "visualize", "visualizer", ["RS.dcm"]),
    (dose.DoseIndexSettings, "dose-indices", "dose_indices", ["CASE"]),
    (transform.TransformSettings, "case-transform", "case_modifier", ["CASE"]),
    (transform.TransformSettings, "ct-transform", "modifier", ["CT"]),
]
NON_SETTINGS = {**structures.NON_SETTINGS, **dose.NON_SETTINGS, **transform.NON_SETTINGS}
SAMPLES = {   # je Klasse Werte abseits der Defaults (alle Felder mit CLI-Flag)
    structures.StructureSettings: [dict(targets=["PTV", "GTV"], oars=["Hirnstamm", "Chiasma"])],
    dose.DoseIndexSettings: [
        dict(target=["PTV_1", "PTV_2"], rx=18.5, isodose="100,80,50,12Gy", eclipse_compat="high",
             write_rs=False, viz=False, label="_X", grid_mm=0.5, dose_interp="cubic",
             volume_model="eclipse", piv_scope="global", iso_contours="field",
             eclipse_values="TV=1.2,PIV=1.45", eclipse_tol_pct=3.0, eclipse_dvh=False,
             append_csv="C:/sammel/indices.csv", include_target=True, viz_ct=False,
             simplify_mm=0.05, transfer_syntax="implicit", max_name_len=32),
        dict(rx_pct_of_max=80.0, eclipse_compat="default"),
    ],
    transform.TransformSettings: [
        dict(tx=-10.5, ty=2.0, tz=-0.25, rx=1.5, ry=-2.0, rz=15.0, center="-10,5,3",
             method="metadata", label="_T1", new_frame_of_reference=True, order=3, verify=True,
             viz=False, viz_ct_surface=True),
        dict(center="marker:HS1", order=0),
    ],
}


def _mapped(cls, tool):
    return {name: meta for name, meta in cls.field_meta().items()
            if meta.cli_flag and meta.applies_to(tool)}


@pytest.mark.parametrize("cls, tool, module, pos", PARITY, ids=[p[1] for p in PARITY])
def test_every_cli_option_is_a_field_or_known_non_setting(cls, tool, module, pos):
    actions = {a.dest: a for a in _parser(module)._actions}
    mapped = {meta.dest: name for name, meta in _mapped(cls, tool).items()}
    assert set(actions) - set(mapped) == NON_SETTINGS[tool]
    assert set(mapped) <= set(actions)


@pytest.mark.parametrize("cls, tool, module, pos", PARITY, ids=[p[1] for p in PARITY])
def test_defaults_and_choices_match_argparse(cls, tool, module, pos):
    actions = {a.dest: a for a in _parser(module)._actions}
    s = cls()
    for name, meta in _mapped(cls, tool).items():
        a = actions[meta.dest]
        assert a.default == s._cli_default(name, meta), name
        if not isinstance(a.default, bool):
            assert tuple(a.choices or ()) == meta.choices, name


@pytest.mark.parametrize("cls, tool, module, pos", PARITY, ids=[p[1] for p in PARITY])
def test_settings_round_trip_through_argparse(cls, tool, module, pos):
    parser = _parser(module)
    for values in SAMPLES[cls]:
        s = cls(**values)
        ns = parser.parse_args(pos + s.to_argv(tool))
        applicable = {k: v for k, v in values.items() if k in _mapped(cls, tool)}
        assert cls.from_namespace(ns, tool) == cls(**applicable), values
    assert cls().to_argv(tool) == ([] if cls is not transform.TransformSettings or tool == "ct-transform"
                                   else ["--center", "volume"])


def test_dose_cli_hands_the_core_the_same_arguments_as_the_api(monkeypatch, tmp_path):
    for values in SAMPLES[dose.DoseIndexSettings]:
        captured = {}
        monkeypatch.setattr(di, "run_dose_indices", lambda *a, **kw: captured.update(kw))
        s = dose.DoseIndexSettings(**values)
        assert di.main([str(tmp_path)] + s.to_argv("dose-indices") + ["--output", "o"]) == 0
        for key in ("output", "rs", "rd", "rp", "eclipse_ref"):
            captured.pop(key)
        assert captured == dose.run_kwargs(s)


def test_settings_dict_round_trip_and_type_rules():
    s = dose.DoseIndexSettings(**SAMPLES[dose.DoseIndexSettings][0])
    assert dose.DoseIndexSettings.from_dict(json.loads(json.dumps(s.to_dict()))) == s
    t = transform.TransformSettings.from_dict({"tx": 3, "order": 3.0})
    assert t.tx == 3.0 and isinstance(t.tx, float) and t.order == 3
    with pytest.raises(ValueError, match="Unbekannte"):
        transform.TransformSettings.from_dict({"tx": 1, "gibtsnicht": 2})
    assert transform.TransformSettings.from_dict({"gibtsnicht": 2}, strict=False) == transform.TransformSettings()
    with pytest.raises(ValueError, match="verify"):
        transform.TransformSettings.from_dict({"verify": "ja"})
    assert transform.TransformSettings(tx=1.0).non_default() == {"tx": 1.0}


def test_validation_ranges_choices_and_locks():
    codes = {i.code: i.field for i in dose.DoseIndexSettings(
        grid_mm=0.3, rx=20.0, rx_pct_of_max=80.0, max_name_len=100, label="_a:b").validate()}
    assert codes == {"SETTINGS.CHOICE": "grid_mm", "SETTINGS.RANGE": "max_name_len",
                     "DOSE.RX_CONFLICT": "rx", "DOSE.LABEL_INVALID": "label"}
    locked = dose.DoseIndexSettings(eclipse_compat="high", write_rs=False).disabled_fields()
    assert set(dose.ECLIPSE_LOCKED) | {"include_target", "transfer_syntax", "max_name_len"} == set(locked)
    assert transform.TransformSettings(method="metadata").disabled_fields() == {"order": mock.ANY}


def test_negative_values_are_passed_with_equals_sign():
    argv = transform.TransformSettings(tx=-5.0, center="-1,2,3").to_argv("case-transform")
    assert "--tx=-5.0" in argv and "--center=-1,2,3" in argv


# -- Inspektion und Vorschau ----------------------------------------------------

def test_inspect_and_preview_print_nothing(demo, capfd):
    s_info = structures.inspect(selection.from_rtstruct(str(demo.rs)))
    structures.preview(s_info, structures.Settings())
    d_info = dose.inspect(selection.for_dose(str(demo.root)))
    dose.preview(d_info, dose.Settings())
    dose.preview(d_info, dose.Settings(eclipse_compat="default"))
    t_info = transform.inspect(selection.for_transform(str(demo.root)))
    transform.preview(t_info, transform.Settings(rz=5.0, center="marker:HS1"))
    transform.preview(transform.inspect(selection.from_ct_dir(str(demo.ct_dir))), transform.Settings(tx=1.0))
    out, err = capfd.readouterr()
    assert out == "" and err == ""
    assert s_info.ok and d_info.ok and t_info.ok


def test_api_import_and_inspection_load_no_plot_modules(demo):
    code = ("import sys\n"
            "from dicom_file_modifier.api import selection as S, structures, dose, transform, results\n"
            f"case = {str(demo.root)!r}\n"
            "structures.preview(structures.inspect(S.from_rtstruct(" + repr(str(demo.rs)) + ")), structures.Settings())\n"
            "dose.preview(dose.inspect(S.for_dose(case)), dose.Settings())\n"
            "transform.preview(transform.inspect(S.for_transform(case)), transform.Settings(tx=1.0))\n"
            "print(sorted(m for m in ('matplotlib.pyplot', 'dicom_file_modifier.visualizer',\n"
            "                         'dicom_file_modifier.dose_viz', 'plotly') if m in sys.modules))")
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=300)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == "[]"


def test_dose_inspection_reports_candidates(demo):
    info = dose.inspect(selection.for_dose(str(demo.root)))
    assert info.targets["default"] == ["PTV_1"] and set(info.targets["ptvs"]) == {"PTV_1", "PTV_2"}
    assert info.prescriptions["default"]["target_prescription_dose_gy"] == 20.0
    assert info.dvh["available"] and info.dvh["trusted"] and info.dvh["body"] == "BODY"
    assert info.rs_export["possible"] and info.ct["n_slices"] == 80
    assert info.eclipse_compat["high"]["grid_mm"] == 2.0
    assert info.issues == []


def test_preview_problems_point_at_the_field(demo):
    info = dose.inspect(selection.for_dose(str(demo.root)))

    def error_fields(**kw):
        return [i.field for i in dose.preview(info, dose.Settings(**kw)).issues if i.level == "error"]

    assert error_fields(target=["PTV"]) == ["target"]          # PTV_1, PTV_2, h_PTV_gesamt
    assert error_fields(isodose="100,abc") == ["isodose"]
    assert error_fields(eclipse_values="TV=x") == ["eclipse_values"]
    pv = dose.preview(info, dose.Settings(grid_mm=0.1))        # passt beim kleinen Demo-Fall
    assert pv.ok and 0 < pv.grid["n_voxels"] <= pv.grid["limit"]


def test_too_large_fine_grid_is_reported_with_its_size(demo, monkeypatch):
    info = dose.inspect(selection.for_dose(str(demo.root)))
    monkeypatch.setattr(di.dm, "MAX_FINE_VOXELS", 1000)
    pv = dose.preview(info, dose.Settings())
    assert not pv.ok and pv.grid["n_voxels"] > 1000
    assert [(i.code, i.field) for i in pv.issues if i.level == "error"] == [("DOSE.GRID_TOO_LARGE", "grid_mm")]


def test_eclipse_mode_preview_matches_resolve_effective(demo):
    info = dose.inspect(selection.for_dose(str(demo.root)))
    for mode in ("high", "default"):
        s = dose.Settings(eclipse_compat=mode)
        values, locked, issues = dose.resolve_effective(s, info)
        pv = dose.preview(info, s)
        assert pv.ok and not issues and set(locked) == set(dose.ECLIPSE_LOCKED)
        assert {k: pv.effective[k] for k in values} == values


def test_eclipse_mode_without_ct_is_an_error(demo):
    sel = dataclasses.replace(selection.for_dose(str(demo.root)), ct_files=[])
    pv = dose.preview(dose.inspect(sel), dose.Settings(eclipse_compat="high"))
    assert [i.code for i in pv.issues if i.level == "error"] == ["DOSE.ECLIPSE_COMPAT_UNAVAILABLE"]


def test_transform_preview(demo):
    info = transform.inspect(selection.for_transform(str(demo.root)))
    pv = transform.preview(info, transform.Settings(tx=10.0, tz=-5.0, rz=15.0, center="marker:HS1"))
    assert pv.ok and pv.center_label == "Marker 'HS1'" and pv.drehpunkt_mm == [20.0, 5.0, -20.5]
    assert pv.description == "10 mm nach links · 5 mm nach inferior · +15° um die Kopf-Fuss-Achse"
    assert pv.planned == {"ct_dir": "CT", "rs": "RS_RB.dcm"} and pv.memory_bytes > 0
    bad = transform.preview(info, transform.Settings(tx=1.0, center="marker:Gibtsnicht"))
    assert [(i.field) for i in bad.issues if i.level == "error"] == ["center"]
    ct_only = transform.inspect(selection.from_ct_dir(str(demo.ct_dir)))
    codes = [i.code for i in transform.preview(ct_only, transform.Settings(center="1,2,3")).issues]
    assert "TRANSFORM.CENTER_NEEDS_RS" in codes


def test_transform_inspection_checks_the_ct_references(demo, tmp_path):
    case = tmp_path / "case"
    shutil.copytree(demo.root, case)
    rs_file = next(case.glob("RS*.dcm"))
    ds = pydicom.dcmread(str(rs_file))
    c = next(c for rc in ds.ROIContourSequence for c in rc.get("ContourSequence", [])
             if "ContourImageSequence" in c)
    c.ContourImageSequence[0].ReferencedSOPInstanceUID = "1.2.826.0.1.3680043.8.498.97"
    ds.save_as(str(rs_file))
    info = transform.inspect(selection.for_transform(str(case)))
    errors = [i for i in info.issues if i.level == "error"]
    assert [i.code for i in errors] == ["INPUT.MISSING_KEY"] and "unbekannte CT-SOP" in errors[0].message_de
    assert info.ct and not transform.preview(info, transform.Settings(tx=1.0)).ok


# -- API = CLI -------------------------------------------------------------------

def test_structures_api_equals_cli(demo, tmp_path, capsys):
    res = structures.run(structures.Settings(), selection.from_rtstruct(str(demo.rs)), tmp_path / "api")
    assert res.status == "ok" and res.exit_code == 0
    cli = tmp_path / "cli"
    from dicom_file_modifier import analyzer, visualizer
    assert analyzer.main([str(demo.rs), "--output", str(cli)]) == 0
    assert visualizer.main([str(demo.rs), "--output", str(cli)]) == 0
    capsys.readouterr()
    api_files = {p.name: _sha(p) for p in (tmp_path / "api").iterdir()}
    cli_files = {p.name: _sha(p) for p in cli.iterdir()}
    assert api_files == cli_files
    assert res.command[0][:2] == ["dfm", "analyze"] and res.command[1][:2] == ["dfm", "visualize"]


def _report(path: Path) -> dict:
    r = json.loads(path.read_text(encoding="utf-8"))
    r["meta"].pop("timestamp")
    r.pop("outputs")
    return r


def _rois(rs_path: Path) -> dict:
    ds = pydicom.dcmread(str(rs_path))
    names = {int(r.ROINumber): str(r.ROIName) for r in ds.StructureSetROISequence}
    return {names[int(rc.ReferencedROINumber)]: [list(c.ContourData) for c in rc.get("ContourSequence", [])]
            for rc in ds.ROIContourSequence}


def test_dose_api_equals_cli_and_preview_names_the_rois(demo, tmp_path, capsys):
    sel = selection.for_dose(str(demo.root))
    s = dose.Settings(viz=False, include_target=True)
    info = dose.inspect(sel)
    pv = dose.preview(info, s)
    res = dose.run(s, sel, tmp_path / "api" / dose.default_folder(s, sel))
    cid = sel.case_id
    assert res.ok and res.outputs["json"] == f"{cid}_indices.json" and res.outputs["rs"] == f"RS_{cid}_IDX.dcm"
    assert di.main([str(demo.root), "--no-viz", "--include-target", "--output", str(tmp_path / "cli")]) == 0
    capsys.readouterr()
    api_dir, cli_dir = tmp_path / "api" / f"{cid}_IDX", tmp_path / "cli" / f"{cid}_IDX"
    assert _report(api_dir / f"{cid}_indices.json") == _report(cli_dir / f"{cid}_indices.json")
    api_rois = _rois(api_dir / f"RS_{cid}_IDX.dcm")
    assert api_rois == _rois(cli_dir / f"RS_{cid}_IDX.dcm")
    assert [name for _, name, _ in pv.roi_names] == list(api_rois)
    assert res.command == [["dfm", "dose-indices", str(demo.root), "--no-viz", "--include-target",
                            "--output", str(tmp_path / "api")]]


def _ct_geometry(ct_dir: Path) -> list:
    out = []
    for p in sorted(ct_dir.glob("*.dcm")):
        ds = pydicom.dcmread(str(p), stop_before_pixels=True)
        out.append((p.name, list(ds.ImagePositionPatient), list(ds.ImageOrientationPatient),
                    int(ds.SeriesNumber)))
    return out


def test_case_transform_api_equals_cli(demo, tmp_path, capsys):
    s = transform.Settings(tx=3.0, ty=-2.0, rz=5.0, method="metadata", viz=False,
                           center="marker:HS1", verify=True)
    sel = selection.for_transform(str(demo.root))
    api_dir = tmp_path / "api" / f"{sel.case_id}_RB"
    res = transform.run(s, sel, api_dir)
    assert res.ok and res.outputs == {"ct_dir": "CT", "rs": "RS_RB.dcm"}
    assert res.summary["verify"]["passed"] and res.summary["rotation_center_label"] == "Marker 'HS1'"
    assert {i.code for i in res.issues} == {"CASE.FOR_KEPT", "CASE.SIBLINGS_NOT_TRANSFORMED"}
    assert res.command == [["dfm", "case-transform", str(demo.root), "--rs", str(demo.rs),
                            "--tx", "3.0", "--ty=-2.0", "--rz", "5.0", "--center", "marker:HS1",
                            "--method", "metadata", "--verify", "--no-viz",
                            "--output", str(tmp_path / "api")]]
    argv = res.command[0][2:-1] + [str(tmp_path / "cli")]      # kopierter Befehl, anderer Ordner
    assert cm.main(argv) == 0
    capsys.readouterr()
    cli_dir = tmp_path / "cli" / f"{sel.case_id}_RB"
    assert _ct_geometry(cli_dir / "CT") == _ct_geometry(api_dir / "CT")
    assert _rois(cli_dir / "RS_RB.dcm") == _rois(api_dir / "RS_RB.dcm")
    assert "Drehpunkt" in _rois(api_dir / "RS_RB.dcm")


def test_ct_only_transform_api_equals_cli(demo, tmp_path, capsys):
    s = transform.Settings(tx=-4.0, rx=2.0, method="metadata", viz=False)
    res = transform.run(s, selection.from_ct_dir(str(demo.ct_dir)), tmp_path / "api")
    assert res.ok and res.outputs == {"ct_dir": "CT"} and len(res.manifest) == 80
    assert res.command == [["dfm", "ct-transform", str(demo.ct_dir), "--tx=-4.0", "--rx", "2.0",
                            "--method", "metadata", "--no-viz", "--output",
                            str(tmp_path / "api" / "CT")]]
    from dicom_file_modifier import modifier
    argv = res.command[0][2:-1] + [str(tmp_path / "cli")]
    assert modifier.main(argv) == 0
    capsys.readouterr()
    geo = _ct_geometry(tmp_path / "api" / "CT")
    assert geo == _ct_geometry(tmp_path / "cli") and geo != _ct_geometry(demo.ct_dir)


def test_flat_export_selection_writes_the_isodose_rtstruct(tmp_path, capsys):
    flat = make_demo_case(tmp_path / "flat", DemoSpec(layout="flat"))
    sel = selection.CaseSelection(case_id="flat", ct_files=[str(p) for p in flat.ct_files],
                                  rs=str(flat.rs), rd=str(flat.rd[0]), rp=str(flat.rp))
    info = dose.inspect(sel)
    assert info.ok and info.rs_export["possible"]
    res = dose.run(dose.Settings(viz=False), sel, tmp_path / "out" / "flat_IDX")
    capsys.readouterr()
    assert res.ok and res.outputs["rs"] == "RS_flat_IDX.dcm"
    assert "DOSE.COMMAND_APPROX" in {i.code for i in res.issues}


# -- Ergebnisse und Fehlerklassen -------------------------------------------------

def test_job_result_round_trip(demo, tmp_path, capsys):
    res = transform.run(transform.Settings(tx=1.0, method="metadata", viz=False),
                        selection.from_ct_dir(str(demo.ct_dir)), tmp_path / "o")
    capsys.readouterr()
    back = JobResult.from_dict(json.loads(json.dumps(res.to_dict())))
    assert back == res and [s["key"] for s in res.timings["stages"]] == ["load", "transform"]


def test_input_error_internal_error_and_cancel(demo, tmp_path, monkeypatch, capsys):
    sel = selection.for_dose(str(demo.root))
    missing = dataclasses.replace(sel, rd=str(tmp_path / "fehlt.dcm"))
    res = dose.run(dose.Settings(viz=False), missing, tmp_path / "a")
    assert (res.status, res.exit_code, res.issues[-1].code) == ("failed", 2, "FILE.NOT_FOUND")

    monkeypatch.setattr(di, "compute_dose_indices", mock.Mock(side_effect=RuntimeError("kaputt")))
    res = dose.run(dose.Settings(viz=False), sel, tmp_path / "b")
    assert (res.status, res.exit_code, res.issues[-1].code) == ("failed", 1, "INTERNAL")
    assert "RuntimeError: kaputt" in res.issues[-1].detail
    monkeypatch.undo()

    class _Cancel(_runtime.JobContext):
        def check_cancel(self):
            raise _runtime.JobCancelled()

    with _runtime.use(_Cancel()):
        res = transform.run(transform.Settings(tx=1.0, method="resample", viz=False),
                            selection.from_ct_dir(str(demo.ct_dir)), tmp_path / "c")
    capsys.readouterr()
    assert (res.status, res.exit_code, res.issues[-1].code) == ("cancelled", 3, "JOB.CANCELLED")


def test_blocked_run_writes_nothing(demo, tmp_path):
    res = dose.run(dose.Settings(rx=20.0, rx_pct_of_max=80.0), selection.for_dose(str(demo.root)),
                   tmp_path / "x")
    assert (res.status, res.exit_code) == ("failed", 2) and not (tmp_path / "x").exists()


@pytest.mark.parametrize("exc, code, exit_code", [
    (FileNotFoundError(2, "No such file", "a.dcm"), "FILE.NOT_FOUND", 2),
    (PermissionError(13, "Permission denied", "indices.csv"), "FILE.PERMISSION", 2),
    (MemoryError(), "SYSTEM.MEMORY", 1),
    (AttributeError("'FileDataset' object has no attribute 'PixelSpacing'"), "DICOM.MISSING_TAG", 2),
    (AttributeError("'NoneType' object has no attribute 'x'"), "INTERNAL", 1),
    (ValueError("schlecht"), "INPUT.INVALID", 2),
    (UserInputError(Issue("error", "CT.X", "Text")), "CT.X", 2),
    (_runtime.JobCancelled(), "JOB.CANCELLED", 3),
])
def test_exceptions_become_plain_german_issues(exc, code, exit_code):
    assert issue_from_exception(exc).code == code and exit_code_for(exc) == exit_code


def test_invalid_dicom_is_an_input_error(tmp_path):
    bad = tmp_path / "kein.dcm"
    bad.write_bytes(b"kein DICOM")
    with pytest.raises(Exception) as e:
        pydicom.dcmread(str(bad))
    assert issue_from_exception(e.value).code == "DICOM.INVALID" and exit_code_for(e.value) == 2


def test_output_spec(tmp_path):
    (tmp_path / "demo_IDX").mkdir()
    (tmp_path / "demo_IDX_2").mkdir()
    assert OutputSpec(str(tmp_path)).target_dir("demo_IDX") == tmp_path / "demo_IDX_3"
    assert OutputSpec(str(tmp_path), policy="overwrite").target_dir("demo_IDX") == tmp_path / "demo_IDX"
    codes = {i.code for i in OutputSpec("relativ", folder="CON", policy="x").validate()}
    assert codes == {"OUTPUT.ROOT_RELATIVE", "OUTPUT.POLICY", "OUTPUT.FOLDER_INVALID"}


def test_selection_matches_cli_discovery(demo):
    sel = selection.for_dose(str(demo.root))
    files = di.discover_dose_case(str(demo.root))
    assert (Path(sel.rs), Path(sel.rd), Path(sel.rp)) == (files["rs"], files["rd"], files["rp"])
    assert [Path(f) for f in sel.ct_files] == sorted(files["ct_dir"].glob("*.dcm"))
    assert selection.CaseSelection.from_dict(sel.to_dict()) == sel
    t = selection.for_transform(str(demo.root))
    assert t.rs == sel.rs and t.ct_files == sel.ct_files and len(t.related) == 2


def test_command_string_quotes_for_windows():
    text = command_string([["dfm", "analyze", "C:/Mit Leerzeichen/RS.dcm"]])
    assert text == 'dfm analyze "C:/Mit Leerzeichen/RS.dcm"'


def test_runs_never_modify_the_inputs(demo):
    """Zuletzt: alle Laeufe dieser Datei haben den Fallordner nur gelesen."""
    assert _hashes(demo.root) == demo.hashes
