"""
dose_indices (Plan P0.4): Kandidaten statt Abbruch, Fixes, Stufen, Fortschritt
und Abbruch.  Alle Laeufe auf dem synthetischen Demo-Fall.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pydicom
import pytest

from dicom_file_modifier import _runtime
from dicom_file_modifier import dose_indices as di
from dicom_file_modifier import rtstruct_writer as rw
from dicom_file_modifier.demo import DemoSpec, make_demo_case
from dicom_file_modifier.dicom_utils import set_sop_instance_uid


@pytest.fixture(scope="module")
def demo(tmp_path_factory):
    return make_demo_case(tmp_path_factory.mktemp("demo") / "case", DemoSpec())


@pytest.fixture(scope="module")
def plan_inputs(demo):
    rs_ds = pydicom.dcmread(str(demo.rs))
    rp_refs = di.dm.prescription_references(pydicom.dcmread(str(demo.rp)))
    dose = di.dm.dose_grid_from_dataset(di.dm.load_rtdose(str(demo.rd[0])), str(demo.rd[0]))
    return rs_ds, rp_refs, dose


def _copy_case(demo, tmp_path: Path) -> Path:
    dst = tmp_path / "case"
    shutil.copytree(demo.root, dst)
    return dst


def _run(case_dir, tmp_path: Path, **kw):
    kw.setdefault("no_viz", True)
    kw.setdefault("quiet", True)
    return di.run_dose_indices_ex(str(case_dir), output=str(tmp_path / "out"), **kw)


class _Recorder(_runtime.JobContext):
    def __init__(self, cancel_after=None):
        self.stages, self.n_progress, self.n_checks = [], 0, 0
        self.cancel_after = cancel_after

    def stage(self, key, label):
        self.stages.append(key)

    def progress(self, done, total, text=""):
        self.n_progress += 1

    def check_cancel(self):
        self.n_checks += 1
        if self.cancel_after is not None and self.n_checks > self.cancel_after:
            raise _runtime.JobCancelled()


# -- Isodosen-Level ----------------------------------------------------------

def test_level_on_rx_given_in_gy_gets_the_canonical_key():
    levels, _ = di.parse_isodose_levels("20Gy", 20.0)
    assert [lv["key"] for lv in levels] == ["100", "50"]
    assert levels[0]["label"] == "20 Gy"
    levels, notes = di.parse_isodose_levels("100,10Gy", 20.0)
    assert [lv["key"] for lv in levels] == ["100", "50"] and not notes


def test_run_with_isodose_equal_to_rx(demo, tmp_path):
    report, art = _run(demo.root, tmp_path, isodose="20Gy")
    assert set(art.levels) == {"100", "50"}
    assert report["targets"]["PTV_1"]["indices"]["ci_paddick"] > 0


# -- Kandidaten statt Abbruch ------------------------------------------------

def test_target_candidates_match_select_targets(plan_inputs):
    rs_ds, rp_refs, _ = plan_inputs
    cand = di.target_candidates(rs_ds, rp_refs)
    assert [r[1] for r in cand["ptvs"]] == ["PTV_1", "PTV_2"]
    assert [n for _, n in cand["default"]] == ["PTV_1"]
    assert "RTPLAN" in cand["reason"]
    assert (cand["default"], cand["notes"]) == di.select_targets(rs_ds, None, rp_refs)
    assert [n for _, n in di.target_candidates(rs_ds)["default"]] == ["PTV_1", "PTV_2"]


def test_rx_candidates_match_resolve_prescription(plan_inputs):
    _, rp_refs, dose = plan_inputs
    one = di.rx_candidates(rp_refs, dose, ["PTV_1"])
    assert one["default"]["target_prescription_dose_gy"] == pytest.approx(20.0)
    assert one["dmax_gy"] == pytest.approx(dose.dmax)
    assert di.resolve_prescription(None, None, rp_refs, dose, ["PTV_1"])[0] == pytest.approx(20.0)
    both = di.rx_candidates(rp_refs, dose, ["PTV_1", "PTV_2"])
    assert both["default"] is None and len(both["refs"]) == 2
    with pytest.raises(ValueError, match="Mehrere Verschreibungen"):
        di.resolve_prescription(None, None, rp_refs, dose, ["PTV_1", "PTV_2"])


def test_roi_table_lists_every_roi_with_its_category(plan_inputs):
    rows = di.roi_table(plan_inputs[0])
    assert len(rows) == 14 and [r[1] for r in rows][:2] == ["BODY", "PTV_1"]
    assert all(r[3] for r in rows)


# -- Eclipse-Referenz --------------------------------------------------------

def test_foreign_target_in_eclipse_ref_is_a_note(demo, tmp_path):
    case = _copy_case(demo, tmp_path)
    (case / "eclipse_ref.json").write_text('{"PTV_1": {"TV": 4.2}, "Fremd": {"TV": 1.0}}',
                                          encoding="utf-8")
    report, art = _run(case, tmp_path)
    assert any("'Fremd'" in n and "nicht ausgewertet" in n for n in report["meta"]["notes"])
    assert art.eclipse["PTV_1"]["n_compared"] >= 1


def test_ambiguous_target_in_eclipse_ref_stays_an_error():
    with pytest.raises(ValueError, match="nicht eindeutig"):
        di.eclipse_reference_from_dict({"PTV": {"TV": 1.0}}, ["PTV_1", "PTV_2"], 20.0, "json")


def test_dvh_of_another_structure_set_is_not_shown(demo, tmp_path):
    case = _copy_case(demo, tmp_path)
    rs_file = next(case.glob("RS*.dcm"))
    ds = pydicom.dcmread(str(rs_file))
    set_sop_instance_uid(ds, "1.2.826.0.1.3680043.8.498.99")
    ds.save_as(str(rs_file))
    report, art = _run(case, tmp_path)
    assert art.eclipse_dvh == {}
    assert any("anderes Structure Set" in n for n in report["meta"]["notes"])


def test_eclipse_ci_and_gi_use_the_whole_isodose_in_component_scope(demo, tmp_path):
    """Eclipse bezieht CI/GI auf die ganze Isodose: bei PIV-Scope component kein Fehlalarm."""
    _report, art = _run(demo.root, tmp_path)
    r, ec = art.targets["PTV_1"].result, art.eclipse["PTV_1"]
    c = r["components"]
    assert c["piv_cm3"] < c["piv_global_cm3"]                     # Teil der Rx-Isodose am PTV_2
    rows = {row["key"]: row for row in ec["rows"]}
    assert rows["ci_paddick"]["tool"] == pytest.approx(c["tv_piv_cm3"] ** 2 / (c["tv_cm3"] * c["piv_global_cm3"]))
    assert rows["gi"]["tool"] == pytest.approx(c["piv50_global_cm3"] / c["piv_global_cm3"])
    assert rows["ci_paddick"]["note"] == "mit globalem PIV" and ec["piv_scope_note"]
    assert r["indices"]["ci_paddick"] > rows["ci_paddick"]["tool"]    # das Tool-Ergebnis bleibt component
    assert not {"ci_paddick", "gi"} & set(ec["flagged"])
    assert not any("CI Paddick" in w or "GI " in w for w in r["warnings"])


# -- --list, --no-rs ---------------------------------------------------------

def test_list_works_without_rtdose(demo, tmp_path, capsys):
    case = _copy_case(demo, tmp_path)
    for f in case.glob("RD*.dcm"):
        f.unlink()
    assert di.main([str(case), "--list"]) == 0
    out = capsys.readouterr().out
    assert "PTV_1" in out and "Verschreibungen" in out


def test_unusable_ct_is_a_warning_without_rs_export(demo, tmp_path):
    case = _copy_case(demo, tmp_path)
    first = sorted((case / "CT").glob("*.dcm"))[0]
    ds = pydicom.dcmread(str(first))
    ds.SeriesInstanceUID = "1.2.826.0.1.3680043.8.498.98"      # zweite Serie im CT-Ordner
    ds.save_as(str(case / "CT" / "zweite_serie.dcm"))
    report, _ = _run(case, tmp_path, write_rs=False)
    assert any("CT-Ordner nicht verwendbar" in w for w in report["warnings"])
    with pytest.raises(ValueError, match="Serien"):
        _run(case, tmp_path / "mit_rs", write_rs=True)


# -- RS-Export ---------------------------------------------------------------

def test_rs_export_leaves_no_tmp_and_keeps_the_full_version(demo, tmp_path):
    report, _ = _run(demo.root, tmp_path)
    rs_path = Path(report["outputs"]["rs_path"])
    assert rs_path.is_file() and not list(rs_path.parent.glob("*.tmp"))
    assert str(pydicom.dcmread(str(rs_path)).SoftwareVersions) == f"{di.TOOL_NAME} {di.TOOL_VERSION}"


def test_failed_rs_check_leaves_no_file(demo, tmp_path, monkeypatch):
    def broken(*args, **kwargs):
        raise RuntimeError("Pruefung nicht moeglich")

    monkeypatch.setattr(rw, "verify_rtstruct", broken)
    report, _ = _run(demo.root, tmp_path)
    out_dir = Path(report["outputs"]["json_path"]).parent
    assert report["outputs"]["rs_path"] is None
    assert not list(out_dir.glob("RS_*"))
    assert any("RS-Export fehlgeschlagen" in w for w in report["warnings"])


# -- Stufen, Fortschritt, Abbruch --------------------------------------------

def test_stages_and_progress_are_reported(demo, tmp_path):
    rec = _Recorder()
    with _runtime.use(rec):
        report, art = _run(demo.root, tmp_path)
    assert rec.stages == ["prepare", "compute", "rs_export", "reports"]
    assert rec.n_progress > 10 and rec.n_checks > 10
    assert art.results is report


def test_cancel_stops_the_run(demo, tmp_path):
    rec = _Recorder(cancel_after=3)
    with _runtime.use(rec), pytest.raises(_runtime.JobCancelled):
        _run(demo.root, tmp_path)
    assert "reports" not in rec.stages
    assert _runtime.current() is not rec


def test_wrapper_returns_the_same_report(demo, tmp_path):
    report = di.run_dose_indices(str(demo.root), output=str(tmp_path / "out"), no_viz=True, quiet=True)
    assert set(report) >= {"meta", "targets", "outputs", "warnings"}
    assert report["meta"]["tool_version"] == di.TOOL_VERSION
