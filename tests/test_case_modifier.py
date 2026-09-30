"""
case_modifier (Plan P0.5): Stufen, Header-Pruefung, nichts Halbes auf der
Platte, Label- und Orientierungspruefung, Issues, Verify, Zentrumswahl.
Alle Laeufe auf dem synthetischen Demo-Fall.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pydicom
import pytest

from dicom_file_modifier import _runtime
from dicom_file_modifier import case_modifier as cm
from dicom_file_modifier import modifier as mod
from dicom_file_modifier.demo import DemoSpec, make_demo_case
from dicom_file_modifier.issues import UserInputError


@pytest.fixture(scope="module")
def demo(tmp_path_factory):
    return make_demo_case(tmp_path_factory.mktemp("demo") / "case", DemoSpec())


def _copy_case(demo, tmp_path: Path) -> Path:
    dst = tmp_path / "case"
    shutil.copytree(demo.root, dst)
    return dst


def _run(case_dir, out_dir: Path, **kw):
    kw.setdefault("method", "metadata")
    kw.setdefault("no_viz", True)
    return cm.run_case_transform(str(case_dir), str(out_dir), 3, -2, 1, 0, 0, 5, **kw)


def _codes(result: dict) -> set:
    return {i["code"] for i in result["issues"]}


# -- Geometrie- und Orientierungspruefung (nur Header) ------------------------

def _fake_slices(iop, step, n=4):
    return [SimpleNamespace(ImageOrientationPatient=list(iop), PixelSpacing=[1.0, 1.0],
                            ImagePositionPatient=list(np.asarray(step, float) * k)) for k in range(n)]


@pytest.mark.parametrize("iop", [(1, 0, 0, 0, 1, 0), (-1, 0, 0, 0, -1, 0)])   # HFS, HFP
def test_head_first_passes_the_orientation_guard(iop):
    cm.validate_ct_geometry(_fake_slices(iop, (0, 0, 2.0)))


@pytest.mark.parametrize("iop, step", [
    ((-1, 0, 0, 0, 1, 0), (0, 0, 2.0)),                                        # FFS
    ((1, 0, 0, 0, 1, 0), (0, 2.0 * np.sin(np.radians(10)), 2.0 * np.cos(np.radians(10)))),  # Kippung 10 Grad
])
def test_feet_first_and_tilt_are_rejected(iop, step):
    with pytest.raises(UserInputError) as err:
        cm.validate_ct_geometry(_fake_slices(iop, step))
    assert err.value.issue.code == "CT.ORIENTATION_UNSUPPORTED"


def test_feet_first_case_is_rejected_before_anything_is_written(demo, tmp_path):
    case = _copy_case(demo, tmp_path)
    for f in (case / "CT").glob("*.dcm"):
        ds = pydicom.dcmread(str(f))
        ds.ImageOrientationPatient = [-1, 0, 0, 0, 1, 0]
        ds.save_as(str(f))
    with pytest.raises(UserInputError, match="Feet-First"):
        _run(case, tmp_path / "out")
    assert not (tmp_path / "out").exists()


# -- Label --------------------------------------------------------------------

@pytest.mark.parametrize("label, case_id", [("_a/b", "case"), ("_x:y", "case"), ("_x.", "case"), ("", "CON")])
def test_invalid_labels_are_rejected(label, case_id):
    with pytest.raises(UserInputError) as err:
        cm.validate_label(label, case_id)
    assert err.value.issue.code == "CASE.LABEL_INVALID" and err.value.issue.field == "label"


def test_valid_labels_pass():
    for label in ("_RB", "_NF-2", "_v1.2", ""):
        cm.validate_label(label, "case")


# -- Dry-Run, Planung vor dem Schreiben ---------------------------------------

def test_dry_run_reads_no_pixels(demo, tmp_path, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("Dry-Run darf keine Pixel laden")

    monkeypatch.setattr(mod, "load_ct_series_files", forbidden)
    monkeypatch.setattr(mod, "slices_to_hu", forbidden)
    res = _run(demo.root, tmp_path / "out", dry_run=True)
    assert res["dry_run"] and res["planned_rs_output_path"].endswith("RS_RB.dcm")
    assert not (tmp_path / "out").exists()


def test_rs_reference_to_unknown_ct_fails_before_writing(demo, tmp_path):
    case = _copy_case(demo, tmp_path)
    rs_file = next(case.glob("RS*.dcm"))
    ds = pydicom.dcmread(str(rs_file))
    for rc in ds.ROIContourSequence:
        for c in rc.get("ContourSequence", []):
            if "ContourImageSequence" in c:
                c.ContourImageSequence[0].ReferencedSOPInstanceUID = "1.2.826.0.1.3680043.8.498.97"
                break
        else:
            continue
        break
    ds.save_as(str(rs_file))
    with pytest.raises(KeyError, match="unbekannte CT-SOPInstanceUID"):
        _run(case, tmp_path / "out")
    assert not (tmp_path / "out").exists()


# -- Bildverweise der Konturen ------------------------------------------------

def _plan(demo, tz=0.0, rx=0.0, method="resample"):
    pre = cm.preflight_case(str(demo.root), quiet=True)
    return pre, cm.plan_transform(pre, 3.0, -2.0, tz, rx, 0.0, 5.0, out_dir=".", method=method, quiet=True)


def _off_plane(plan):
    return next((i for i in plan.issues if i.code == "CASE.CONTOURS_OFF_PLANE"), None)


def test_resample_points_each_contour_to_the_slice_at_its_new_z(demo):
    pre, plan = _plan(demo, tz=2.0)                              # 2 Schichten nach superior
    z_old = {str(h.SOPInstanceUID): float(h.ImagePositionPatient[2]) for h in pre.ct_headers}
    z_of_new = {new: z_old[old] for old, new in plan.sop_map.items()}   # das Raster bleibt
    lo, hi = min(z_old.values()), max(z_old.values())
    pairs = []
    for rc in plan.new_rs.ROIContourSequence:
        for c in rc.get("ContourSequence", []):
            if "ContourImageSequence" in c and str(c.ContourGeometricType) != "POINT":
                z = float(np.mean(np.asarray(c.ContourData, dtype=float).reshape(-1, 3)[:, 2]))
                if lo <= z <= hi:
                    pairs.append((z, z_of_new[str(c.ContourImageSequence[0].ReferencedSOPInstanceUID)]))
    assert len(pairs) > 100 and all(abs(z - zr) < 1e-6 for z, zr in pairs)
    assert _off_plane(plan) is None


@pytest.mark.parametrize("tz, rx, tilted", [(0.5, 0.0, False), (0.0, 2.0, True)])
def test_resample_warns_when_contours_leave_the_slice_planes(demo, tz, rx, tilted):
    _, plan = _plan(demo, tz=tz, rx=rx)
    issue = _off_plane(plan)
    assert issue is not None and issue.level == "warning"
    assert ("gekippt" in issue.message_de) == tilted
    assert ("metadata" in issue.hint_de) == tilted and ("Vielfaches" in issue.hint_de) != tilted


def test_metadata_keeps_the_slice_of_each_contour(demo):
    pre, plan = _plan(demo, tz=0.5, rx=2.0, method="metadata")  # die Schichten wandern mit
    assert _off_plane(plan) is None
    n = 0
    for rc_old, rc_new in zip(pre.rs_ds.ROIContourSequence, plan.new_rs.ROIContourSequence,
                              strict=False):                    # neu: zusaetzlich der Drehpunkt
        for c_old, c_new in zip(rc_old.get("ContourSequence", []), rc_new.get("ContourSequence", []),
                                strict=True):
            if "ContourImageSequence" in c_old:
                n += 1
                assert str(c_new.ContourImageSequence[0].ReferencedSOPInstanceUID) == \
                    plan.sop_map[str(c_old.ContourImageSequence[0].ReferencedSOPInstanceUID)]
    assert n > 100


# -- Ausfuehren ---------------------------------------------------------------

def test_metadata_run_writes_series_number_offset_in_one_pass(demo, tmp_path):
    res = _run(demo.root, tmp_path / "out", verify=True)
    first = sorted(Path(res["ct_output_dir"]).glob("*.dcm"))[0]
    orig = pydicom.dcmread(str(demo.ct_files[0]), stop_before_pixels=True)
    assert int(pydicom.dcmread(str(first)).SeriesNumber) == int(orig.SeriesNumber) + 1000
    assert res["verify"]["passed"] and res["verify"]["threshold_mm"] == cm.VERIFY_THRESHOLD_MM
    assert {"CASE.SIBLINGS_NOT_TRANSFORMED", "CASE.FOR_KEPT"} <= _codes(res)
    assert next(i for i in res["issues"] if i["code"] == "CASE.FOR_KEPT")["level"] == "info"   # gewollt
    assert res["for_strategy"] == "keep" and res["clipping"] == []


def test_new_frame_of_reference_has_no_for_warning(demo, tmp_path):
    res = _run(demo.root, tmp_path / "out", new_frame_of_reference=True)
    assert "CASE.FOR_KEPT" not in _codes(res) and res["for_strategy"] == "new"


def test_discover_case_reports_siblings_quietly(demo, capsys):
    ct_dir, rs, sib = cm.discover_case(str(demo.root))
    assert ct_dir.is_dir() and rs.is_file()
    assert {p.name[:2] for p in sib} == {"RP", "RD"} and capsys.readouterr().out == ""


# -- Zentrum ------------------------------------------------------------------

def test_resolve_center(demo, monkeypatch):
    pre_rs = pydicom.dcmread(str(demo.rs))
    vol = np.zeros(3)
    assert cm.resolve_center("volume", pre_rs, vol) == (None, "Volumenmitte")
    pos, label = cm.resolve_center("marker:HS1", pre_rs, vol)
    assert label == "Marker 'HS1'" and np.allclose(pos, (10.0, 5.0, -15.5))
    pos, label = cm.resolve_center("1,2,3", pre_rs, vol)
    assert label == "manuell" and np.allclose(pos, (1, 2, 3))
    monkeypatch.setattr("sys.stdin", None)            # z.B. pythonw / GUI-EXE
    assert cm.resolve_center(None, pre_rs, vol) == (None, "Volumenmitte")
    assert cm.resolve_center(None, pre_rs, vol, interactive=False) == (None, "Volumenmitte")


# -- Stufen und Abbruch -------------------------------------------------------

class _Recorder(_runtime.JobContext):
    def __init__(self, cancel_after=None):
        self.stages, self.n_checks, self.cancel_after = [], 0, cancel_after

    def stage(self, key, label):
        self.stages.append(key)

    def check_cancel(self):
        self.n_checks += 1
        if self.cancel_after is not None and self.n_checks > self.cancel_after:
            raise _runtime.JobCancelled()


def test_stages_are_reported(demo, tmp_path):
    rec = _Recorder()
    with _runtime.use(rec):
        _run(demo.root, tmp_path / "out", verify=True)
    assert rec.stages == ["preflight", "plan", "ct_transform", "rs_write", "verify"]


def test_cancel_before_execution_writes_nothing(demo, tmp_path):
    rec = _Recorder(cancel_after=1)                   # Planung laeuft, Ausfuehrung nicht
    with _runtime.use(rec), pytest.raises(_runtime.JobCancelled):
        _run(demo.root, tmp_path / "out")
    assert rec.stages == ["preflight", "plan"] and not (tmp_path / "out").exists()
