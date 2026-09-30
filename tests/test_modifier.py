"""
modifier (Plan P0.6): CT-Loader mit Modality-/Serienfilter, run_ct_transform,
Fortschritt/Abbruch im Resampling, kein Browserfenster, Exit-Codes, Hilfetexte.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pydicom
import pytest

from dicom_file_modifier import _runtime
from dicom_file_modifier import modifier as mod
from dicom_file_modifier.demo import DemoSpec, make_demo_case
from dicom_file_modifier.issues import UserInputError


@pytest.fixture(scope="module")
def flat(tmp_path_factory):
    """Flacher Eclipse-Export: CT, RS, RP und RD in einem Ordner."""
    return make_demo_case(tmp_path_factory.mktemp("flat") / "case", DemoSpec(layout="flat"))


class _Recorder(_runtime.JobContext):
    def __init__(self, cancel_after=None):
        self.stages, self.n_progress, self.n_checks, self.cancel_after = [], 0, 0, cancel_after

    def stage(self, key, label):
        self.stages.append(key)

    def progress(self, done, total, text=""):
        self.n_progress += 1

    def check_cancel(self):
        self.n_checks += 1
        if self.cancel_after is not None and self.n_checks > self.cancel_after:
            raise _runtime.JobCancelled()


def test_loader_keeps_only_the_ct_series(flat):
    files = mod.ct_dir_files(flat.root)
    slices = mod.load_ct_series_files(files)
    assert len(slices) == 80 < len(files)
    assert {str(s.Modality) for s in slices} == {"CT"}
    z = [float(s.ImagePositionPatient[2]) for s in slices]
    assert z == sorted(z)
    assert len(mod.load_ct_headers(flat.root)) == 80


def test_two_series_need_a_choice(flat, tmp_path):
    for f in mod.ct_dir_files(flat.root):
        shutil.copy(f, tmp_path)
    second = pydicom.dcmread(str(flat.ct_files[0]))
    second.SeriesInstanceUID = "1.2.826.0.1.3680043.8.498.96"
    second.SOPInstanceUID = second.file_meta.MediaStorageSOPInstanceUID = "1.2.826.0.1.3680043.8.498.95"
    second.save_as(str(tmp_path / "zweite_serie.dcm"))
    with pytest.raises(UserInputError) as err:
        mod.load_ct_headers(tmp_path)
    assert err.value.issue.code == "CT.MULTIPLE_SERIES"
    first_uid = str(pydicom.dcmread(str(flat.ct_files[1]), stop_before_pixels=True).SeriesInstanceUID)
    assert len(mod.load_ct_headers(tmp_path, series_uid=first_uid)) == 80


def test_run_ct_transform_on_a_flat_export(flat, tmp_path):
    res = mod.run_ct_transform(flat.root, tmp_path / "out", tx=3, method="metadata", viz=False)
    assert res["n_slices"] == 80 and len(list((tmp_path / "out").glob("*.dcm"))) == 80
    assert res["viz_html_path"] is None and len(res["sop_map"]) == 80


def test_resample_reports_progress_and_can_be_cancelled(flat, tmp_path):
    rec = _Recorder()
    with _runtime.use(rec):
        mod.run_ct_transform(flat.root, tmp_path / "a", rz=5, viz=False)
    assert rec.stages == ["load", "transform"] and rec.n_progress >= 5   # 80 Schichten / 20 je Chunk + Ende

    rec = _Recorder(cancel_after=2)                   # im Resampling, vor dem Speichern
    with _runtime.use(rec), pytest.raises(_runtime.JobCancelled):
        mod.run_ct_transform(flat.root, tmp_path / "b", rz=5, viz=False)
    assert not (tmp_path / "b").exists()


def test_visualize_3d_never_opens_a_browser(monkeypatch):
    import plotly.graph_objects as go

    def forbidden(self, *args, **kwargs):
        raise AssertionError("fig.show() im API-Pfad")

    monkeypatch.setattr(go.Figure, "show", forbidden)
    vol = np.full((12, 24, 24), -1000.0, dtype=np.float32)
    vol[3:9, 6:18, 6:18] = 40.0
    geom = {"affine": np.eye(4), "shape": vol.shape}
    fig = mod.visualize_3d(vol, vol, geom, np.eye(4), method="metadata", output_html=None)
    assert isinstance(fig, go.Figure) and len(fig.data) >= 2


def test_feet_first_input_is_exit_2(flat, tmp_path, capsys):
    ct_dir = tmp_path / "ct"
    ct_dir.mkdir()
    for f in flat.ct_files:
        ds = pydicom.dcmread(str(f))
        ds.ImageOrientationPatient = [-1, 0, 0, 0, 1, 0]
        ds.save_as(str(ct_dir / Path(f).name))
    assert mod.main([str(ct_dir), "--no-viz", "-o", str(tmp_path / "out")]) == 2
    assert "Feet-First" in capsys.readouterr().err


def test_help_names_yaw_and_roll_for_a_supine_patient(capsys):
    with pytest.raises(SystemExit):
        mod.main(["--help"])
    lines = capsys.readouterr().out.splitlines()

    def help_for(opt: str) -> str:            # Optionszeile (nicht die Usage) plus Umbruch
        i = next(k for k, ln in enumerate(lines) if ln.strip().startswith(opt + " "))
        return " ".join(lines[i:i + 2])

    assert "Pitch" in help_for("--rx") and "Yaw" in help_for("--ry") and "Roll" in help_for("--rz")
