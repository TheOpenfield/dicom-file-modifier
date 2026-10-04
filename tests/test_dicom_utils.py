"""
dicom_utils.contour_points: ContourData direkt aus den Rohbytes, bitgleich zu
pydicom und mit Rueckfall auf pydicom bei allem Unerwarteten.
"""

from __future__ import annotations

import warnings

import numpy as np
import pydicom
import pytest
from pydicom.dataelem import RawDataElement
from pydicom.dataset import Dataset
from pydicom.tag import Tag

from dicom_file_modifier.demo import DemoSpec, make_demo_case
from dicom_file_modifier.dicom_utils import CONTOUR_DATA_TAG, contour_points


@pytest.fixture(scope="module")
def demo(tmp_path_factory):
    return make_demo_case(tmp_path_factory.mktemp("demo") / "case", DemoSpec())


def _contours(ds):
    return [c for rc in ds.ROIContourSequence for c in rc.get("ContourSequence", [])]


def _bits(a: np.ndarray) -> np.ndarray:
    return np.ascontiguousarray(a, dtype=np.float64).view(np.int64)


def test_raw_bytes_give_the_pydicom_values_and_keep_the_dataset_raw(demo):
    fast, ref = pydicom.dcmread(str(demo.rs)), pydicom.dcmread(str(demo.rs))
    pairs = list(zip(_contours(fast), _contours(ref), strict=True))
    assert len(pairs) > 100
    for c_fast, c_ref in pairs:
        got = contour_points(c_fast)
        want = np.asarray(c_ref.ContourData, dtype=np.float64).reshape(-1, 3)
        assert got.dtype == np.float64 and np.array_equal(_bits(got), _bits(want))
        assert isinstance(c_fast.get_item(CONTOUR_DATA_TAG), RawDataElement)     # nichts umgewandelt


def _contour(raw: bytes) -> Dataset:
    c = Dataset()
    tag = Tag(CONTOUR_DATA_TAG)
    c[tag] = RawDataElement(tag, "DS", len(raw), raw, 0, False, True)
    return c


def _outcome(fn):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            return "ok", fn()
        except Exception as exc:  # noqa: BLE001 - Verhalten vergleichen
            return "error", type(exc)


@pytest.mark.parametrize("raw", [
    b"1.5\\-2.25\\3 ",                     # Fuellzeichen am Ende
    b" 1e1\\+2.5\\-0.000000",              # Leerzeichen, Exponent, Vorzeichen
    b"1.5\\2.5",                           # keine Punkt-Tripel
    b"1.5\\abc\\3",                        # ungueltiger Wert
    b"1.5\\\\3",                           # leerer Wert
    b"",                                   # leeres Element
])
def test_unusual_values_behave_like_pydicom(raw):
    got = _outcome(lambda: contour_points(_contour(raw)))
    want = _outcome(lambda: np.asarray(_contour(raw).ContourData, dtype=np.float64).reshape(-1, 3))
    assert got[0] == want[0]
    if got[0] == "ok":
        assert np.array_equal(_bits(got[1]), _bits(want[1]))
    else:
        assert got[1] is want[1]


def test_converted_and_missing_elements(demo):
    ds = pydicom.dcmread(str(demo.rs))
    c = _contours(ds)[0]
    want = np.asarray(c.ContourData, dtype=np.float64).reshape(-1, 3)    # jetzt umgewandelt
    assert np.array_equal(_bits(contour_points(c)), _bits(want))
    with pytest.raises(AttributeError):
        contour_points(Dataset())
