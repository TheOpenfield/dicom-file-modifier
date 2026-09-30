"""
demo.py - Synthetischer Demo- und Testfall (CT + RTSTRUCT + RTPLAN + RTDOSE).

Erzeugt einen vollstaendig kuenstlichen Kopf-Phantom-Fall (Head-First-Supine):
Wasserzylinder mit Knochenschale, zwei kugelfoermige PTVs mit getrennten
Isodosen, Risikoorgane in bekannten Abstaenden, Hilfsstrukturen, POINT-Marker
und eine leere ROI; dazu ein RTPLAN mit zwei Verschreibungen und eine RTDOSE
mit analytischer Dosis (Hill-Profil je Ziel) und Eclipse-artiger DVHSequence.
Kein Patient, keine echten Daten: der Fall darf in Tests, CI und Doku stehen.

Die UIDs sind deterministisch (gleicher ``seed`` -> identische Dateien), die
analytischen Sollwerte landen in ``demo_expected.json``.

Verwendung:
  python -m dicom_file_modifier.demo <ausgabeordner>
  python -m dicom_file_modifier.demo <ausgabeordner> --layout flat --rd-set plan+2beams
  python -m dicom_file_modifier.demo <ausgabeordner> --rle --no-dvh --eclipse-ref both
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import pydicom
from pydicom.dataset import Dataset, FileMetaDataset
from pydicom.sequence import Sequence
from pydicom.tag import Tag
from pydicom.uid import (PYDICOM_IMPLEMENTATION_UID, PYDICOM_ROOT_UID, ExplicitVRLittleEndian,
                         ImplicitVRLittleEndian, RLELossless, generate_uid)

from . import analyzer as ana
from . import phantom as ph

DEMO_VERSION = "1"

CT_SOP_CLASS = "1.2.840.10008.5.1.4.1.1.2"
RTSTRUCT_SOP_CLASS = "1.2.840.10008.5.1.4.1.1.481.3"
RTPLAN_SOP_CLASS = "1.2.840.10008.5.1.4.1.1.481.5"
RTDOSE_SOP_CLASS = "1.2.840.10008.5.1.4.1.1.481.2"
STUDY_REF_SOP_CLASS = "1.2.840.10008.3.1.2.3.1"     # Detached Study Management (wie Eclipse)

DATE, TIME = "20260101", "120000"
N_POLY = 360

# Geometrie (LPS, mm) -------------------------------------------------------
BODY_R, WATER_R, BRAIN_R = 85.0, 80.0, 70.0
PTV1_C, PTV1_R, GTV1_R, RX1 = (20.0, -10.0, 0.0), 10.0, 8.0, 20.0
PTV2_C, PTV2_R, RX2 = (-25.0, 15.0, 10.0), 6.0, 18.0
HIRNSTAMM = dict(cx=20.0, cy=14.0, r=7.0, z_min=-39.5, z_max=9.5)
CHIASMA = dict(c=(5.0, -30.0, 5.0), r=4.0)
RUECKENMARK = dict(cx=0.0, cy=40.0, r=5.0, z_min=-39.5, z_max=-20.5)
HIRN_HOLE = dict(cx=0.0, cy=0.0, r=8.0, z_min=-5.5, z_max=5.5)
AUGEN = ((32.0, -62.0, 0.0), (-32.0, -62.0, 0.0), 8.0)
RING = dict(r_out=20.0, r_in=13.0, z_min=-11.5, z_max=11.5)
HS1_POS = (10.0, 5.0, -15.5)
ISO_POS = PTV1_C

SIZES = {                       # (Pixel, Pixelabstand mm)
    "small": (128, 2.0),
    "medium": (256, 1.0),
}
N_SLICES, SLICE_DZ = 80, 1.0
DOSE_RES, DOSE_NX, DOSE_NY = 1.0, 110, 90   # RTDOSE-Gitter 1 mm, x -54.5..54.5, y -44.5..44.5
DOSE_SCALING = 1e-6                          # Gy pro Pixelwert (uint32)
DVH_BIN_GY = 0.01


@dataclass
class DemoSpec:
    """Varianten des Demo-Falls (Default = Standardfall mit CT/-Unterordner)."""
    layout: str = "ct_subdir"          # ct_subdir | flat (Eclipse-Export: alles in einem Ordner)
    rd_set: str = "plan"               # plan | plan+2beams | 2plans
    dvh: bool = True                   # DVHSequence in der RTDOSE
    explicit_vr: bool = False          # Explicit statt Implicit VR Little Endian
    compress_ct: Optional[str] = None  # None | "rle"
    eclipse_ref: Optional[str] = None  # None | "ptv1" | "both" (eclipse_ref.json im Fallordner)
    empty_roi: bool = True             # ROI ohne Konturen
    size: str = "small"                # small (128^2 x 2 mm) | medium (256^2 x 1 mm)
    seed: Optional[str] = "dfm-demo"   # None -> zufaellige UIDs

    def validate(self) -> None:
        if self.layout not in ("ct_subdir", "flat"):
            raise ValueError(f"layout {self.layout!r}: erlaubt ct_subdir | flat")
        if self.rd_set not in ("plan", "plan+2beams", "2plans"):
            raise ValueError(f"rd_set {self.rd_set!r}: erlaubt plan | plan+2beams | 2plans")
        if self.compress_ct not in (None, "rle"):
            raise ValueError(f"compress_ct {self.compress_ct!r}: erlaubt None | rle")
        if self.eclipse_ref not in (None, "ptv1", "both"):
            raise ValueError(f"eclipse_ref {self.eclipse_ref!r}: erlaubt None | ptv1 | both")
        if self.size not in SIZES:
            raise ValueError(f"size {self.size!r}: erlaubt {', '.join(SIZES)}")


@dataclass
class DemoCase:
    """Pfade und Eckdaten eines erzeugten Demo-Falls."""
    root: Path
    ct_dir: Path
    ct_files: list
    rs: Path
    rp: Path
    rd: list
    eclipse_ref: Optional[Path]
    expected_json: Path
    spec: DemoSpec
    uids: dict = field(default_factory=dict)


@dataclass
class _Roi:
    number: int
    name: str
    rt_type: str
    color: tuple
    contours: list
    point: bool = False
    algorithm: str = "MANUAL"


# ---------------------------------------------------------------------------
# Hilfen
# ---------------------------------------------------------------------------

class _Uids:
    """Deterministische UIDs aus ``seed`` + Rolle (oder zufaellig ohne seed)."""

    def __init__(self, seed: Optional[str]):
        self.seed = seed

    def __call__(self, *parts) -> str:
        if self.seed is None:
            return str(generate_uid())
        return str(generate_uid(prefix=PYDICOM_ROOT_UID,
                                entropy_srcs=[str(self.seed)] + [str(p) for p in parts]))


def _ds(v: float) -> str:
    """DS-String (<= 16 Zeichen), identisch auf allen pydicom-Versionen."""
    s = format(float(v), ".8g")
    return "0" if s in ("-0", "0") else s


def _flat_ds(pts: np.ndarray) -> list:
    return [_ds(v) for v in np.asarray(pts, float).ravel()]


def _base(sop_class: str, sop_uid: str, modality: str, ctx: dict) -> Dataset:
    ds = Dataset()
    fm = FileMetaDataset()
    fm.MediaStorageSOPClassUID = sop_class
    fm.MediaStorageSOPInstanceUID = sop_uid
    fm.TransferSyntaxUID = ctx["ts"]
    fm.ImplementationClassUID = PYDICOM_IMPLEMENTATION_UID
    fm.ImplementationVersionName = f"DFM_DEMO_{DEMO_VERSION}"
    ds.file_meta = fm
    ds.SpecificCharacterSet = "ISO_IR 100"
    ds.InstanceCreationDate, ds.InstanceCreationTime = DATE, TIME
    ds.SOPClassUID = sop_class
    ds.SOPInstanceUID = sop_uid
    ds.StudyDate, ds.StudyTime = DATE, TIME
    ds.AccessionNumber = ""
    ds.Modality = modality
    ds.Manufacturer = "dicom_file_modifier"
    ds.ReferringPhysicianName = ""
    ds.StudyDescription = "SYNTHETISCH - kein Patient"
    ds.ManufacturerModelName = "demo"
    ds.PatientName = "DFM^DEMO"
    ds.PatientID = "DFM-DEMO"
    ds.PatientBirthDate = ""
    ds.PatientSex = "O"
    ds.SoftwareVersions = f"dfm demo {DEMO_VERSION}"
    ds.StudyInstanceUID = ctx["study_uid"]
    ds.StudyID = "1"
    return ds


def _write(ds: Dataset, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    pydicom.dcmwrite(str(path), ds, enforce_file_format=True)
    return path


# ---------------------------------------------------------------------------
# Geometrie und Dosis
# ---------------------------------------------------------------------------

def _planes() -> np.ndarray:
    return (np.arange(N_SLICES) - (N_SLICES - 1) / 2.0) * SLICE_DZ     # -39.5 .. 39.5


def _fields() -> tuple:
    f1 = ph.HillField.for_target(PTV1_C, PTV1_R, RX1)
    f2 = ph.HillField.for_target(PTV2_C, PTV2_R, RX2)
    return f1, f2


def _build_rois(planes: np.ndarray, empty_roi: bool) -> list:
    sph = ph.sphere_contours
    ptv1 = sph(PTV1_C, PTV1_R, planes, N_POLY)
    ptv2 = sph(PTV2_C, PTV2_R, planes, N_POLY)
    hirn = ph.cylinder_contours(0.0, 0.0, BRAIN_R, planes, -34.5, 34.5, N_POLY)
    hirn += ph.cylinder_contours(HIRN_HOLE["cx"], HIRN_HOLE["cy"], HIRN_HOLE["r"], planes,
                                 HIRN_HOLE["z_min"], HIRN_HOLE["z_max"], N_POLY)
    augen = sph(AUGEN[0], AUGEN[2], planes, N_POLY) + sph(AUGEN[1], AUGEN[2], planes, N_POLY)
    hs, rm, ch = HIRNSTAMM, RUECKENMARK, CHIASMA
    rois = [
        _Roi(1, "BODY", "EXTERNAL", (0, 255, 0),
             ph.cylinder_contours(0.0, 0.0, BODY_R, planes, planes[0], planes[-1], N_POLY),
             algorithm="AUTOMATIC"),
        _Roi(2, "PTV_1", "PTV", (255, 0, 0), ptv1),
        _Roi(3, "GTV_1", "GTV", (255, 128, 0), sph(PTV1_C, GTV1_R, planes, N_POLY)),
        _Roi(4, "PTV_2", "PTV", (255, 0, 128), ptv2),
        _Roi(5, "Hirnstamm", "ORGAN", (0, 0, 255),
             ph.cylinder_contours(hs["cx"], hs["cy"], hs["r"], planes, hs["z_min"], hs["z_max"], N_POLY)),
        _Roi(6, "Chiasma", "ORGAN", (0, 255, 255), sph(ch["c"], ch["r"], planes, N_POLY)),
        _Roi(7, "R\u00fcckenmark", "ORGAN", (255, 255, 0),
             ph.cylinder_contours(rm["cx"], rm["cy"], rm["r"], planes, rm["z_min"], rm["z_max"], N_POLY)),
        _Roi(8, "Hirn", "ORGAN", (200, 200, 200), hirn),
        _Roi(9, "Augen", "ORGAN", (0, 128, 255), augen),
        _Roi(10, "h_PTV_gesamt", "PTV", (128, 0, 0), ptv1 + ptv2, algorithm="AUTOMATIC"),
        _Roi(11, "Ring_PTV", "CONTROL", (128, 128, 0),
             ph.ring_contours(PTV1_C[0], PTV1_C[1], RING["r_out"], RING["r_in"], planes,
                              RING["z_min"], RING["z_max"], N_POLY), algorithm="AUTOMATIC"),
        _Roi(12, "HS1", "MARKER", (255, 0, 255), [np.asarray([HS1_POS], float)], point=True),
        _Roi(13, "Iso", "ISOCENTER", (255, 255, 255), [np.asarray([ISO_POS], float)], point=True),
    ]
    if empty_roi:
        rois.append(_Roi(14, "Leer", "ORGAN", (100, 100, 100), []))
    return rois


def _ct_image(n: int, ps: float, z: float) -> np.ndarray:
    """HU-Bild einer Schicht: Luft, Knochenschale, Wasser, Hirnbereich (HU 35)."""
    c = (np.arange(n) - (n - 1) / 2.0) * ps
    X, Y = np.meshgrid(c, c)
    r = np.hypot(X, Y)
    img = np.full((n, n), -1000.0)
    img[r <= BODY_R] = 1000.0
    img[r <= WATER_R] = 0.0
    if abs(z) <= 34.5:
        img[r <= BRAIN_R] = 35.0
    return img


# ---------------------------------------------------------------------------
# Objekte bauen
# ---------------------------------------------------------------------------

def _ct_slices(ctx: dict, spec: DemoSpec) -> list:
    n, ps = SIZES[spec.size]
    x0 = -(n - 1) / 2.0 * ps
    out = []
    for k, z in enumerate(ctx["planes"]):
        sop = ctx["uid"]("ct", k)
        ds = _base(CT_SOP_CLASS, sop, "CT", ctx)
        ds.ImageType = ["ORIGINAL", "PRIMARY", "AXIAL"]
        ds.SeriesInstanceUID = ctx["ct_series"]
        ds.SeriesNumber = 1
        ds.SeriesDescription = "DEMO CT"
        ds.FrameOfReferenceUID = ctx["for_uid"]
        ds.PositionReferenceIndicator = ""
        ds.InstanceNumber = k + 1
        ds.PatientPosition = "HFS"
        ds.KVP = "120"
        ds.SliceThickness = _ds(SLICE_DZ)
        ds.ImagePositionPatient = [_ds(x0), _ds(x0), _ds(z)]
        ds.ImageOrientationPatient = ["1", "0", "0", "0", "1", "0"]
        ds.SliceLocation = _ds(z)
        ds.SamplesPerPixel = 1
        ds.PhotometricInterpretation = "MONOCHROME2"
        ds.Rows = ds.Columns = n
        ds.PixelSpacing = [_ds(ps), _ds(ps)]
        ds.BitsAllocated, ds.BitsStored, ds.HighBit, ds.PixelRepresentation = 16, 12, 11, 0
        ds.WindowCenter, ds.WindowWidth = "40", "400"
        ds.RescaleIntercept, ds.RescaleSlope, ds.RescaleType = "-1024", "1", "HU"
        stored = (_ct_image(n, ps, float(z)) + 1024.0).astype(np.uint16)
        ds.PixelData = stored.tobytes()
        if spec.compress_ct == "rle":
            ds.file_meta.TransferSyntaxUID = ExplicitVRLittleEndian
            ds.compress(RLELossless, arr=stored)
        out.append(ds)
    return out


def _ct_plane_sop(ctx: dict, z: float) -> Optional[str]:
    k = int(np.argmin(np.abs(ctx["planes"] - z)))
    if abs(ctx["planes"][k] - z) < 1e-3:
        return ctx["ct_sops"][k]
    return None


def _rtstruct(ctx: dict, rois: list) -> Dataset:
    ds = _base(RTSTRUCT_SOP_CLASS, ctx["rs_sop"], "RTSTRUCT", ctx)
    ds.SeriesInstanceUID = ctx["rs_series"]
    ds.SeriesNumber = 2
    ds.SeriesDescription = "DEMO RS"
    ds.StructureSetLabel = "DEMO"
    ds.StructureSetName = "DEMO"
    ds.StructureSetDate, ds.StructureSetTime = DATE, TIME

    ci_all = Sequence()
    for sop in ctx["ct_sops"]:
        ci = Dataset()
        ci.ReferencedSOPClassUID = CT_SOP_CLASS
        ci.ReferencedSOPInstanceUID = sop
        ci_all.append(ci)
    se = Dataset()
    se.SeriesInstanceUID = ctx["ct_series"]
    se.ContourImageSequence = ci_all
    st = Dataset()
    st.ReferencedSOPClassUID = STUDY_REF_SOP_CLASS
    st.ReferencedSOPInstanceUID = ctx["study_uid"]
    st.RTReferencedSeriesSequence = Sequence([se])
    rf = Dataset()
    rf.FrameOfReferenceUID = ctx["for_uid"]
    rf.RTReferencedStudySequence = Sequence([st])
    ds.ReferencedFrameOfReferenceSequence = Sequence([rf])

    ssr_seq, rc_seq, obs_seq = Sequence(), Sequence(), Sequence()
    for roi in rois:
        ssr = Dataset()
        ssr.ROINumber = roi.number
        ssr.ReferencedFrameOfReferenceUID = ctx["for_uid"]
        ssr.ROIName = roi.name
        ssr.ROIGenerationAlgorithm = roi.algorithm
        ssr_seq.append(ssr)

        rc = Dataset()
        rc.ROIDisplayColor = [int(v) for v in roi.color]
        if roi.contours:
            cseq = Sequence()
            for pts in roi.contours:
                c = Dataset()
                sop = _ct_plane_sop(ctx, float(pts[0, 2]))
                if sop is not None:
                    ci = Dataset()
                    ci.ReferencedSOPClassUID = CT_SOP_CLASS
                    ci.ReferencedSOPInstanceUID = sop
                    c.ContourImageSequence = Sequence([ci])
                c.ContourGeometricType = "POINT" if roi.point else "CLOSED_PLANAR"
                c.NumberOfContourPoints = int(len(pts))
                c.ContourData = _flat_ds(pts)
                cseq.append(c)
            rc.ContourSequence = cseq
        rc.ReferencedROINumber = roi.number
        rc_seq.append(rc)

        obs = Dataset()
        obs.ObservationNumber = roi.number
        obs.ReferencedROINumber = roi.number
        obs.ROIObservationLabel = roi.name[:16]
        obs.RTROIInterpretedType = roi.rt_type
        obs.ROIInterpreter = ""
        obs_seq.append(obs)
    ds.StructureSetROISequence = ssr_seq
    ds.ROIContourSequence = rc_seq
    ds.RTROIObservationsSequence = obs_seq
    ds.ApprovalStatus = "UNAPPROVED"
    return ds


def _rtplan(ctx: dict, with_beams: bool) -> Dataset:
    ds = _base(RTPLAN_SOP_CLASS, ctx["rp_sop"], "RTPLAN", ctx)
    ds.SeriesInstanceUID = ctx["rp_series"]
    ds.SeriesNumber = 3
    ds.FrameOfReferenceUID = ctx["for_uid"]
    ds.PositionReferenceIndicator = ""
    ds.RTPlanLabel = "DEMO_PLAN"
    ds.RTPlanName = "DEMO"
    ds.RTPlanDate, ds.RTPlanTime = DATE, TIME
    ds.RTPlanGeometry = "PATIENT"

    rs_ref = Dataset()
    rs_ref.ReferencedSOPClassUID = RTSTRUCT_SOP_CLASS
    rs_ref.ReferencedSOPInstanceUID = ctx["rs_sop"]
    ds.ReferencedStructureSetSequence = Sequence([rs_ref])

    drs = Sequence()
    for num, (name, roi_num, rx) in enumerate((("PTV_1", 2, RX1), ("PTV_2", 4, RX2)), start=1):
        dr = Dataset()
        dr.DoseReferenceNumber = num
        dr.DoseReferenceUID = ctx["uid"]("doseref", num)
        dr.DoseReferenceStructureType = "VOLUME"
        dr.DoseReferenceDescription = name
        dr.ReferencedROINumber = roi_num
        dr.DoseReferenceType = "TARGET"
        dr.TargetPrescriptionDose = _ds(rx)
        drs.append(dr)
    ds.DoseReferenceSequence = drs

    fg = Dataset()
    fg.FractionGroupNumber = 1
    fg.NumberOfFractionsPlanned = 1
    fg.NumberOfBeams = 2 if with_beams else 0
    fg.NumberOfBrachyApplicationSetups = 0
    if with_beams:
        refs = Sequence()
        for b in (1, 2):
            rb = Dataset()
            rb.ReferencedBeamNumber = b
            refs.append(rb)
        fg.ReferencedBeamSequence = refs
        beams = Sequence()
        for b in (1, 2):
            be = Dataset()
            be.BeamNumber = b
            be.BeamName = f"Feld {b}"
            be.BeamType = "STATIC"
            be.RadiationType = "PHOTON"
            be.TreatmentMachineName = "DEMO"
            be.NumberOfWedges = 0
            be.NumberOfCompensators = 0
            be.NumberOfBoli = 0
            be.NumberOfBlocks = 0
            be.NumberOfControlPoints = 0
            be.ReferencedPatientSetupNumber = 1
            beams.append(be)
        ds.BeamSequence = beams
    ds.FractionGroupSequence = Sequence([fg])

    ps_item = Dataset()
    ps_item.PatientPosition = "HFS"
    ps_item.PatientSetupNumber = 1
    ds.PatientSetupSequence = Sequence([ps_item])
    ds.ApprovalStatus = "UNAPPROVED"
    return ds


def _dvh_item(roi_number: int, dvh: dict) -> Dataset:
    it = Dataset()
    ref = Dataset()
    ref.ReferencedROINumber = roi_number
    ref.DVHROIContributionType = "INCLUDED"
    it.DVHReferencedROISequence = Sequence([ref])
    it.DVHType = "CUMULATIVE"
    it.DoseUnits = "GY"
    it.DoseType = "PHYSICAL"
    it.DVHDoseScaling = "1"
    it.DVHVolumeUnits = "CM3"
    it.DVHNumberOfBins = int(len(dvh["edges"]))
    width = _ds(DVH_BIN_GY)
    data = []
    for v in dvh["volumes"]:
        data.extend([width, _ds(v)])
    it.DVHData = data
    it.DVHMinimumDose = _ds(dvh["dmin"])
    it.DVHMaximumDose = _ds(dvh["dmax"])
    it.DVHMeanDose = _ds(dvh["dmean"])
    return it


def _rtdose(ctx: dict, sop: str, series: str, dose_fn, summation: str,
            dvh_items: Optional[list], beam_number: Optional[int] = None) -> Dataset:
    ds = _base(RTDOSE_SOP_CLASS, sop, "RTDOSE", ctx)
    ds.SeriesInstanceUID = series
    ds.SeriesNumber = 4
    ds.FrameOfReferenceUID = ctx["for_uid"]
    ds.PositionReferenceIndicator = ""
    xs = (np.arange(DOSE_NX) - (DOSE_NX - 1) / 2.0) * DOSE_RES
    ys = (np.arange(DOSE_NY) - (DOSE_NY - 1) / 2.0) * DOSE_RES
    zs = ctx["planes"]
    Z, Y, X = np.meshgrid(zs, ys, xs, indexing="ij")
    dose = dose_fn(np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])).reshape(Z.shape)
    px = np.round(dose / DOSE_SCALING).astype(np.uint32)
    ds.ImagePositionPatient = [_ds(xs[0]), _ds(ys[0]), _ds(zs[0])]
    ds.ImageOrientationPatient = ["1", "0", "0", "0", "1", "0"]
    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = "MONOCHROME2"
    ds.NumberOfFrames = int(len(zs))
    ds.FrameIncrementPointer = Tag(0x3004, 0x000C)
    ds.Rows, ds.Columns = DOSE_NY, DOSE_NX
    ds.PixelSpacing = [_ds(DOSE_RES), _ds(DOSE_RES)]
    ds.BitsAllocated, ds.BitsStored, ds.HighBit, ds.PixelRepresentation = 32, 32, 31, 0
    ds.DoseUnits, ds.DoseType, ds.DoseSummationType = "GY", "PHYSICAL", summation
    ds.GridFrameOffsetVector = [_ds(z - zs[0]) for z in zs]
    ds.DoseGridScaling = _ds(DOSE_SCALING)

    rp_ref = Dataset()
    rp_ref.ReferencedSOPClassUID = RTPLAN_SOP_CLASS
    rp_ref.ReferencedSOPInstanceUID = ctx["rp_sop"]
    if beam_number is not None:
        rb = Dataset()
        rb.ReferencedBeamNumber = beam_number
        rfg = Dataset()
        rfg.ReferencedFractionGroupNumber = 1
        rfg.ReferencedBeamSequence = Sequence([rb])
        rp_ref.ReferencedFractionGroupSequence = Sequence([rfg])
    ds.ReferencedRTPlanSequence = Sequence([rp_ref])
    rs_ref = Dataset()
    rs_ref.ReferencedSOPClassUID = RTSTRUCT_SOP_CLASS
    rs_ref.ReferencedSOPInstanceUID = ctx["rs_sop"]
    ds.ReferencedStructureSetSequence = Sequence([rs_ref])
    if dvh_items:
        ds.DVHSequence = Sequence(dvh_items)
    ds.PixelData = px.tobytes()
    return ds


# ---------------------------------------------------------------------------
# Sollwerte
# ---------------------------------------------------------------------------

def _expected(rois: list, planes: np.ndarray, fields: tuple) -> dict:
    f1, f2 = fields
    t1 = ph.SphereTarget("PTV_1", np.asarray(PTV1_C), PTV1_R, RX1, f1,
                         np.array([z for z in planes if abs(z - PTV1_C[2]) < PTV1_R]), SLICE_DZ, N_POLY)
    t2 = ph.SphereTarget("PTV_2", np.asarray(PTV2_C), PTV2_R, RX2, f2,
                         np.array([z for z in planes if abs(z - PTV2_C[2]) < PTV2_R]), SLICE_DZ, N_POLY)
    targets = {}
    for t, other in ((t1, (f2,)), (t2, (f1,))):
        targets[t.name] = {
            "slab": t.expected(planes, other, "slab"),
            "eclipse": t.expected(planes, other, "eclipse"),
        }
    # PTV_2 mit der Verschreibung von PTV_1 (Standardlauf wertet nur PTV_1 aus;
    # mit --target PTV_1,PTV_2 --rx 20 gilt fuer beide Rx = 20 Gy)
    targets["PTV_2"]["slab_rx20"] = t2.expected(planes, (f1,), "slab", rx=RX1)

    volumes = {}
    for roi in rois:
        if roi.point or not roi.contours:
            volumes[roi.name] = 0.0
            continue
        volumes[roi.name] = float(ana.compute_volume(roi.contours))
    hs = HIRNSTAMM
    ptv1_hs = math.hypot(PTV1_C[0] - hs["cx"], PTV1_C[1] - hs["cy"]) - PTV1_R - hs["r"]
    return {
        "targets": targets,
        "structure_volumes_cm3": volumes,
        "markers": {"HS1": list(HS1_POS), "Iso": list(ISO_POS)},
        "surface_distance_mm": {"PTV_1-Hirnstamm": ptv1_hs},
        "fields": [{"center": f.center.tolist(), "dmax_gy": f.dmax, "r0_mm": f.r0, "p": f.p}
                   for f in fields],
        "notes": [
            "Zielwerte in geschlossener Form (Scheibenstapel); Uebersprechen des anderen "
            "Feldes vernachlaessigt (< 0.02 Gy).",
            "Strukturvolumina planimetrisch wie analyzer.compute_volume (XOR x nominale Schichtdicke).",
            "PIV/PIV50 'component' = nur das eigene Feld, 'global' = beide Felder.",
        ],
    }


def _eclipse_ref_payload(expected: dict, which: str) -> dict:
    t1 = expected["targets"]["PTV_1"]["slab"]
    out = {"_meta": {"source": "Demo (synthetisch)", "date": "2026-01-01"},
           "PTV_1": {"PIV": round(t1["piv_global_cm3"], 4), "V10Gy": round(t1["piv50_global_cm3"], 4)}}
    if which == "both":
        t2 = expected["targets"]["PTV_2"]["slab"]
        out["PTV_2"] = {"TV": round(t2["tv_cm3"], 4)}
    return out


# ---------------------------------------------------------------------------
# Hauptfunktion
# ---------------------------------------------------------------------------

def make_demo_case(out_dir, spec: Optional[DemoSpec] = None) -> DemoCase:
    """Schreibt den Demo-Fall nach ``out_dir`` (wird angelegt) und liefert die Pfade."""
    spec = spec or DemoSpec()
    spec.validate()
    root = Path(out_dir)
    root.mkdir(parents=True, exist_ok=True)
    ct_dir = root / "CT" if spec.layout == "ct_subdir" else root
    uid = _Uids(spec.seed)
    planes = _planes()
    ctx = {
        "uid": uid, "planes": planes,
        "ts": ExplicitVRLittleEndian if spec.explicit_vr else ImplicitVRLittleEndian,
        "study_uid": uid("study"), "for_uid": uid("for"),
        "ct_series": uid("ct_series"), "rs_series": uid("rs_series"),
        "rp_series": uid("rp_series"), "rs_sop": uid("rs"), "rp_sop": uid("rp"),
    }
    ctx["ct_sops"] = [uid("ct", k) for k in range(len(planes))]

    ct_files = []
    for ds in _ct_slices(ctx, spec):
        ct_files.append(_write(ds, ct_dir / f"CT.{ds.SOPInstanceUID}.dcm"))

    rois = _build_rois(planes, spec.empty_roi)
    rs_path = _write(_rtstruct(ctx, rois), root / f"RS.{ctx['rs_sop']}.dcm")
    rp_path = _write(_rtplan(ctx, with_beams=(spec.rd_set == "plan+2beams")),
                     root / f"RP.{ctx['rp_sop']}.dcm")

    fields = _fields()
    dose_fn = ph.summed_dose(fields)
    dvh_items = None
    if spec.dvh:
        by_num = {r.number: r for r in rois}
        hx, hy = DOSE_NX / 2.0 * DOSE_RES, DOSE_NY / 2.0 * DOSE_RES
        box = (-hx, hx, -hy, hy)                 # Dosisgitter; ausserhalb Dosis 0
        dvh_items = []
        for num, res in ((2, 0.25), (4, 0.25), (1, 0.5), (5, 0.5)):
            cont = by_num[num].contours
            dvh = ph.cumulative_dvh(cont, dose_fn, res, DVH_BIN_GY, SLICE_DZ, box=box,
                                    total_volume=float(ana.compute_volume(cont)))
            dvh_items.append(_dvh_item(num, dvh))

    rd_paths = []
    rd_plan_sop = uid("rd", "plan")
    rd_paths.append(_write(_rtdose(ctx, rd_plan_sop, uid("rd_series", "plan"), dose_fn, "PLAN", dvh_items),
                           root / f"RD.{rd_plan_sop}.dcm"))
    if spec.rd_set == "plan+2beams":
        for b in (1, 2):
            sop = uid("rd", "beam", b)
            half = (lambda p, _f=dose_fn: 0.5 * _f(p))
            rd_paths.append(_write(_rtdose(ctx, sop, uid("rd_series", "beam", b), half, "BEAM", None,
                                           beam_number=b), root / f"RD.{sop}.dcm"))
    elif spec.rd_set == "2plans":
        sop = uid("rd", "plan2")
        scaled = (lambda p, _f=dose_fn: 0.9 * _f(p))
        rd_paths.append(_write(_rtdose(ctx, sop, uid("rd_series", "plan2"), scaled, "PLAN", None),
                               root / f"RD.{sop}.dcm"))

    expected = _expected(rois, planes, fields)
    expected["_meta"] = {"generator": "dicom_file_modifier.demo", "demo_version": DEMO_VERSION,
                         "spec": asdict(spec), "patient": "DFM^DEMO (synthetisch)"}
    exp_path = root / "demo_expected.json"
    exp_path.write_text(json.dumps(expected, indent=2, ensure_ascii=False), encoding="utf-8")

    ecl_path = None
    if spec.eclipse_ref:
        ecl_path = root / "eclipse_ref.json"
        ecl_path.write_text(json.dumps(_eclipse_ref_payload(expected, spec.eclipse_ref), indent=2),
                            encoding="utf-8")

    return DemoCase(root=root, ct_dir=ct_dir, ct_files=ct_files, rs=rs_path, rp=rp_path, rd=rd_paths,
                    eclipse_ref=ecl_path, expected_json=exp_path, spec=spec,
                    uids={k: v for k, v in ctx.items() if k.endswith(("_uid", "_series", "_sop"))})


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Synthetischen Demo-Fall (CT + RTSTRUCT + RTPLAN + RTDOSE, kein Patient) erzeugen.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Beispiel:\n  python -m dicom_file_modifier.demo output/demo --layout flat",
    )
    p.add_argument("out_dir", help="Ausgabeordner (wird angelegt)")
    p.add_argument("--layout", choices=("ct_subdir", "flat"), default="ct_subdir",
                   help="ct_subdir: CT/ als Unterordner (Default); flat: alles in einem Ordner wie ein Eclipse-Export")
    p.add_argument("--rd-set", choices=("plan", "plan+2beams", "2plans"), default="plan",
                   help="RTDOSE-Dateien: nur Plansumme (Default), zusaetzlich 2 Felddosen, oder 2 Plansummen")
    p.add_argument("--no-dvh", action="store_true", help="RTDOSE ohne DVHSequence")
    p.add_argument("--explicit-vr", action="store_true", help="Explicit statt Implicit VR Little Endian")
    p.add_argument("--rle", action="store_true", help="CT-Schichten RLE-komprimiert")
    p.add_argument("--eclipse-ref", choices=("ptv1", "both"), default=None,
                   help="eclipse_ref.json mit PTV_1 (ptv1) oder PTV_1 + PTV_2 (both) anlegen")
    p.add_argument("--no-empty-roi", action="store_true", help="Keine leere ROI anlegen")
    p.add_argument("--size", choices=tuple(SIZES), default="small",
                   help="small: 128x128 bei 2 mm (Default); medium: 256x256 bei 1 mm")
    p.add_argument("--seed", default="dfm-demo", help="Seed fuer deterministische UIDs (Default dfm-demo)")
    p.add_argument("--random-uids", action="store_true", help="Zufaellige statt deterministischer UIDs")
    return p


def main(argv: Optional[list] = None) -> int:
    args = _build_parser().parse_args(argv)
    spec = DemoSpec(layout=args.layout, rd_set=args.rd_set, dvh=not args.no_dvh,
                    explicit_vr=args.explicit_vr, compress_ct="rle" if args.rle else None,
                    eclipse_ref=args.eclipse_ref, empty_roi=not args.no_empty_roi, size=args.size,
                    seed=None if args.random_uids else args.seed)
    try:
        case = make_demo_case(args.out_dir, spec)
    except (ValueError, OSError) as e:
        print(f"Fehler: {e}", file=sys.stderr)
        return 2
    print(f"Demo-Fall geschrieben: {case.root}")
    print(f"  CT : {len(case.ct_files)} Schichten in {case.ct_dir}")
    print(f"  RS : {case.rs.name}")
    print(f"  RP : {case.rp.name}")
    for rd in case.rd:
        print(f"  RD : {rd.name}")
    if case.eclipse_ref:
        print(f"  Eclipse-Referenz: {case.eclipse_ref.name}")
    print(f"  Sollwerte: {case.expected_json.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
