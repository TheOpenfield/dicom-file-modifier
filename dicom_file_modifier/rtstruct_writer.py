"""
rtstruct_writer.py - Separate RTSTRUCT mit Isodosen-ROIs und Hilfskonturen
(Schnitt, Unterdosierung, Spill) fuer die Dosisindex-Berechnung.

Die Datei wird NEU aufgebaut (kein Deepcopy des Original-RS): Patient/Studie/
Frame-of-Reference werden aus dem Original uebernommen, damit sie in Eclipse/
ARIA neben der Original-CT-Serie landet; alle ROIs, UIDs und Labels sind neu.
Jede Kontur referenziert die CT-Schicht ihrer z-Position (ContourImageSequence).

Kein eigenes CLI ausser ``--self-test``:
  python -m dicom_file_modifier.rtstruct_writer --self-test
"""

from __future__ import annotations

import copy
import datetime as _dt
import sys
import threading
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pydicom
from pydicom import config as pdconfig
from pydicom.dataset import Dataset, FileMetaDataset, validate_file_meta
from pydicom.sequence import Sequence
from pydicom.uid import (PYDICOM_IMPLEMENTATION_UID, ExplicitVRLittleEndian,
                         ImplicitVRLittleEndian, generate_uid)

from . import dose as dm
from . import analyzer as ana
from .dicom_utils import _label_with_suffix, _truncate, get_rs_frame_of_references, set_sop_instance_uid
from .dose_constants import HELPER_COLORS, TOOL_NAME, TOOL_VERSION

RTSTRUCT_SOP_CLASS = "1.2.840.10008.5.1.4.1.1.481.3"
CT_IMAGE_SOP_CLASS = "1.2.840.10008.5.1.4.1.1.2"

# Aus dem Original uebernommene Top-Level-Attribute (Patient/Studie/Kontext)
_COPY_KEYWORDS = (
    "SpecificCharacterSet", "StudyDate", "StudyTime", "AccessionNumber", "StudyDescription",
    "StudyID", "StudyInstanceUID", "ReferringPhysicianName", "FrameOfReferenceUID",
    "PositionReferenceIndicator", "InstitutionName", "InstitutionAddress",
    "InstitutionalDepartmentName", "PatientName", "PatientID", "PatientBirthDate",
    "PatientSex", "PatientAge", "PatientWeight", "OtherPatientIDs",
    "CodingSchemeIdentificationSequence", "ContextGroupIdentificationSequence",
    "MappingResourceIdentificationSequence",
)


@dataclass
class RoiSpec:
    """Beschreibung einer zu schreibenden ROI."""
    name: str
    color: tuple
    contours_by_z: Optional[dict] = None        # {z: [(N,3), ...]}
    copy_from_roi: Optional[int] = None         # Zielkopie: ContourSequence verbatim
    interpreted_type: str = "CONTROL"
    generation_algorithm: str = "AUTOMATIC"
    description: str = ""
    identification_code: Optional[Dataset] = None
    kind: str = "isodose"                       # isodose | inter | under | spill | target


# ---------------------------------------------------------------------------
# 1. CT-Schichtindex
# ---------------------------------------------------------------------------

def build_ct_slice_index(ct_dir) -> dict:
    """
    Liest die CT-Header (ohne Pixel) und liefert
    ``{'series_uid', 'study_uid', 'for_uid', 'sop_class_uid', 'z_values',
    'z_to_sop', 'sops', 'pixel_spacing', 'ipp_xy', 'n_slices', 'z_to_path',
    'paths'}`` (``z_to_path``/``paths``: Dateipfade je Schicht, damit
    ``dose_viz`` nur die Schichten im Bereich des Feingitters mit Pixeln laedt).
    ``ct_dir`` ist ein Ordner (``*.dcm``) oder eine Liste von Dateien.
    """
    import glob
    import os

    if isinstance(ct_dir, (str, os.PathLike)):
        files = sorted(glob.glob(os.path.join(glob.escape(str(ct_dir)), "*.dcm")))
        where = repr(ct_dir)
    else:
        files = [str(f) for f in ct_dir]
        where = f"den {len(files)} CT-Dateien"
    slices = []
    for f in files:
        try:
            ds = pydicom.dcmread(f, stop_before_pixels=True)
        except Exception:
            continue
        if str(ds.get("Modality", "")).upper() == "CT" and "ImagePositionPatient" in ds:
            slices.append(ds)
    if not slices:
        raise ValueError(f"Keine CT-Schichten in {where} gefunden.")
    slices.sort(key=lambda s: float(s.ImagePositionPatient[2]))
    series = {str(s.SeriesInstanceUID) for s in slices}
    fors = {str(s.get("FrameOfReferenceUID", "")) for s in slices}
    if len(series) != 1:
        raise ValueError(f"CT-Ordner enthaelt {len(series)} Serien; genau eine erwartet.")
    if len(fors) != 1:
        raise ValueError("CT-Schichten haben unterschiedliche FrameOfReferenceUIDs.")
    z_values = np.array([float(s.ImagePositionPatient[2]) for s in slices])
    z_to_sop = {round(float(z), 3): str(s.SOPInstanceUID) for z, s in zip(z_values, slices)}
    z_to_path = {round(float(z), 3): str(getattr(s, "filename", "") or "")
                 for z, s in zip(z_values, slices)}
    s0 = slices[0]
    return {
        "series_uid": series.pop(),
        "study_uid": str(s0.get("StudyInstanceUID", "")),
        "for_uid": fors.pop(),
        "sop_class_uid": str(s0.get("SOPClassUID", CT_IMAGE_SOP_CLASS)),
        "z_values": z_values,
        "z_to_sop": z_to_sop,
        "sops": set(z_to_sop.values()),
        "pixel_spacing": (float(s0.PixelSpacing[0]), float(s0.PixelSpacing[1])),
        "ipp_xy": (float(s0.ImagePositionPatient[0]), float(s0.ImagePositionPatient[1])),
        "n_slices": len(slices),
        "z_to_path": z_to_path,
        "paths": [z_to_path[round(float(z), 3)] for z in z_values],
    }


def validate_index_against_rs(ct_index: dict, orig_rs: pydicom.Dataset) -> list:
    """
    FoR- oder Serien-Mismatch zwischen CT-Ordner und Original-RS -> ``ValueError``;
    fehlende referenzierte Schichten -> Warnungstexte.
    """
    warnings = []
    rs_fors = get_rs_frame_of_references(orig_rs)
    if rs_fors and ct_index["for_uid"] not in rs_fors:
        raise ValueError(
            f"FrameOfReferenceUID des CT ({ct_index['for_uid']}) kommt im RTSTRUCT nicht vor."
        )
    ref_series, missing = set(), 0
    for rf in orig_rs.get("ReferencedFrameOfReferenceSequence", []):
        for st in rf.get("RTReferencedStudySequence", []):
            for se in st.get("RTReferencedSeriesSequence", []):
                ref_series.add(str(se.SeriesInstanceUID))
                for ci in se.get("ContourImageSequence", []):
                    if str(ci.ReferencedSOPInstanceUID) not in ct_index["sops"]:
                        missing += 1
    if ref_series and ct_index["series_uid"] not in ref_series:
        raise ValueError(
            f"RTSTRUCT referenziert die CT-Serie(n) {sorted(ref_series)}, der CT-Ordner "
            f"enthaelt {ct_index['series_uid']}."
        )
    if missing:
        warnings.append(f"{missing} im RTSTRUCT referenzierte CT-Schichten fehlen im CT-Ordner.")
    return warnings


def _sop_for_z(ct_index: dict, z: float, z_tol_mm: float) -> str:
    sop = ct_index["z_to_sop"].get(round(float(z), 3))
    if sop is not None:
        return sop
    zv = ct_index["z_values"]
    k = int(np.argmin(np.abs(zv - z)))
    if abs(zv[k] - z) <= z_tol_mm:
        return ct_index["z_to_sop"][round(float(zv[k]), 3)]
    raise ValueError(f"Keine CT-Schicht bei z={z:.3f} mm (naechste: {zv[k]:.3f} mm).")


# ---------------------------------------------------------------------------
# 2. ROI-Spezifikationen aus den Artefakten (Namen, Farben)
# ---------------------------------------------------------------------------

def iso_name(pct: float, gy: float) -> str:
    return f"ISO_{pct:g}%_{gy:.1f}Gy"


def _fit_name(prefix: str, target: str, suffix: str, max_len: int) -> str:
    """Kuerzt nur den Zielanteil, damit ``prefix + target + suffix`` in ``max_len`` passt."""
    room = max_len - len(prefix) - len(suffix)
    if room < 1:
        return (prefix + suffix)[:max_len]
    return prefix + target[:room] + suffix


def _unique_names(specs: list) -> None:
    """Haengt ``_2``, ``_3`` ... an case-insensitive Namenskollisionen an (in place)."""
    seen = {}
    for sp in specs:
        base = sp.name
        key = base.lower()
        n = seen.get(key, 0)
        if n:
            cand = f"{base}_{n + 1}"
            while cand.lower() in seen:
                n += 1
                cand = f"{base}_{n + 1}"
            sp.name = cand
            seen[cand.lower()] = 1
        seen[key] = n + 1


def _dose_code_item(orig_rs: pydicom.Dataset) -> Optional[Dataset]:
    """RTROIIdentificationCodeSequence-Item einer Eclipse-Dosis-ROI (CodeValue 'Dose')."""
    for obs in orig_rs.get("RTROIObservationsSequence", []):
        for item in obs.get("RTROIIdentificationCodeSequence", []):
            if str(item.get("CodeValue", "")).lower() == "dose":
                return copy.deepcopy(item)
    return None


def build_roi_specs(art, include_target: bool = False, max_name_len: int = 64,
                    simplify_mm: float = 0.1, vertex_mode: str = "edge",
                    iso_contours: str = "mask") -> list:
    """
    ``[RoiSpec, ...]`` aus ``DoseIndexArtifacts``: Zielkopien (optional),
    Isodosen absteigend, je Ziel Schnitt / Unterdosierung / Spill (Rx-Level).
    Setzt ``art.levels[*].roi_name`` und ``art.targets[*].helper_names``.
    """
    grid = art.grid
    specs = []
    code = _dose_code_item(art.rs_ds) if art.rs_ds is not None else None
    types = ana.get_structure_type(art.rs_ds) if art.rs_ds is not None else {}

    if include_target:
        for tm in art.targets.values():
            specs.append(RoiSpec(
                name=tm.name[:max_name_len], color=tm.color, copy_from_roi=tm.roi_number,
                interpreted_type=types.get(tm.roi_number, "PTV") or "PTV",
                generation_algorithm="MANUAL", description="Zielkopie (Original-Konturen)",
                kind="target",
            ))
    for lv in art.levels.values():
        lv.roi_name = iso_name(lv.pct, lv.gy)[:max_name_len]
        if iso_contours == "field":
            iso_cz = dm.field_to_contours(art.dose_fine, lv.gy, grid, simplify_mm=simplify_mm)
        else:
            iso_cz = dm.mask_to_contours(lv.mask, grid, simplify_mm=simplify_mm,
                                         vertex_mode=vertex_mode)
        specs.append(RoiSpec(
            name=lv.roi_name, color=lv.color, contours_by_z=iso_cz,
            interpreted_type="CONTROL", description=f"Isodose {lv.label} = {lv.gy:.2f} Gy "
            f"(Rx {art.rx_gy:.2f} Gy), {art.settings.get('grid_mm', grid.res_xy):g} mm Raster",
            identification_code=code, kind="isodose",
        ))
    for tm in art.targets.values():
        names = {
            "inter": _fit_name("", tm.name, "_x_ISO100", max_name_len),
            "under": _fit_name("", tm.name, "_minus_ISO100", max_name_len),
            "spill": _fit_name("ISO100_minus_", tm.name, "", max_name_len),
        }
        masks = {"inter": tm.inter, "under": tm.under, "spill": tm.spill}
        descs = {"inter": "Schnitt Ziel & Rx-Isodose (TV&PIV)",
                 "under": "Ziel ausserhalb der Rx-Isodose (unterdosiert)",
                 "spill": "Rx-Isodose ausserhalb des Ziels (Spill)"}
        for kind in ("inter", "under", "spill"):
            specs.append(RoiSpec(
                name=names[kind], color=HELPER_COLORS[kind],
                contours_by_z=dm.mask_to_contours(masks[kind], grid, simplify_mm=simplify_mm,
                                                  vertex_mode=vertex_mode),
                interpreted_type="CONTROL", description=f"{descs[kind]}: {tm.name}",
                kind=kind,
            ))
        tm.helper_names = {"intersection": names["inter"], "underdosed": names["under"],
                           "spill": names["spill"]}
    _unique_names(specs)
    # Namen nach der Eindeutigkeits-Passe zurueckschreiben
    by_kind = {(sp.kind, sp.description): sp.name for sp in specs}
    for tm in art.targets.values():
        for kind, key in (("inter", "intersection"), ("under", "underdosed"), ("spill", "spill")):
            for sp in specs:
                if sp.kind == kind and sp.description.endswith(f": {tm.name}"):
                    tm.helper_names[key] = sp.name
    for lv in art.levels.values():
        for sp in specs:
            if sp.kind == "isodose" and sp.description.startswith(f"Isodose {lv.label} "):
                lv.roi_name = sp.name
    return specs


def planned_roi_names(level_specs: list, target_names: list, include_target: bool = False,
                      max_name_len: int = 64) -> list:
    """
    ROI-Namen, die ``build_roi_specs`` schreiben wuerde, ohne Rechnung (Vorschau):
    ``[(kind, name, bezug), ...]`` in derselben Reihenfolge und mit derselben
    Eindeutigkeits-Passe; ``bezug`` ist der Level-``label`` bzw. der Zielname.
    """
    stubs, refs = [], []
    if include_target:
        for t in target_names:
            stubs.append(RoiSpec(name=t[:max_name_len], color=(0, 0, 0), kind="target"))
            refs.append(t)
    for lv in level_specs:
        stubs.append(RoiSpec(name=iso_name(lv["pct"], lv["gy"])[:max_name_len], color=(0, 0, 0),
                             kind="isodose"))
        refs.append(lv["label"])
    for t in target_names:
        for kind, prefix, suffix in (("inter", "", "_x_ISO100"), ("under", "", "_minus_ISO100"),
                                     ("spill", "ISO100_minus_", "")):
            stubs.append(RoiSpec(name=_fit_name(prefix, t, suffix, max_name_len), color=(0, 0, 0),
                                 kind=kind))
            refs.append(t)
    _unique_names(stubs)
    return [(sp.kind, sp.name, ref) for sp, ref in zip(stubs, refs)]


def summary_description(art, max_len: int = 200) -> str:
    """Kurzfassung der Indizes fuer ``StructureSetDescription`` (ST)."""
    parts = []
    for name, tm in art.targets.items():
        ix = tm.result["indices"]
        parts.append(f"{name}: CI {ix['ci_paddick']:.3f} GI {ix['gi']:.2f} HI {ix['hi_icru83']:.3f}")
    s = " | ".join(parts) + f" | Rx {art.rx_gy:.2f} Gy | dicom_file_modifier dose_indices"
    return _truncate(s, max_len)


# ---------------------------------------------------------------------------
# 3. RTSTRUCT schreiben
# ---------------------------------------------------------------------------

_WRITING_VALIDATION_LOCK = threading.RLock()


@contextmanager
def _writing_validation(mode):
    """
    Setzt ``writing_validation_mode`` fuer die Dauer des Blocks.  Der Schalter
    gilt prozessweit; das Lock verhindert, dass zwei Threads ihn gleichzeitig
    umstellen und sich den alten Wert gegenseitig ueberschreiben.
    """
    with _WRITING_VALIDATION_LOCK:
        old = pdconfig.settings.writing_validation_mode
        pdconfig.settings.writing_validation_mode = mode
        try:
            yield
        finally:
            pdconfig.settings.writing_validation_mode = old


def _now_strings() -> tuple:
    now = _dt.datetime.now()
    return now.strftime("%Y%m%d"), now.strftime("%H%M%S")


def _tool_version() -> str:
    # SoftwareVersions hat VR LO (64 Zeichen); bis P0.4 auf 16 gekuerzt
    return f"{TOOL_NAME} {TOOL_VERSION}"[:64]


def write_isodose_rtstruct(orig_rs: pydicom.Dataset, ct_index: dict, rois: list,
                           out_path, label: str = "_IDX", description: str = "",
                           series_number_offset: int = 1000, transfer_syntax: str = "explicit",
                           z_tol_mm: float = 0.01, coord_decimals: int = 3,
                           max_name_len: int = 64) -> pydicom.Dataset:
    """
    Schreibt eine neue RTSTRUCT-Datei mit den ``rois`` (``RoiSpec``) und gibt das
    Dataset zurueck.  Patient/Studie/FoR aus ``orig_rs``; Konturen referenzieren
    die CT-Schicht ihrer z-Position aus ``ct_index``.
    """
    out_path = Path(out_path)
    ds = Dataset()
    for kw in _COPY_KEYWORDS:
        if kw in orig_rs:
            ds[kw] = copy.deepcopy(orig_rs[kw])

    date, time = _now_strings()
    ds.SOPClassUID = RTSTRUCT_SOP_CLASS
    ds.Modality = "RTSTRUCT"
    ds.Manufacturer = "dicom_file_modifier"
    ds.ManufacturerModelName = "dose_indices"
    ds.SoftwareVersions = _tool_version()
    ds.InstanceCreationDate, ds.InstanceCreationTime = date, time
    ds.SeriesDate, ds.SeriesTime = date, time
    ds.StructureSetDate, ds.StructureSetTime = date, time
    ds.SeriesInstanceUID = generate_uid()
    ds.SeriesDescription = _truncate(f"Dosisindizes{label}", 64)
    try:
        ds.SeriesNumber = int(orig_rs.get("SeriesNumber", 0) or 0) + series_number_offset
    except (TypeError, ValueError):
        ds.SeriesNumber = series_number_offset
    ds.StructureSetLabel = _label_with_suffix(str(orig_rs.get("StructureSetLabel", "RS")), label, 16)
    ds.StructureSetName = _truncate(f"Isodosen/Indexkonturen{label}", 64)
    ds.StructureSetDescription = _truncate(description or "dicom_file_modifier dose_indices", 1024)
    ds.ApprovalStatus = "UNAPPROVED"

    # Frame-of-Reference-Referenz auf die CT-Serie (Kopie, sonst Neuaufbau)
    for_uid = ct_index["for_uid"]
    rfor_seq = orig_rs.get("ReferencedFrameOfReferenceSequence")
    if rfor_seq:
        ds.ReferencedFrameOfReferenceSequence = copy.deepcopy(rfor_seq)
    else:
        ci_seq = Sequence()
        for z in ct_index["z_values"]:
            ci = Dataset()
            ci.ReferencedSOPClassUID = ct_index["sop_class_uid"]
            ci.ReferencedSOPInstanceUID = ct_index["z_to_sop"][round(float(z), 3)]
            ci_seq.append(ci)
        se = Dataset()
        se.SeriesInstanceUID = ct_index["series_uid"]
        se.ContourImageSequence = ci_seq
        st = Dataset()
        st.ReferencedSOPClassUID = "1.2.840.10008.3.1.2.3.1"   # Detached Study Management
        st.ReferencedSOPInstanceUID = ct_index["study_uid"] or generate_uid()
        st.RTReferencedSeriesSequence = Sequence([se])
        rf = Dataset()
        rf.FrameOfReferenceUID = for_uid
        rf.RTReferencedStudySequence = Sequence([st])
        ds.ReferencedFrameOfReferenceSequence = Sequence([rf])
    ds.FrameOfReferenceUID = for_uid

    ssr_seq, rc_seq, obs_seq = Sequence(), Sequence(), Sequence()
    for num, spec in enumerate(rois, start=1):
        name = spec.name[:max_name_len]
        ssr = Dataset()
        ssr.ROINumber = num
        ssr.ReferencedFrameOfReferenceUID = for_uid
        ssr.ROIName = name
        ssr.ROIGenerationAlgorithm = spec.generation_algorithm
        ssr_seq.append(ssr)

        rc = Dataset()
        rc.ReferencedROINumber = num
        rc.ROIDisplayColor = [int(v) for v in spec.color]
        contours = Sequence()
        if spec.copy_from_roi is not None:
            src = [r for r in orig_rs.ROIContourSequence
                   if int(r.ReferencedROINumber) == int(spec.copy_from_roi)]
            if not src:
                raise ValueError(f"ROI {spec.copy_from_roi} fuer die Zielkopie nicht im Original.")
            for c in src[0].get("ContourSequence", []):
                cc = copy.deepcopy(c)
                for ci in cc.get("ContourImageSequence", []):
                    if str(ci.ReferencedSOPInstanceUID) not in ct_index["sops"]:
                        raise ValueError(
                            f"Zielkopie {name!r} referenziert eine CT-Schicht ausserhalb des CT-Ordners."
                        )
                contours.append(cc)
            rc.ROIDisplayColor = list(src[0].get("ROIDisplayColor", rc.ROIDisplayColor))
        else:
            for z in sorted(spec.contours_by_z or {}):
                sop = _sop_for_z(ct_index, z, z_tol_mm)
                for ring in spec.contours_by_z[z]:
                    pts = np.asarray(ring, dtype=float)
                    if len(pts) < 3:
                        continue
                    c = Dataset()
                    ci = Dataset()
                    ci.ReferencedSOPClassUID = ct_index["sop_class_uid"]
                    ci.ReferencedSOPInstanceUID = sop
                    c.ContourImageSequence = Sequence([ci])
                    c.ContourGeometricType = "CLOSED_PLANAR"
                    c.NumberOfContourPoints = int(len(pts))
                    c.ContourData = [f"{v:.{coord_decimals}f}" for v in pts.ravel()]
                    contours.append(c)
        rc.ContourSequence = contours
        rc_seq.append(rc)

        obs = Dataset()
        obs.ObservationNumber = num
        obs.ReferencedROINumber = num
        obs.ROIObservationLabel = _truncate(name, 16)
        if spec.description:
            obs.ROIObservationDescription = _truncate(spec.description, 1024)
        if spec.identification_code is not None:
            obs.RTROIIdentificationCodeSequence = Sequence([copy.deepcopy(spec.identification_code)])
        obs.RTROIInterpretedType = spec.interpreted_type
        obs.ROIInterpreter = ""
        obs_seq.append(obs)
    ds.StructureSetROISequence = ssr_seq
    ds.ROIContourSequence = rc_seq
    ds.RTROIObservationsSequence = obs_seq

    # File-Meta + Transfer-Syntax
    ts = ExplicitVRLittleEndian if transfer_syntax == "explicit" else ImplicitVRLittleEndian
    fm = FileMetaDataset()
    fm.MediaStorageSOPClassUID = RTSTRUCT_SOP_CLASS
    fm.MediaStorageSOPInstanceUID = generate_uid()
    fm.TransferSyntaxUID = ts
    fm.ImplementationClassUID = PYDICOM_IMPLEMENTATION_UID
    fm.ImplementationVersionName = f"PYDICOM {pydicom.__version__}"[:16]
    validate_file_meta(fm, enforce_standard=True)
    ds.file_meta = fm
    set_sop_instance_uid(ds, str(fm.MediaStorageSOPInstanceUID))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with _writing_validation(pdconfig.RAISE):
        pydicom.dcmwrite(str(out_path), ds, enforce_file_format=True)   # Kodierung aus der Transfer Syntax
    return ds


# ---------------------------------------------------------------------------
# 4. Pruefung der geschriebenen Datei
# ---------------------------------------------------------------------------

def verify_rtstruct(path, ct_index: dict, z_tol_mm: float = 0.01) -> list:
    """Liest die Datei zurueck und prueft Konsistenz; leere Liste = OK."""
    problems = []
    ds = pydicom.dcmread(str(path))
    fm = getattr(ds, "file_meta", None)
    if fm is None or str(fm.get("MediaStorageSOPInstanceUID", "")) != str(ds.SOPInstanceUID):
        problems.append("file_meta.MediaStorageSOPInstanceUID != SOPInstanceUID")
    try:
        validate_file_meta(fm, enforce_standard=True)
    except Exception as e:
        problems.append(f"file_meta ungueltig: {e}")
    if str(ds.get("Modality", "")) != "RTSTRUCT" or str(ds.SOPClassUID) != RTSTRUCT_SOP_CLASS:
        problems.append("Modality/SOPClassUID nicht RTSTRUCT")
    if str(ds.get("FrameOfReferenceUID", "")) != ct_index["for_uid"]:
        problems.append("FrameOfReferenceUID != CT-FoR")
    if len(str(ds.get("StructureSetLabel", ""))) > 16:
        problems.append("StructureSetLabel > 16 Zeichen")
    ref_series = {str(se.SeriesInstanceUID)
                  for rf in ds.get("ReferencedFrameOfReferenceSequence", [])
                  for st in rf.get("RTReferencedStudySequence", [])
                  for se in st.get("RTReferencedSeriesSequence", [])}
    if ct_index["series_uid"] not in ref_series:
        problems.append("RTReferencedSeriesSequence referenziert nicht die CT-Serie")

    roi_nums = [int(r.ROINumber) for r in ds.get("StructureSetROISequence", [])]
    if len(roi_nums) != len(set(roi_nums)):
        problems.append("ROINumber nicht eindeutig")
    names = [str(r.ROIName) for r in ds.get("StructureSetROISequence", [])]
    if len({n.lower() for n in names}) != len(names):
        problems.append("ROI-Namen nicht eindeutig (case-insensitiv)")
    if any(len(n) > 64 for n in names):
        problems.append("ROI-Name > 64 Zeichen")
    for seq_name in ("ROIContourSequence", "RTROIObservationsSequence"):
        for item in ds.get(seq_name, []):
            if int(item.ReferencedROINumber) not in roi_nums:
                problems.append(f"{seq_name}: ReferencedROINumber {item.ReferencedROINumber} unbekannt")
    zv = ct_index["z_values"]
    sop_to_z = {v: k for k, v in ct_index["z_to_sop"].items()}
    n_contours = 0
    for rc in ds.get("ROIContourSequence", []):
        for c in rc.get("ContourSequence", []):
            n_contours += 1
            data = list(c.ContourData)
            npts = int(c.NumberOfContourPoints)
            if npts * 3 != len(data):
                problems.append(f"Kontur mit NumberOfContourPoints*3 != len(ContourData) (ROI {rc.ReferencedROINumber})")
            if any(len(str(v)) > 16 for v in data):
                problems.append("ContourData-Wert laenger als 16 Zeichen (DS)")
            pts = np.asarray(data, dtype=float).reshape(-1, 3)
            if np.ptp(pts[:, 2]) > 1e-6:
                problems.append("Kontur nicht planar (mehrere z)")
            cis = c.get("ContourImageSequence", [])
            if not cis:
                problems.append("Kontur ohne ContourImageSequence")
                continue
            sop = str(cis[0].ReferencedSOPInstanceUID)
            if sop not in sop_to_z:
                problems.append("Kontur referenziert eine SOP ausserhalb des CT-Ordners")
            elif abs(sop_to_z[sop] - float(pts[0, 2])) > z_tol_mm:
                problems.append(f"Kontur-z {pts[0, 2]:.3f} != z der referenzierten CT-Schicht")
    if n_contours == 0:
        problems.append("Keine Konturen geschrieben")
    # Duplikate zusammenfassen
    seen, uniq = set(), []
    for p in problems:
        if p not in seen:
            uniq.append(p)
            seen.add(p)
    return uniq


# ---------------------------------------------------------------------------
# 5. Self-Test (synthetische Masken, kein CT noetig)
# ---------------------------------------------------------------------------

def _run_self_test() -> int:
    import tempfile

    results = []

    def check(name, ok, detail=""):
        results.append(bool(ok))
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))

    print("\nRTSTRUCT-Writer Self-Test (synthetische Masken)")
    print("-" * 70)
    res = 0.25
    gx = -10 + (np.arange(80) + 0.5) * res
    gy = -10 + (np.arange(80) + 0.5) * res
    gz = np.arange(-4.5, 5.0, 1.0)
    grid = dm.FineGrid(gx=gx, gy=gy, gz=gz, res_xy=res, dz=1.0)
    Z, Y, X = np.meshgrid(gz, gy, gx, indexing="ij")

    def signed_area(ring):
        x, y = ring[:, 0], ring[:, 1]
        return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))

    # 1) Kugel
    sphere = (X ** 2 + Y ** 2 + Z ** 2) <= 6.0 ** 2
    d = dm.contours_roundtrip_dice(sphere, grid)
    check("Kugel: Roundtrip-Dice >= 0.99", d >= 0.99, f"Dice={d:.4f}")
    # 2) Annulus: zwei Ringe je Ebene, Aussen +, Loch -, XOR-Volumen = aussen - innen
    ann = ((X ** 2 + Y ** 2) <= 7.0 ** 2) & ((X ** 2 + Y ** 2) > 3.0 ** 2)
    cz = dm.mask_to_contours(ann, grid)
    rings = cz[float(gz[0])]
    areas = sorted(signed_area(r) for r in rings)
    vol_back = dm.mask_volume_cm3(ana.rasterize_contours(dm.contours_flat(cz), gx, gy, gz), grid)
    check("Annulus: 2 Ringe je Ebene, Aussenring +, Loch -", len(rings) == 2 and areas[0] < 0 < areas[1],
          f"Flaechen {np.round(areas, 1).tolist()}")
    check("Annulus: XOR-Volumen der Konturen = Maskenvolumen (0.5 %)",
          abs(vol_back - dm.mask_volume_cm3(ann, grid)) < 0.005 * dm.mask_volume_cm3(ann, grid),
          f"{vol_back:.4f} vs {dm.mask_volume_cm3(ann, grid):.4f}")
    # 3) Zwei Inseln
    two = (((X - 4) ** 2 + Y ** 2) <= 2.5 ** 2) | (((X + 4) ** 2 + Y ** 2) <= 2.5 ** 2)
    cz2 = dm.mask_to_contours(two, grid)
    check("Zwei Inseln: 2 Ringe je Ebene, Roundtrip-Dice >= 0.99",
          all(len(v) == 2 for v in cz2.values()) and dm.contours_roundtrip_dice(two, grid) >= 0.99)
    # 4) Randberuehrender Block schliesst
    edge = np.zeros_like(sphere)
    edge[:, :20, :20] = True
    cz3 = dm.mask_to_contours(edge, grid)
    a = abs(signed_area(cz3[float(gz[0])][0]))
    check("Randberuehrender Block: geschlossener Ring, Flaeche = 20x20 Voxel",
          abs(a - (20 * res) ** 2) < 0.01 * (20 * res) ** 2, f"A={a:.3f} mm2")
    # 5) Einzelvoxel: Raute um das Voxelzentrum (keine Simplifizierung)
    single = np.zeros_like(sphere)
    single[0, 40, 40] = True
    cz4 = dm.mask_to_contours(single, grid, simplify_mm=0.0, min_area_mm2=0.0)
    ring = cz4.get(float(gz[0]), [np.zeros((0, 3))])[0]
    ok5 = (len(ring) == 4 and abs(ring[:, 0].mean() - gx[40]) < 1e-9
           and abs(ring[:, 1].mean() - gy[40]) < 1e-9)
    check("Einzelvoxel: Raute (4 Punkte) zentriert auf dem Voxelmittelpunkt", ok5,
          f"n={len(ring)}, mitte=({ring[:, 0].mean():.4f},{ring[:, 1].mean():.4f}) vs ({gx[40]:.4f},{gy[40]:.4f})")
    # 6) Vertex-Modus 'center': Vertices auf Zellmittelpunkten
    cz5 = dm.mask_to_contours(sphere, grid, simplify_mm=0.0, vertex_mode="center")
    ring = cz5[float(gz[4])][0]
    on_grid = np.all(np.abs((ring[:, 0] - gx[0]) / res - np.round((ring[:, 0] - gx[0]) / res)) < 1e-6)
    check("Vertex-Modus center: alle Vertices auf Zellmittelpunkten", bool(on_grid), f"n={len(ring)}")

    # 7) Schreiben + Pruefen mit synthetischem Original-RS und CT-Index
    orig = Dataset()
    orig.SpecificCharacterSet = "ISO_IR 192"
    orig.PatientName, orig.PatientID = "Test^Phantom", "PH-001"
    orig.StudyInstanceUID = generate_uid()
    orig.FrameOfReferenceUID = "1.2.826.0.1.3680043.8.498.1.2.3"
    orig.StructureSetLabel = "Phantom_Original"
    orig.SeriesNumber = 7
    orig.StructureSetROISequence = Sequence()
    orig.ROIContourSequence = Sequence()
    orig.RTROIObservationsSequence = Sequence()
    sops = {round(float(z), 3): generate_uid() for z in gz}
    ct_index = {"series_uid": generate_uid(), "study_uid": orig.StudyInstanceUID,
                "for_uid": orig.FrameOfReferenceUID, "sop_class_uid": CT_IMAGE_SOP_CLASS,
                "z_values": np.array(gz), "z_to_sop": sops, "sops": set(sops.values()),
                "pixel_spacing": (0.78125, 0.78125), "ipp_xy": (-10.0, -10.0), "n_slices": len(gz)}
    specs = [
        RoiSpec(name="ISO_100%_20.0Gy", color=(255, 0, 255), contours_by_z=dm.mask_to_contours(sphere, grid)),
        RoiSpec(name="iso_100%_20.0gy", color=(0, 255, 255), contours_by_z=cz, kind="inter"),
    ]
    _unique_names(specs)
    with tempfile.TemporaryDirectory(prefix="rtstruct_writer_selftest_") as tmp:
        out = Path(tmp) / "RS_test.dcm"
        ds = write_isodose_rtstruct(orig, ct_index, specs, out, label="_IDX", description="Self-Test")
        problems = verify_rtstruct(out, ct_index)
        back = pydicom.dcmread(str(out))
        check("Schreiben/Pruefen: verify_rtstruct ohne Befund", not problems, "; ".join(problems))
        check("Schreiben: Label 'Phantom_Original' + '_IDX' auf 16 Zeichen gekuerzt",
              back.StructureSetLabel == "Phantom_Orig_IDX", back.StructureSetLabel)
        check("Schreiben: Namenskollision case-insensitiv -> Suffix _2",
              [str(r.ROIName) for r in back.StructureSetROISequence] == ["ISO_100%_20.0Gy", "iso_100%_20.0gy_2"])
        n_c = sum(len(rc.ContourSequence) for rc in back.ROIContourSequence)
        check("Schreiben: Konturen vorhanden, SOPInstanceUID == file_meta", n_c > 10
              and str(back.file_meta.MediaStorageSOPInstanceUID) == str(back.SOPInstanceUID), f"{n_c} Konturen")
        vol = ana.compute_volume(ana.extract_contours(back, 1))
        check("Analyzer-Volumen der Kugel-ROI ~ Maskenvolumen (2 %)",
              abs(vol - dm.mask_volume_cm3(sphere, grid)) < 0.02 * dm.mask_volume_cm3(sphere, grid),
              f"{vol:.4f} vs {dm.mask_volume_cm3(sphere, grid):.4f} cm3")

    n_fail = results.count(False)
    print("-" * 70)
    print(f"Gesamt: {'PASS' if n_fail == 0 else 'FAIL'}  ({len(results) - n_fail}/{len(results)} Pruefungen bestanden)")
    return 0 if n_fail == 0 else 1


def main(argv: Optional[list] = None) -> int:
    import argparse
    p = argparse.ArgumentParser(description="RTSTRUCT-Writer fuer Isodosen-/Indexkonturen (nur Self-Test)")
    p.add_argument("--self-test", action="store_true", help="Synthetische Masken schreiben und pruefen")
    args = p.parse_args(argv)
    if args.self_test:
        return _run_self_test()
    p.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
