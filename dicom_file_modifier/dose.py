"""
dose.py - Numerik-Bibliothek fuer Dosisindizes (RTDOSE-Gitter, Feinraster,
Konturrasterung, Dosis-Sampling, Isodosenmasken, DVH-Statistik).

Kein eigenstaendiges Werkzeug (ausser ``--info``); wird von ``dose_indices``
und ``rtstruct_writer`` benutzt.

Konventionen:
  - Dosisarray ``(k, j, i)`` = (Frame, Zeile, Spalte), Affine wie
    ``modifier.extract_geometry``:  P_patient = A @ [k, j, i, 1]^T
  - Feingitter-Achsen ``gx, gy, gz`` sind Voxel-MITTELPUNKTE in mm (LPS),
    Masken haben die Achsenreihenfolge ``(z, y, x)``.
  - Konturen werden mit ``analyzer.rasterize_contours`` (XOR je Ring, Even-Odd)
    gerastert; dieses Modul bringt keinen eigenen Rasterizer mit.
  - Alle Konsolenausgaben sind ASCII (Windows-cp1252-Konsole).

Verwendung (nur Gitterinfo):
  python -m dicom_file_modifier.dose <RD.dcm>
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pydicom
from scipy import ndimage
from scipy.integrate import trapezoid

from . import analyzer as ana


# ---------------------------------------------------------------------------
# 1. RTDOSE laden
# ---------------------------------------------------------------------------

@dataclass
class DoseGrid:
    """RTDOSE als Gy-Array ``(k, j, i)`` mit Voxel->Patient-Affine."""
    array: np.ndarray                  # float32 (nk, nj, ni), pixel * DoseGridScaling
    affine: np.ndarray                 # 4x4, P = affine @ [k, j, i, 1]
    spacing: tuple                     # (dz, dr, dc) in mm, alle > 0
    origin: np.ndarray                 # (3,) = affine[:3, 3]
    units: str                         # "GY"
    dose_type: str                     # "PHYSICAL" | "EFFECTIVE" | ...
    summation_type: str                # "PLAN" | "BEAM" | "FRACTION" | ...
    dmax: float                        # globales Maximum in Gy (natives Gitter)
    frame_of_reference_uid: str
    sop_instance_uid: str
    referenced_plan_uid: Optional[str] = None
    referenced_rtstruct_uid: Optional[str] = None
    gfov_mode: str = "relative"        # "relative" | "absolute"
    source_path: Optional[str] = None
    warnings: list = field(default_factory=list)

    @property
    def shape(self) -> tuple:
        return tuple(self.array.shape)

    @property
    def normal(self) -> np.ndarray:
        """Einheitsvektor der k-Achse (Schichtnormale)."""
        v = self.affine[:3, 0]
        return v / np.linalg.norm(v)

    def patient_to_index(self, pts_xyz: np.ndarray) -> np.ndarray:
        """(N,3) Patienten-mm -> (N,3) fraktionale Voxelindizes (k, j, i)."""
        pts = np.atleast_2d(np.asarray(pts_xyz, dtype=float))
        homo = np.column_stack([pts, np.ones(len(pts))])
        return (np.linalg.inv(self.affine) @ homo.T).T[:, :3]

    def index_to_patient(self, kji: np.ndarray) -> np.ndarray:
        """(N,3) Voxelindizes (k, j, i) -> (N,3) Patienten-mm."""
        idx = np.atleast_2d(np.asarray(kji, dtype=float))
        homo = np.column_stack([idx, np.ones(len(idx))])
        return (self.affine @ homo.T).T[:, :3]

    def plane_z_values(self) -> np.ndarray:
        """z-Koordinaten der Dosisebenen (aufsteigend); nur fuer axiale Gitter."""
        if abs(abs(self.normal[2]) - 1.0) > 1e-3:
            raise ValueError("Dosisgitter ist nicht axial (Schichtnormale nicht +-z).")
        nk = self.shape[0]
        zs = self.affine[2, 3] + np.arange(nk) * self.affine[2, 0]
        return np.sort(zs)


def load_rtdose(path: str) -> pydicom.Dataset:
    """Liest eine RTDOSE-Datei (inkl. Pixel) und prueft die Modalitaet."""
    ds = pydicom.dcmread(str(path))
    if str(ds.get("Modality", "")).upper() != "RTDOSE":
        raise ValueError(f"Datei ist keine RTDOSE (Modalitaet: {ds.get('Modality')})")
    return ds


def load_rtplan(path: str) -> pydicom.Dataset:
    """Liest eine RTPLAN-Datei (ohne Pixel) und prueft die Modalitaet."""
    ds = pydicom.dcmread(str(path), stop_before_pixels=True)
    if str(ds.get("Modality", "")).upper() != "RTPLAN":
        raise ValueError(f"Datei ist keine RTPLAN (Modalitaet: {ds.get('Modality')})")
    return ds


def dose_grid_from_dataset(ds: pydicom.Dataset, source_path: Optional[str] = None) -> DoseGrid:
    """
    Baut den ``DoseGrid`` aus IOP/IPP/PixelSpacing/GridFrameOffsetVector/
    DoseGridScaling.  Validiert Einheiten, Frame-Anzahl und GFOV-Schritt
    (``ValueError``); Warnungen (Summationstyp, GFOV-Konvention, Orientierung)
    landen in ``DoseGrid.warnings``.
    """
    warnings = []
    units = str(ds.get("DoseUnits", "")).upper()
    if units != "GY":
        raise ValueError(
            f"DoseUnits '{units}' wird nicht unterstuetzt (nur GY; RELATIVE-Dosen "
            "muessten erst mit der Verschreibung skaliert werden)."
        )
    dose_type = str(ds.get("DoseType", "")).upper()
    summation = str(ds.get("DoseSummationType", "")).upper()
    if summation not in ("PLAN", "MULTI_PLAN"):
        warnings.append(
            f"DoseSummationType '{summation}': Indizes beziehen sich auf eine "
            "Teildosis (kein Plan-Summendosis)."
        )

    n_frames = int(ds.get("NumberOfFrames", 1))
    gfov = np.asarray(ds.get("GridFrameOffsetVector", [0.0]), dtype=float).reshape(-1)
    if len(gfov) != n_frames:
        raise ValueError(
            f"GridFrameOffsetVector hat {len(gfov)} Eintraege, NumberOfFrames ist {n_frames}."
        )
    if n_frames < 2:
        raise ValueError("RTDOSE mit nur einem Frame wird nicht unterstuetzt.")

    steps = np.diff(gfov)
    step = float(np.mean(steps))
    if step == 0.0 or np.any(np.abs(steps - step) > 1e-3 * abs(step)):
        raise ValueError(
            "GridFrameOffsetVector ist nicht aequidistant "
            f"(Schritte {np.unique(np.round(steps, 4)).tolist()} mm)."
        )

    iop = np.asarray([float(v) for v in ds.ImageOrientationPatient], dtype=float)
    ipp = np.asarray([float(v) for v in ds.ImagePositionPatient], dtype=float)
    row_dir, col_dir = iop[:3], iop[3:]
    normal = np.cross(row_dir, col_dir)
    ps = ds.PixelSpacing
    dr, dc = float(ps[0]), float(ps[1])

    # GFOV-Konvention: relativ zu IPP (erster Wert 0) oder absolute z-Werte.
    if abs(gfov[0]) < 1e-6:
        gfov_mode = "relative"
    elif abs(gfov[0] - ipp[2]) < 1e-3 and abs(abs(normal[2]) - 1.0) < 1e-6:
        gfov_mode = "absolute"
        gfov = gfov - gfov[0]
    else:
        gfov_mode = "relative"
        warnings.append(
            f"GridFrameOffsetVector beginnt bei {gfov[0]:g} (weder 0 noch IPP-z); "
            "wird als relativer Offset interpretiert."
        )
    if np.any(np.abs(np.round(iop) - iop) > 1e-6) or np.any(np.abs(iop) % 1 > 1e-6):
        warnings.append("ImageOrientationPatient ist nicht achsparallel; Sampling ueber die Affine.")

    A = np.eye(4)
    A[:3, 0] = normal * step
    A[:3, 1] = col_dir * dr
    A[:3, 2] = row_dir * dc
    A[:3, 3] = ipp + normal * gfov[0]

    scaling = float(ds.DoseGridScaling)
    arr = ds.pixel_array
    if arr.ndim == 2:
        arr = arr[None, :, :]
    array = arr.astype(np.float32) * np.float32(scaling)

    def _ref_uid(seq_name: str) -> Optional[str]:
        seq = ds.get(seq_name)
        if seq:
            return str(seq[0].get("ReferencedSOPInstanceUID", "")) or None
        return None

    return DoseGrid(
        array=array,
        affine=A,
        spacing=(abs(step), dr, dc),
        origin=A[:3, 3].copy(),
        units=units,
        dose_type=dose_type,
        summation_type=summation,
        dmax=float(array.max()),
        frame_of_reference_uid=str(ds.get("FrameOfReferenceUID", "")),
        sop_instance_uid=str(ds.get("SOPInstanceUID", "")),
        referenced_plan_uid=_ref_uid("ReferencedRTPlanSequence"),
        referenced_rtstruct_uid=_ref_uid("ReferencedStructureSetSequence"),
        gfov_mode=gfov_mode,
        source_path=str(source_path) if source_path else None,
        warnings=warnings,
    )


def validate_dose_against_rtstruct(dose: DoseGrid, rs_ds: pydicom.Dataset,
                                   rp_ds: Optional[pydicom.Dataset] = None) -> list:
    """
    FoR-Mismatch RD<->RS -> ``ValueError``; abweichende Referenz-UIDs
    (RD->RS, RD->RP, RP->RS) -> Warnungstexte.
    """
    from .case_modifier import get_rs_frame_of_references

    warnings = []
    rs_fors = get_rs_frame_of_references(rs_ds)
    if dose.frame_of_reference_uid and rs_fors and dose.frame_of_reference_uid not in rs_fors:
        raise ValueError(
            "FrameOfReferenceUID-Mismatch: RTDOSE "
            f"{dose.frame_of_reference_uid} vs RTSTRUCT {sorted(rs_fors)}"
        )
    rs_uid = str(rs_ds.get("SOPInstanceUID", ""))
    if dose.referenced_rtstruct_uid and rs_uid and dose.referenced_rtstruct_uid != rs_uid:
        warnings.append(
            "RTDOSE referenziert ein anderes Structure Set "
            f"({dose.referenced_rtstruct_uid[-16:]}... vs {rs_uid[-16:]}...)."
        )
    if rp_ds is not None:
        rp_uid = str(rp_ds.get("SOPInstanceUID", ""))
        if dose.referenced_plan_uid and rp_uid and dose.referenced_plan_uid != rp_uid:
            warnings.append("RTDOSE referenziert einen anderen RTPLAN als die uebergebene RP-Datei.")
        rp_rs = rp_ds.get("ReferencedStructureSetSequence")
        if rp_rs and rs_uid and str(rp_rs[0].get("ReferencedSOPInstanceUID", "")) != rs_uid:
            warnings.append("RTPLAN referenziert ein anderes Structure Set als die uebergebene RS-Datei.")
    return warnings


def prescription_references(rp_ds: pydicom.Dataset) -> list:
    """
    ``DoseReferenceSequence`` -> Liste von Dicts mit ``number, description,
    structure_type, reference_type, target_prescription_dose_gy,
    referenced_roi_number`` (fehlende Felder = None).
    """
    out = []
    for dr in rp_ds.get("DoseReferenceSequence", []):
        tp = dr.get("TargetPrescriptionDose")
        out.append({
            "number": int(dr.get("DoseReferenceNumber", 0) or 0),
            "description": str(dr.get("DoseReferenceDescription", "") or ""),
            "structure_type": str(dr.get("DoseReferenceStructureType", "") or ""),
            "reference_type": str(dr.get("DoseReferenceType", "") or ""),
            "target_prescription_dose_gy": float(tp) if tp is not None else None,
            "referenced_roi_number": (int(dr.ReferencedROINumber)
                                      if "ReferencedROINumber" in dr else None),
        })
    return out


def fractions_planned(rp_ds: pydicom.Dataset) -> Optional[int]:
    """Anzahl geplanter Fraktionen aus der ersten FractionGroup (oder None)."""
    for fg in rp_ds.get("FractionGroupSequence", []):
        n = fg.get("NumberOfFractionsPlanned")
        if n is not None:
            return int(n)
    return None


# ---------------------------------------------------------------------------
# 1b. Eclipse-DVH (DVHSequence der RTDOSE)
# ---------------------------------------------------------------------------

@dataclass
class EclipseDVH:
    """
    Kumulatives DVH einer ROI aus der RTDOSE-``DVHSequence`` (Eclipse-Export).

    ``dose_gy[i]`` ist die LINKE Kante von Bin i (Bin 0 bei 0 Gy),
    ``volume_cm3[i]`` das Volumen mit D >= ``dose_gy[i]``.  Linke Kante,
    weil Bin 0 das Gesamtvolumen traegt (V(D >= 0) = TV) und der letzte Bin
    mit Volumen dann exakt bei ``DVHMaximumDose`` liegt; bei den ueblichen
    0.01-Gy-Bins wuerde die rechte Kante V(Rx) um ~0.0003 cm3 verschieben.
    """
    roi_number: int
    dose_gy: np.ndarray
    volume_cm3: np.ndarray
    total_volume_cm3: float
    dvh_type: str
    dose_units: str
    volume_units: str
    dmin_gy: Optional[float]
    dmax_gy: Optional[float]
    dmean_gy: Optional[float]
    n_bins: int
    bin_width_gy: float

    def v_at(self, dose_gy: float) -> float:
        """Volumen (cm3) mit D >= ``dose_gy`` (linear zwischen den Bins)."""
        return float(np.interp(float(dose_gy), self.dose_gy, self.volume_cm3,
                               left=self.total_volume_cm3, right=0.0))

    def d_at(self, volume_pct: float) -> float:
        """
        D_x: Dosis, die mindestens ``volume_pct`` % des Volumens erhaelt
        (Plateau-Konvention wie ``dose_at_volume_fraction``: letzter Bin, der
        das Volumen noch erreicht, dann linear zum naechsten Bin).
        """
        vol, dose = self.volume_cm3, self.dose_gy
        if len(vol) == 0:
            return float("nan")
        v = volume_pct / 100.0 * self.total_volume_cm3
        if v <= 0:
            nz = np.nonzero(vol > 0)[0]
            return float(dose[nz[-1]]) if len(nz) else float(dose[0])
        idx = np.nonzero(vol >= v)[0]
        if len(idx) == 0:
            return float(dose[0])
        i = int(idx[-1])
        if i + 1 >= len(dose) or vol[i] <= v:
            return float(dose[i])
        v0, v1, d0, d1 = vol[i], vol[i + 1], dose[i], dose[i + 1]
        if v0 == v1:
            return float(d0)
        return float(d0 + (v0 - v) / (v0 - v1) * (d1 - d0))


def read_dvh_sequence(rd_ds: pydicom.Dataset) -> tuple:
    """
    ``DVHSequence`` der RTDOSE -> ``({roi_number: EclipseDVH}, hinweise)``.
    Tolerant: keine Sequenz -> ``({}, [])``; Items, die nicht CUMULATIVE/GY/CM3
    sind, mehrere ROIs kombinieren (``DVHReferencedROISequence`` != 1 Eintrag),
    ``EXCLUDED`` beitragen oder leere/ungerade ``DVHData`` haben, werden mit
    Hinweis uebersprungen.  ``DVHDoseScaling`` wird auf die Bin-Breiten
    angewendet; ``DVHMinimumDose``/``DVHMaximumDose``/``DVHMeanDose`` werden
    unveraendert uebernommen (None, wenn nicht vorhanden).
    """
    out, notes = {}, []
    seq = rd_ds.get("DVHSequence", None)
    if not seq:
        return out, notes
    for n, item in enumerate(seq, start=1):
        refs = item.get("DVHReferencedROISequence", []) or []
        if len(refs) != 1:
            notes.append(f"DVH #{n}: {len(refs)} referenzierte ROIs (nur genau eine wird "
                         "unterstuetzt); uebersprungen.")
            continue
        ref = refs[0]
        roi = int(ref.get("ReferencedROINumber", 0) or 0)
        contrib = str(ref.get("DVHROIContributionType", "INCLUDED") or "INCLUDED").upper()
        dvh_type = str(item.get("DVHType", "") or "").upper()
        units = str(item.get("DoseUnits", "") or "").upper()
        vunits = str(item.get("DVHVolumeUnits", "") or "").upper()
        if contrib != "INCLUDED":
            notes.append(f"DVH #{n} (ROI {roi}): Beitragstyp {contrib}; uebersprungen.")
            continue
        if dvh_type != "CUMULATIVE":
            notes.append(f"DVH #{n} (ROI {roi}): Typ {dvh_type or '?'} statt CUMULATIVE; uebersprungen.")
            continue
        if units != "GY":
            notes.append(f"DVH #{n} (ROI {roi}): DoseUnits {units or '?'} statt GY; uebersprungen.")
            continue
        if vunits != "CM3":
            notes.append(f"DVH #{n} (ROI {roi}): Volumeneinheit {vunits or '?'} statt CM3; uebersprungen.")
            continue
        data = np.asarray(item.get("DVHData", []) or [], dtype=float).reshape(-1)
        if data.size == 0 or data.size % 2:
            notes.append(f"DVH #{n} (ROI {roi}): DVHData leer oder ungerade ({data.size} Werte); uebersprungen.")
            continue
        if roi in out:
            notes.append(f"DVH #{n}: zweites DVH fuer ROI {roi} ignoriert.")
            continue
        scaling = float(item.get("DVHDoseScaling", 1.0) or 1.0)
        pairs = data.reshape(-1, 2)
        widths = pairs[:, 0] * scaling
        volumes = pairs[:, 1].astype(float)
        dose = np.concatenate([[0.0], np.cumsum(widths)[:-1]])

        def _opt(tag):
            v = item.get(tag)
            return float(v) if v is not None else None

        out[roi] = EclipseDVH(
            roi_number=roi, dose_gy=dose, volume_cm3=volumes,
            total_volume_cm3=float(volumes[0]), dvh_type=dvh_type, dose_units=units,
            volume_units=vunits, dmin_gy=_opt("DVHMinimumDose"), dmax_gy=_opt("DVHMaximumDose"),
            dmean_gy=_opt("DVHMeanDose"), n_bins=int(len(widths)),
            bin_width_gy=float(np.median(widths)) if len(widths) else 0.0,
        )
    return out, notes


def dvh_statistics(dvh: EclipseDVH, rx_gy: float) -> dict:
    """
    Kennwerte eines Eclipse-DVH: ``total_cm3, v_rx_cm3, v_half_rx_cm3, d2_gy,
    d50_gy, d95_gy, d98_gy, dmin_gy, dmax_gy, dmean_gy``.  Dmin/Dmax/Dmean
    kommen aus den DICOM-Attributen, sonst aus der Kurve (D100, D0, Trapez).
    """
    tot = dvh.total_volume_cm3
    dmean = dvh.dmean_gy
    if dmean is None and tot > 0 and len(dvh.dose_gy) > 1:
        dmean = float(trapezoid(dvh.volume_cm3, dvh.dose_gy) / tot)
    return {
        "total_cm3": tot,
        "v_rx_cm3": dvh.v_at(rx_gy),
        "v_half_rx_cm3": dvh.v_at(0.5 * rx_gy),
        "d2_gy": dvh.d_at(2.0), "d50_gy": dvh.d_at(50.0),
        "d95_gy": dvh.d_at(95.0), "d98_gy": dvh.d_at(98.0),
        "dmin_gy": dvh.dmin_gy if dvh.dmin_gy is not None else dvh.d_at(100.0),
        "dmax_gy": dvh.dmax_gy if dvh.dmax_gy is not None else dvh.d_at(0.0),
        "dmean_gy": dmean,
    }


# ---------------------------------------------------------------------------
# 2. Feingitter
# ---------------------------------------------------------------------------

@dataclass
class FineGrid:
    """Achsparalleles Auswertegitter; ``gx, gy, gz`` = Voxelmittelpunkte (mm)."""
    gx: np.ndarray
    gy: np.ndarray
    gz: np.ndarray
    res_xy: float
    dz: float

    @property
    def shape(self) -> tuple:
        return (len(self.gz), len(self.gy), len(self.gx))

    @property
    def n_voxels(self) -> int:
        return int(np.prod(self.shape))

    @property
    def voxel_volume_mm3(self) -> float:
        return float(self.res_xy * self.res_xy * self.dz)

    @property
    def affine(self) -> np.ndarray:
        """4x4 mit P = A @ [k, j, i, 1] (kompatibel mit ``modifier._extract_surface``)."""
        A = np.zeros((4, 4))
        A[3, 3] = 1.0
        A[2, 0] = self.dz          # k -> z
        A[1, 1] = self.res_xy      # j -> y
        A[0, 2] = self.res_xy      # i -> x
        A[:3, 3] = [self.gx[0], self.gy[0], self.gz[0]]
        return A

    @property
    def bbox(self) -> tuple:
        """((xmin, ymin, zmin), (xmax, ymax, zmax)) der Voxel-Aussenkanten."""
        h, hz = 0.5 * self.res_xy, 0.5 * self.dz
        lo = (self.gx[0] - h, self.gy[0] - h, self.gz[0] - hz)
        hi = (self.gx[-1] + h, self.gy[-1] + h, self.gz[-1] + hz)
        return lo, hi

    def plane_index(self, z: float, tol: Optional[float] = None) -> Optional[int]:
        """Index der Ebene mit |gz - z| <= tol (Default 0.05*dz), sonst None."""
        if tol is None:
            tol = 0.05 * self.dz
        k = int(np.argmin(np.abs(self.gz - z)))
        return k if abs(self.gz[k] - z) <= tol else None

    def plane_points(self, k: int) -> np.ndarray:
        """(ny*nx, 3) Patientenkoordinaten aller Voxelmittelpunkte der Ebene k."""
        X, Y = np.meshgrid(self.gx, self.gy)          # (ny, nx)
        return np.column_stack([X.ravel(), Y.ravel(), np.full(X.size, self.gz[k])])


def contours_bbox(contour_sets: list) -> tuple:
    """Gemeinsame BBox mehrerer Konturlisten -> (lo(3,), hi(3,))."""
    pts = np.vstack([c for cs in contour_sets for c in cs])
    return pts.min(axis=0), pts.max(axis=0)


def native_level_bbox(dose: DoseGrid, level_gy: float, margin_voxels: int = 2) -> Optional[tuple]:
    """
    BBox (Patienten-mm) aller nativen Voxel mit D >= level, um ``margin_voxels``
    erweitert; ``None`` wenn kein Voxel das Level erreicht.
    """
    m = dose.array >= level_gy
    if not m.any():
        return None
    lo_idx, hi_idx = [], []
    for axis in range(3):
        other = tuple(a for a in range(3) if a != axis)
        nz = np.flatnonzero(m.any(axis=other))
        lo_idx.append(max(int(nz[0]) - margin_voxels, 0))
        hi_idx.append(min(int(nz[-1]) + margin_voxels, m.shape[axis] - 1))
    corners = np.array([[k, j, i] for k in (lo_idx[0], hi_idx[0])
                        for j in (lo_idx[1], hi_idx[1]) for i in (lo_idx[2], hi_idx[2])], float)
    pts = dose.index_to_patient(corners)
    return pts.min(axis=0), pts.max(axis=0)


def build_fine_grid(dose: DoseGrid, bbox_lo, bbox_hi, res_xy: float,
                    contour_z: Optional[np.ndarray] = None,
                    restrict_z_to: Optional[np.ndarray] = None,
                    margin_mm: float = 1.0, max_voxels: int = 40_000_000,
                    align: Optional[tuple] = None) -> FineGrid:
    """
    Feingitter ueber die BBox: In-Plane-Achsen an 1-mm-Vielfache gesnappt
    (Raster 1.0/0.5/0.25/0.1 nisten ineinander), Mittelpunkte bei
    ``start + (i + 0.5) * res``; z = native Dosisebenen im Bereich.
    ``align=(ax, ay)``: Voxelmittelpunkte stattdessen auf das Gitter
    ``ax + i*res`` / ``ay + j*res`` legen (z.B. CT-Pixelzentren, Eclipse-Modus).

    ``contour_z``: Konturebenen, die auf Gitterebenen liegen muessen; sonst wird
    die z-Achse ganzzahlig verfeinert (dz/n, n <= 8) oder ``ValueError``.
    ``restrict_z_to``: nur Dosisebenen behalten, die (innerhalb 0.05*dz) in
    dieser Liste vorkommen (z.B. die CT-Schichtpositionen).
    """
    lo = np.asarray(bbox_lo, dtype=float) - margin_mm
    hi = np.asarray(bbox_hi, dtype=float) + margin_mm

    def _axis(a0: float, a1: float, origin: Optional[float]) -> np.ndarray:
        if origin is None:
            start = math.floor(a0)
            end = math.ceil(a1)
            n = int(round((end - start) / res_xy))
            return start + (np.arange(n) + 0.5) * res_xy
        i0 = math.floor((a0 - origin) / res_xy)
        i1 = math.ceil((a1 - origin) / res_xy)
        return origin + np.arange(i0, i1 + 1) * res_xy

    gx = _axis(lo[0], hi[0], align[0] if align else None)
    gy = _axis(lo[1], hi[1], align[1] if align else None)

    zs = dose.plane_z_values()
    dz = float(dose.spacing[0])
    if restrict_z_to is not None and len(restrict_z_to):
        allowed = np.asarray(restrict_z_to, dtype=float)
        keep = np.array([np.min(np.abs(allowed - z)) <= 0.05 * dz for z in zs])
        zs = zs[keep]
    sel = (zs >= lo[2] - 0.5 * dz) & (zs <= hi[2] + 0.5 * dz)
    gz = zs[sel]
    if len(gz) == 0:
        raise ValueError("Keine Dosisebene im Bereich des Zielvolumens.")

    # Konturebenen muessen Gitterebenen treffen; sonst ganzzahlige Verfeinerung.
    if contour_z is not None and len(contour_z):
        cz = np.asarray(contour_z, dtype=float)
        gz_use, dz_use = gz, dz
        for n in range(1, 9):
            dz_try = dz / n
            gz_try = gz[0] + np.arange(int(round((gz[-1] - gz[0]) / dz_try)) + 1) * dz_try
            if all(np.min(np.abs(gz_try - z)) <= 0.05 * dz_try for z in cz):
                gz_use, dz_use = gz_try, dz_try
                break
        else:
            raise ValueError(
                "Konturebenen liegen nicht auf den Dosisebenen (auch nicht nach "
                "Verfeinerung bis 1/8): z.B. Kontur z="
                f"{cz[0]:.3f}, Dosisebenen ab {gz[0]:.3f} mit dz={dz:.3f}."
            )
        gz, dz = gz_use, dz_use

    grid = FineGrid(gx=gx, gy=gy, gz=gz, res_xy=float(res_xy), dz=float(dz))
    if grid.n_voxels > max_voxels:
        raise ValueError(
            f"Feingitter zu gross ({grid.n_voxels / 1e6:.1f} M Voxel > "
            f"{max_voxels / 1e6:.0f} M). Groeberes Raster waehlen (--grid 0.5 / 1.0)."
        )
    return grid


# ---------------------------------------------------------------------------
# 3. Konturen rastern (XOR via analyzer.rasterize_contours) + Slab-Gewichte
# ---------------------------------------------------------------------------

@dataclass
class StructureMask:
    """Gerasterte Konturstruktur auf dem Feingitter."""
    name: str
    roi_number: Optional[int]
    mask: np.ndarray                   # bool (nz, ny, nx), Even-Odd (Loecher leer)
    slab_w: np.ndarray                 # float (nz,): 1.0 / 0.5 (eclipse-Endschichten) / 0.0
    planes: np.ndarray                 # int: Gitterebenen mit Struktur
    n_contours: int
    n_slices: int
    n_holes: int
    n_gaps: int
    slice_spacing_mm: Optional[float]
    volume_model: str

    def weighted(self) -> np.ndarray:
        """float32 (nz, ny, nx) = mask * slab_w[:, None, None]."""
        return self.mask.astype(np.float32) * self.slab_w.astype(np.float32)[:, None, None]

    def volume_cm3(self, grid: FineGrid) -> float:
        return mask_volume_cm3(self.weighted(), grid)


def closed_planar_contours(rs_ds: pydicom.Dataset, roi_number: int) -> list:
    """Wie ``analyzer.extract_contours``, aber nur CLOSED_PLANAR(_XOR), float64, >= 3 Punkte."""
    out = []
    for rc in rs_ds.ROIContourSequence:
        if int(rc.ReferencedROINumber) != int(roi_number):
            continue
        for c in rc.get("ContourSequence", []):
            gt = str(c.get("ContourGeometricType", "")).upper()
            if not gt.startswith("CLOSED_PLANAR"):
                continue
            pts = np.asarray(c.ContourData, dtype=float).reshape(-1, 3)
            if len(pts) >= 3:
                out.append(pts)
    return out


def _slice_runs(slice_z: np.ndarray, dz: Optional[float]) -> list:
    """Zusammenhaengende z-Laeufe (Indizes in ``slice_z``); Luecke = Abstand > 1.5*dz."""
    if len(slice_z) == 0:
        return []
    if dz is None or len(slice_z) == 1:
        return [list(range(len(slice_z)))]
    runs, cur = [], [0]
    for i in range(1, len(slice_z)):
        if slice_z[i] - slice_z[i - 1] > 1.5 * dz:
            runs.append(cur)
            cur = [i]
        else:
            cur.append(i)
    runs.append(cur)
    return runs


def slab_weights(slice_z: np.ndarray, slice_planes: list, nz: int,
                 dz: Optional[float], volume_model: str) -> np.ndarray:
    """
    (nz,) Gewichte: 1.0 auf allen Gitterebenen einer Kontur-Schicht; unter
    ``volume_model="eclipse"`` 0.5 auf den Ebenen der ersten und letzten
    Schicht jedes zusammenhaengenden z-Laufs (Endschichten halb); 0.0 sonst.
    ``slice_planes[i]`` = Liste der Gitterebenen-Indizes der Schicht i.
    """
    if volume_model not in ("slab", "eclipse"):
        raise ValueError(f"Unbekanntes Volumenmodell: {volume_model!r}")
    w = np.zeros(nz, dtype=float)
    for planes in slice_planes:
        w[planes] = 1.0
    if volume_model == "eclipse":
        for run in _slice_runs(np.asarray(slice_z, float), dz):
            for end in (run[0], run[-1]):
                w[slice_planes[end]] = 0.5
    return w


def rasterize_structure(contours: list, grid: FineGrid, volume_model: str = "slab",
                        name: str = "", roi_number: Optional[int] = None) -> StructureMask:
    """
    Rastert eine Kontur-Liste auf das Feingitter (``analyzer.rasterize_contours``,
    XOR je Ring), extrudiert jede Kontur-Schicht auf alle Gitterebenen ihres
    Slabs (nur relevant bei verfeinerter z-Achse) und berechnet Slab-Gewichte.
    Konturebenen ohne Gitterebene -> ``ValueError``.
    """
    slices = ana._group_slices(contours)
    slice_z = np.array([z for z, _ in slices], dtype=float)
    dz_c, n_gaps = ana._nominal_slice_spacing(slice_z)
    nz = len(grid.gz)

    slice_planes = []
    for z in slice_z:
        k = grid.plane_index(z)
        if k is None:
            raise ValueError(
                f"Konturebene z={z:.3f} mm von '{name}' liegt auf keiner Gitterebene "
                f"(dz={grid.dz:.3f} mm)."
            )
        if dz_c is not None and grid.dz < dz_c - 1e-6:
            half = 0.5 * dz_c
            ks = np.flatnonzero((grid.gz > z - half - 1e-6) & (grid.gz <= z + half + 1e-6))
            slice_planes.append(sorted(set(ks.tolist()) | {k}))
        else:
            slice_planes.append([k])

    mask = ana.rasterize_contours(contours, grid.gx, grid.gy, grid.gz)
    # Extrusion auf alle Ebenen des Slabs (bei verfeinerter z-Achse)
    for z, planes in zip(slice_z, slice_planes):
        if len(planes) > 1:
            k0 = grid.plane_index(z)
            for k in planes:
                if k != k0:
                    mask[k] |= mask[k0]

    geoms = ana._slice_geometries(contours)
    n_holes = 0
    for _, g in geoms:
        parts = getattr(g, "geoms", [g])
        n_holes += sum(len(p.interiors) for p in parts if p.geom_type == "Polygon")

    w = slab_weights(slice_z, slice_planes, nz, dz_c, volume_model)
    planes = np.array(sorted({k for pl in slice_planes for k in pl}), dtype=int)
    return StructureMask(
        name=name, roi_number=roi_number, mask=mask, slab_w=w, planes=planes,
        n_contours=len(contours), n_slices=len(slices), n_holes=n_holes,
        n_gaps=n_gaps, slice_spacing_mm=dz_c, volume_model=volume_model,
    )


def planimetric_volume_cm3(contours: list, volume_model: str = "slab") -> float:
    """
    Exaktes Even-Odd-Planimetrievolumen (Kreuzcheck fuer das Rastervolumen):
    ``slab`` = ``analyzer.compute_volume``; ``eclipse`` = Endschichten jedes
    z-Laufs halb gewichtet.
    """
    geoms = ana._slice_geometries(contours)
    if volume_model == "slab":
        return ana.compute_volume(contours, geoms=geoms)
    if volume_model != "eclipse":
        raise ValueError(f"Unbekanntes Volumenmodell: {volume_model!r}")
    slice_z = np.array([z for z, _ in ana._group_slices(contours)], dtype=float)
    dz, _ = ana._nominal_slice_spacing(slice_z)
    if dz is None or len(geoms) < 2:
        return 0.0
    area_by_z = {z: g.area for z, g in geoms}
    w = np.ones(len(slice_z))
    for run in _slice_runs(slice_z, dz):
        w[run[0]] = 0.5
        w[run[-1]] = 0.5
    total = sum(w[i] * area_by_z.get(z, 0.0) for i, z in enumerate(slice_z))
    return total * dz / 1000.0


# ---------------------------------------------------------------------------
# 4. Dosis-Sampling, Isodosen, Komponenten, Volumina
# ---------------------------------------------------------------------------

def _crop_for_sampling(dose: DoseGrid, pts_xyz: np.ndarray, order: int):
    """Crop des nativen Arrays um die Punkte (+ Rand) -> (vol, offset_kji)."""
    idx = dose.patient_to_index(pts_xyz)
    margin = 2 if order <= 1 else 10
    lo = np.floor(idx.min(axis=0)).astype(int) - margin
    hi = np.ceil(idx.max(axis=0)).astype(int) + margin
    lo = np.clip(lo, 0, np.array(dose.shape) - 1)
    hi = np.clip(hi, 0, np.array(dose.shape) - 1)
    crop = dose.array[lo[0]:hi[0] + 1, lo[1]:hi[1] + 1, lo[2]:hi[2] + 1]
    if order > 1:
        vol = ndimage.spline_filter(crop.astype(np.float64), order=order, mode="mirror")
    else:
        vol = crop.astype(np.float64)
    return vol, lo


def _sample(vol: np.ndarray, offset: np.ndarray, dose: DoseGrid,
            pts_xyz: np.ndarray, order: int) -> np.ndarray:
    idx = dose.patient_to_index(pts_xyz)
    shape = np.array(dose.shape, dtype=float)
    outside = np.any((idx < -0.5) | (idx > shape - 0.5), axis=1)
    local = (idx - offset).T
    vals = ndimage.map_coordinates(vol, local, order=order, mode="nearest", prefilter=False)
    vals = vals.astype(np.float32)
    vals[outside] = np.nan
    return vals


def sample_dose_on_grid(dose: DoseGrid, grid: FineGrid, order: int = 1) -> np.ndarray:
    """
    Dosis (Gy, float32, ``(nz, ny, nx)``) an den Feingitter-Mittelpunkten;
    NaN ausserhalb des nativen Gitters.  ``order`` 1 = trilinear, 3 = kubischer
    B-Spline (Prefilter einmal auf dem Crop, wie ``modifier.resample_volume``).
    """
    if order not in (1, 3):
        raise ValueError("order muss 1 (linear) oder 3 (kubisch) sein.")
    lo, hi = grid.bbox
    corners = np.array([[x, y, z] for x in (lo[0], hi[0]) for y in (lo[1], hi[1])
                        for z in (lo[2], hi[2])], dtype=float)
    vol, offset = _crop_for_sampling(dose, corners, order)
    out = np.empty(grid.shape, dtype=np.float32)
    ny, nx = len(grid.gy), len(grid.gx)
    for k in range(len(grid.gz)):
        out[k] = _sample(vol, offset, dose, grid.plane_points(k), order).reshape(ny, nx)
    return out


def sample_dose_at_points(dose: DoseGrid, pts_xyz: np.ndarray, order: int = 1) -> np.ndarray:
    """Dosis (Gy) an beliebigen Punkten (N,), NaN ausserhalb des Gitters."""
    pts = np.atleast_2d(np.asarray(pts_xyz, dtype=float))
    if len(pts) == 0:
        return np.empty(0, dtype=np.float32)
    vol, offset = _crop_for_sampling(dose, pts, order)
    return _sample(vol, offset, dose, pts, order)


def isodose_mask(dose_fine: np.ndarray, level_gy: float) -> np.ndarray:
    """bool-Maske ``dose_fine >= level`` (NaN -> False)."""
    with np.errstate(invalid="ignore"):
        return np.nan_to_num(dose_fine, nan=-np.inf) >= level_gy


def label_components(mask: np.ndarray) -> tuple:
    """``scipy.ndimage.label`` mit 6er-Nachbarschaft -> (labels int32, n)."""
    labels, n = ndimage.label(mask)
    return labels.astype(np.int32), int(n)


def scope_mask_to_structure(iso_mask: np.ndarray, struct_mask: np.ndarray,
                            labels: Optional[np.ndarray] = None,
                            grid: Optional[FineGrid] = None) -> tuple:
    """
    Beschraenkt eine Isodosenmaske auf die Zusammenhangskomponenten, die die
    Struktur ueberlappen.  Ohne Ueberlapp: naechste Komponente (Schwerpunkt-
    abstand) mit ``fallback=True``.  Liefert ``(mask, info)``.
    """
    if labels is None:
        labels, _ = label_components(iso_mask)
    n_total = int(labels.max())
    info = {"n_components_total": n_total, "components_used": [], "fallback": False,
            "fallback_distance_mm": None}
    if n_total == 0:
        return np.zeros_like(iso_mask), info
    ids = np.unique(labels[struct_mask & iso_mask])
    ids = [int(i) for i in ids if i != 0]
    if not ids:
        if not struct_mask.any():
            return np.zeros_like(iso_mask), info
        sc = np.array(np.nonzero(struct_mask), dtype=float).mean(axis=1)
        best, best_d = None, np.inf
        for cid in range(1, n_total + 1):
            cc = np.array(np.nonzero(labels == cid), dtype=float).mean(axis=1)
            d = cc - sc
            if grid is not None:
                d = d * np.array([grid.dz, grid.res_xy, grid.res_xy])
            dist = float(np.linalg.norm(d))
            if dist < best_d:
                best, best_d = cid, dist
        ids = [best]
        info["fallback"] = True
        info["fallback_distance_mm"] = best_d
    info["components_used"] = ids
    return np.isin(labels, ids), info


def mask_volume_cm3(weighted: np.ndarray, grid: FineGrid) -> float:
    """``sum(weighted) * Voxelvolumen / 1000`` (bool- oder Gewichtsarray)."""
    return float(np.sum(weighted, dtype=np.float64)) * grid.voxel_volume_mm3 / 1000.0


def intersection_volume_cm3(weighted_struct: np.ndarray, iso_mask: np.ndarray,
                            grid: FineGrid) -> float:
    return float(np.sum(weighted_struct * iso_mask, dtype=np.float64)) * grid.voxel_volume_mm3 / 1000.0


def dice(a: np.ndarray, b: np.ndarray) -> float:
    """Dice-Koeffizient 2|a&b| / (|a|+|b|) fuer bool- oder Gewichtsarrays (0 wenn leer)."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    denom = a.sum() + b.sum()
    if denom <= 0:
        return 0.0
    return float(2.0 * np.sum(np.minimum(a, b)) / denom)


# ---------------------------------------------------------------------------
# 5. Gewichtete DVH-Statistik
# ---------------------------------------------------------------------------

def dose_at_volume_fraction(sorted_desc: np.ndarray, cum_w: np.ndarray, frac: float) -> float:
    """D_x: Dosis, die mindestens ``frac`` (0..1) des gewichteten Volumens erhaelt."""
    if len(sorted_desc) == 0:
        return float("nan")
    if frac <= cum_w[0]:
        return float(sorted_desc[0])
    if frac >= cum_w[-1]:
        return float(sorted_desc[-1])
    return float(np.interp(frac, cum_w, sorted_desc))


def volume_fraction_at_dose(dose_vals: np.ndarray, weights: np.ndarray, level_gy: float) -> float:
    """V_D: gewichteter Volumenanteil (0..1) mit Dosis >= level."""
    tot = float(np.sum(weights))
    if tot <= 0:
        return float("nan")
    return float(np.sum(weights[dose_vals >= level_gy]) / tot)


def weighted_dose_statistics(dose_vals: np.ndarray, weights: np.ndarray, rx_gy: float) -> dict:
    """
    NaN-Werte werden verworfen (Anteil in ``outside_fraction``).  Liefert
    ``dmin, dmax, dmean, d2, d50, d95, d98`` (Gy), ``v95_pct, v100_pct``
    (% des gewichteten Volumens mit D >= 0.95 Rx bzw. >= Rx), ``n_samples``.
    """
    d = np.asarray(dose_vals, dtype=np.float64).ravel()
    w = np.asarray(weights, dtype=np.float64).ravel()
    valid = ~np.isnan(d)
    w_tot = float(w.sum())
    outside = float(w[~valid].sum() / w_tot) if w_tot > 0 else 0.0
    d, w = d[valid], w[valid]
    nan = float("nan")
    if len(d) == 0 or w.sum() <= 0:
        empty = {k: nan for k in ("dmin", "dmax", "dmean", "d2", "d50", "d95", "d98",
                                  "v95_pct", "v100_pct")}
        empty.update({"n_samples": 0, "outside_fraction": outside})
        return empty
    order = np.argsort(-d, kind="stable")
    ds, ws = d[order], w[order]
    cw = np.cumsum(ws) / ws.sum()
    return {
        "dmin": float(ds[-1]),
        "dmax": float(ds[0]),
        "dmean": float(np.sum(ds * ws) / ws.sum()),
        "d2": dose_at_volume_fraction(ds, cw, 0.02),
        "d50": dose_at_volume_fraction(ds, cw, 0.50),
        "d95": dose_at_volume_fraction(ds, cw, 0.95),
        "d98": dose_at_volume_fraction(ds, cw, 0.98),
        "v95_pct": 100.0 * volume_fraction_at_dose(ds, ws, 0.95 * rx_gy),
        "v100_pct": 100.0 * volume_fraction_at_dose(ds, ws, rx_gy),
        "n_samples": int(len(ds)),
        "outside_fraction": outside,
    }


# ---------------------------------------------------------------------------
# 6. Maske -> Konturen (Marching Squares) und Roundtrip-Check
# ---------------------------------------------------------------------------

def _trace_rings(values: np.ndarray, level: float, grid: FineGrid, k: int,
                 simplify_mm: float, min_area_mm2: float, snap_to_centres: bool = False,
                 pad_value: float = 0.0) -> list:
    """
    Marching Squares (``skimage.measure.find_contours``) einer Ebene -> Liste
    geschlossener (N,3)-Ringe in mm.  Die Ebene wird mit einem Rand aus
    ``pad_value`` gepolstert, damit randberuehrende Ringe schliessen.
    Aussenringe haben positive, Loecher negative Flaeche (Orientierung wird
    nicht veraendert).  ``snap_to_centres`` rundet die Vertices auf
    Zellmittelpunkte und entfernt Doppelpunkte.
    """
    from skimage.measure import find_contours
    from shapely.geometry import LinearRing

    res = grid.res_xy
    padded = np.pad(values.astype(np.float64), 1, constant_values=pad_value)
    rings = []
    for rc in find_contours(padded, level):
        if len(rc) > 1 and np.allclose(rc[0], rc[-1]):
            rc = rc[:-1]
        if snap_to_centres:
            rc = np.round(rc)
            keep = np.ones(len(rc), dtype=bool)
            keep[1:] = np.any(rc[1:] != rc[:-1], axis=1)
            rc = rc[keep]
            if len(rc) > 1 and np.all(rc[0] == rc[-1]):
                rc = rc[:-1]
        if len(rc) < 3:
            continue
        x = grid.gx[0] + (rc[:, 1] - 1.0) * res
        y = grid.gy[0] + (rc[:, 0] - 1.0) * res
        xy = np.column_stack([x, y])
        signed = 0.5 * float(np.sum(xy[:, 0] * np.roll(xy[:, 1], -1)
                                    - np.roll(xy[:, 0], -1) * xy[:, 1]))
        if abs(signed) < min_area_mm2:
            continue
        if simplify_mm > 0:
            try:
                simp = LinearRing(xy).simplify(simplify_mm, preserve_topology=True)
                coords = np.asarray(simp.coords)[:-1]
                if len(coords) >= 3:
                    xy = coords
            except Exception:
                pass
        rings.append(np.column_stack([xy, np.full(len(xy), float(grid.gz[k]))]))
    return rings


def mask_to_contours(mask: np.ndarray, grid: FineGrid, simplify_mm: float = 0.1,
                     min_area_mm2: float = 0.05, vertex_mode: str = "edge") -> dict:
    """
    ``{z_mm: [(N,3) Ring, ...]}`` aus einer bool-Maske: pro Ebene Marching
    Squares bei 0.5 (Vertices auf den Kantenmitten zwischen Innen- und
    Aussenvoxeln).  ``vertex_mode="center"`` legt die Vertices stattdessen auf
    die Mittelpunkte der inneren Randvoxel.  Vereinfachung per Douglas-Peucker
    (``simplify_mm`` muss < res/2 bleiben, damit kein Voxelmittelpunkt die
    Seite wechselt).  Ringe mit < 3 Punkten oder |Flaeche| < ``min_area_mm2``
    werden verworfen.
    """
    if vertex_mode not in ("edge", "center"):
        raise ValueError(f"vertex_mode muss 'edge' oder 'center' sein, nicht {vertex_mode!r}")
    level = 0.5 if vertex_mode == "edge" else 1.0 - 1e-3
    out = {}
    for k in range(len(grid.gz)):
        plane = mask[k]
        if not plane.any():
            continue
        rings = _trace_rings(plane.astype(np.float64), level, grid, k, simplify_mm,
                             min_area_mm2, snap_to_centres=(vertex_mode == "center"))
        if rings:
            out[float(grid.gz[k])] = rings
    return out


def field_to_contours(field: np.ndarray, level: float, grid: FineGrid,
                      simplify_mm: float = 0.0, min_area_mm2: float = 0.05) -> dict:
    """
    Isolinie eines kontinuierlichen Feldes (z.B. der Dosis) bei ``level``:
    Vertices liegen auf den Gitterlinien -- eine Koordinate exakt auf dem
    Zellmittelpunkt-Gitter, die andere linear interpoliert.  Das ist die
    Vertex-Konvention der Eclipse-Konturen ("High"/"Default" Resolution auf dem
    CT-Pixelgitter).  NaN (ausserhalb des Dosisgitters) zaehlt als unter Level.
    """
    out = {}
    low = -1.0e30
    for k in range(len(grid.gz)):
        plane = np.nan_to_num(field[k].astype(np.float64), nan=low)
        if not (plane >= level).any():
            continue
        rings = _trace_rings(plane, level, grid, k, simplify_mm, min_area_mm2, pad_value=low)
        if rings:
            out[float(grid.gz[k])] = rings
    return out


def contours_flat(contours_by_z: dict) -> list:
    """``{z: [rings]}`` -> flache Liste von (N,3)-Arrays (z aufsteigend)."""
    return [r for z in sorted(contours_by_z) for r in contours_by_z[z]]


def contours_roundtrip_dice(mask: np.ndarray, grid: FineGrid, **kw) -> float:
    """Dice(mask, analyzer.rasterize_contours(mask_to_contours(mask)))."""
    flat = contours_flat(mask_to_contours(mask, grid, **kw))
    if not flat:
        return 0.0 if mask.any() else 1.0
    back = ana.rasterize_contours(flat, grid.gx, grid.gy, grid.gz)
    return dice(mask, back)


# ---------------------------------------------------------------------------
# 7. Mini-CLI: Gitterinfo
# ---------------------------------------------------------------------------

def main(argv: Optional[list] = None) -> int:
    import argparse
    p = argparse.ArgumentParser(description="RTDOSE-Gitterinfo (Dosisindex-Bibliothek)")
    p.add_argument("file", help="Pfad zur RTDOSE-Datei")
    args = p.parse_args(argv)
    try:
        ds = load_rtdose(args.file)
        dose = dose_grid_from_dataset(ds, args.file)
    except (FileNotFoundError, ValueError, KeyError) as e:
        print(f"Fehler: {e}", file=sys.stderr)
        return 2
    print(f"RTDOSE      : {args.file}")
    print(f"  Typ       : {dose.summation_type} / {dose.dose_type} / {dose.units}")
    print(f"  Gitter    : {dose.shape[0]} x {dose.shape[1]} x {dose.shape[2]} Voxel "
          f"(k x j x i) @ {dose.spacing[0]:.3f}/{dose.spacing[1]:.3f}/{dose.spacing[2]:.3f} mm")
    print(f"  Ursprung  : ({dose.origin[0]:.3f}, {dose.origin[1]:.3f}, {dose.origin[2]:.3f}) mm, "
          f"GFOV {dose.gfov_mode}")
    print(f"  Dmax      : {dose.dmax:.4f} Gy")
    print(f"  FoR       : {dose.frame_of_reference_uid}")
    print(f"  Ref. RP   : {dose.referenced_plan_uid}")
    print(f"  Ref. RS   : {dose.referenced_rtstruct_uid}")
    for w in dose.warnings:
        print(f"  Hinweis: {w}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
