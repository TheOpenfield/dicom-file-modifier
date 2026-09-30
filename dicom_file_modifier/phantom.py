"""
phantom.py - Analytische Geometrie- und Dosismodelle fuer synthetische Faelle.

Grundlage des Demo-/Testfalls (``demo.py``): Konturen von Kugeln, Zylindern und
Ringen als regelmaessige n-Ecke, ein glattes radiales Dosisfeld je Ziel
(Hill-Profil wie im Self-Test von ``dose_indices``) und Erwartungswerte in
geschlossener Form (Scheibenstapel, Slab-Modell).  Keine Datei-I/O.

Konventionen wie im restlichen Paket: Patientenkoordinaten LPS in mm,
Volumina in cm3, Dosis in Gy.  Konsolenausgaben gibt es hier keine.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Optional, Sequence

import numpy as np

from .analyzer import _synthetic_circle, rasterize_contours


# ---------------------------------------------------------------------------
# 1. Konturen
# ---------------------------------------------------------------------------

def circle(z: float, r: float, cx: float = 0.0, cy: float = 0.0, n: int = 360) -> np.ndarray:
    """Regelmaessiges n-Eck (gegen den Uhrzeigersinn) als (n,3)-Kontur bei Hoehe z."""
    return _synthetic_circle(z, r, cx, cy, n)


def sphere_contours(center, radius: float, z_planes: Sequence[float], n: int = 360) -> list:
    """Kreisschnitte einer Kugel auf allen Ebenen mit |z - cz| < R."""
    cx, cy, cz = (float(v) for v in center)
    out = []
    for z in z_planes:
        dz = float(z) - cz
        if abs(dz) < radius:
            out.append(circle(float(z), math.sqrt(radius * radius - dz * dz), cx, cy, n))
    return out


def cylinder_contours(cx: float, cy: float, radius: float, z_planes: Sequence[float],
                      z_min: float, z_max: float, n: int = 360) -> list:
    """Achsparalleler Zylinder: gleicher Kreis auf allen Ebenen in [z_min, z_max]."""
    return [circle(float(z), radius, cx, cy, n) for z in z_planes if z_min <= float(z) <= z_max]


def ring_contours(cx: float, cy: float, r_out: float, r_in: float, z_planes: Sequence[float],
                  z_min: float, z_max: float, n: int = 360) -> list:
    """Ring als zwei verschachtelte Konturen je Ebene (XOR -> Loch, Eclipse-Konvention)."""
    out = []
    for z in z_planes:
        if z_min <= float(z) <= z_max:
            out.append(circle(float(z), r_out, cx, cy, n))
            out.append(circle(float(z), r_in, cx, cy, n))
    return out


def polygon_area_factor(n: int) -> float:
    """Flaeche eines regelmaessigen n-Ecks relativ zum Umkreis: (n / 2pi) sin(2pi / n)."""
    return n / (2.0 * math.pi) * math.sin(2.0 * math.pi / n)


def lens_area(a: float, b: float, d: float) -> float:
    """Schnittflaeche zweier Kreise (Radien a, b, Mittelpunktsabstand d)."""
    if a <= 0 or b <= 0 or d >= a + b:
        return 0.0
    if d <= abs(a - b):
        return math.pi * min(a, b) ** 2
    t1 = a * a * math.acos((d * d + a * a - b * b) / (2 * d * a))
    t2 = b * b * math.acos((d * d + b * b - a * a) / (2 * d * b))
    t3 = 0.5 * math.sqrt((-d + a + b) * (d + a - b) * (d - a + b) * (d + a + b))
    return t1 + t2 - t3


def _trapezoid(y: np.ndarray, x: np.ndarray) -> float:
    """Trapezregel ohne ``np.trapz`` (in numpy 2 veraltet, ``np.trapezoid`` fehlt in 1.x)."""
    y, x = np.asarray(y, float), np.asarray(x, float)
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))


# ---------------------------------------------------------------------------
# 2. Dosisfeld
# ---------------------------------------------------------------------------

@dataclass
class HillField:
    """
    Glattes radiales Dosisfeld ``D(r) = Dmax / (1 + (r / r0)^p)`` um ``center``.
    Bewusst C-unendlich: ein Knick am Rx-Level wuerde die lineare Interpolation
    systematisch verzerren (siehe Self-Test in ``dose_indices``).
    """
    center: np.ndarray
    dmax: float
    r0: float
    p: float = 6.0

    @classmethod
    def for_target(cls, target_center, target_radius: float, rx: float,
                   dmax_factor: float = 1.25, r100_margin: float = 0.5, p: float = 6.0,
                   offset=(2.0, 0.0, 0.0)) -> "HillField":
        """Feld mit Rx-Isodose bei ``target_radius + r100_margin`` um ``center + offset``."""
        dmax = dmax_factor * rx
        r100 = target_radius + r100_margin
        r0 = r100 / (dmax / rx - 1.0) ** (1.0 / p)
        center = np.asarray(target_center, float) + np.asarray(offset, float)
        return cls(center=center, dmax=dmax, r0=r0, p=p)

    def dose_at(self, pts) -> np.ndarray:
        r = np.linalg.norm(np.atleast_2d(np.asarray(pts, float)) - self.center, axis=1)
        return self.dmax / (1.0 + (r / self.r0) ** self.p)

    def level_radius(self, level: float) -> float:
        """Radius der Isodosenkugel ``D >= level`` (0, wenn ``level >= Dmax``)."""
        if level >= self.dmax:
            return 0.0
        return self.r0 * (self.dmax / level - 1.0) ** (1.0 / self.p)

    def isodose_volume_cm3(self, level: float, dose_planes: Sequence[float], dz: float) -> float:
        """Isodosenvolumen als Scheibenstapel auf den Dosisebenen (Slab-Modell)."""
        rl = self.level_radius(level)
        cz = float(self.center[2])
        return sum(math.pi * max(rl * rl - (float(z) - cz) ** 2, 0.0) for z in dose_planes) * dz / 1000.0


def summed_dose(fields: Sequence[HillField]) -> Callable[[np.ndarray], np.ndarray]:
    """Dosisfunktion (N,3) -> (N,) als Summe der Felder."""
    def _f(pts):
        pts = np.atleast_2d(np.asarray(pts, float))
        out = np.zeros(len(pts))
        for fld in fields:
            out += fld.dose_at(pts)
        return out
    return _f


# ---------------------------------------------------------------------------
# 3. Erwartungswerte fuer ein kugelfoermiges Ziel (geschlossene Form)
# ---------------------------------------------------------------------------

@dataclass
class SphereTarget:
    """
    Kugelziel mit eigenem Hill-Feld.  Erwartungswerte als Scheibenstapel auf den
    Konturebenen (Slab-Modell; ``eclipse`` halbiert die Endscheiben), Polygon-
    korrektur fuer das n-Eck wie im Self-Test.  Uebersprechen anderer Felder
    wird vernachlaessigt (bei > 40 mm Abstand < 0.02 Gy).
    """
    name: str
    center: np.ndarray
    radius: float
    rx: float
    field: HillField
    contour_planes: np.ndarray
    dz: float
    n_polygon: int = 360

    def _disc_radius(self, z: float) -> float:
        return math.sqrt(max(self.radius ** 2 - (z - float(self.center[2])) ** 2, 0.0))

    def _weights(self, model: str) -> np.ndarray:
        w = np.ones(len(self.contour_planes))
        if model == "eclipse" and len(w):
            w[0] = w[-1] = 0.5
        return w

    def target_volume(self, model: str = "slab") -> float:
        pf, w = polygon_area_factor(self.n_polygon), self._weights(model)
        return sum(wk * math.pi * self._disc_radius(float(z)) ** 2 * pf
                   for wk, z in zip(w, self.contour_planes)) * self.dz / 1000.0

    def cumulative_volume(self, level: float, model: str = "slab") -> float:
        """V(D >= level) des Ziels in cm3 (Kreis-Kreis-Linsen je Scheibe)."""
        if level <= 0:
            return self.target_volume(model)
        pf, w = polygon_area_factor(self.n_polygon), self._weights(model)
        d = float(np.linalg.norm(self.center[:2] - self.field.center[:2]))
        rl = self.field.level_radius(level)
        fz = float(self.field.center[2])
        tot = 0.0
        for wk, z in zip(w, self.contour_planes):
            a = self._disc_radius(float(z)) * math.sqrt(pf)
            b = math.sqrt(max(rl * rl - (float(z) - fz) ** 2, 0.0))
            tot += wk * lens_area(a, b, d)
        return tot * self.dz / 1000.0

    def expected(self, dose_planes: Sequence[float], other_fields: Sequence[HillField] = (),
                 model: str = "slab", rx: Optional[float] = None) -> dict:
        """
        Kennzahlen fuer Verschreibung ``rx`` (Default: ``self.rx``).  PIV/PIV50
        einmal nur fuer das eigene Feld (``component``) und einmal inklusive der
        ``other_fields`` (``global``; Felder disjunkt angenommen).
        """
        rx = float(self.rx if rx is None else rx)
        tv = self.target_volume(model)
        inter = self.cumulative_volume(rx, model)
        piv_c = self.field.isodose_volume_cm3(rx, dose_planes, self.dz)
        piv50_c = self.field.isodose_volume_cm3(0.5 * rx, dose_planes, self.dz)
        piv_g = piv_c + sum(f.isodose_volume_cm3(rx, dose_planes, self.dz) for f in other_fields)
        piv50_g = piv50_c + sum(f.isodose_volume_cm3(0.5 * rx, dose_planes, self.dz)
                                for f in other_fields)

        def d_at(frac: float) -> float:
            lo, hi = 0.01, self.field.dmax
            for _ in range(80):
                mid = 0.5 * (lo + hi)
                if self.cumulative_volume(mid, model) / tv >= frac:
                    lo = mid
                else:
                    hi = mid
            return 0.5 * (lo + hi)

        levels = np.linspace(0.0, self.field.dmax, 801)
        vcum = np.array([self.cumulative_volume(float(L), model) if L > 0 else tv for L in levels])
        dmean = _trapezoid(vcum, levels) / tv
        d2, d50, d98, d95 = d_at(0.02), d_at(0.50), d_at(0.98), d_at(0.95)

        def r_eq_cm(v_cm3: float) -> float:
            return (3.0 * v_cm3 / (4.0 * math.pi)) ** (1.0 / 3.0)

        return {
            "rx_gy": rx, "volume_model": model,
            "tv_cm3": tv, "tv_piv_cm3": inter,
            "piv_component_cm3": piv_c, "piv50_component_cm3": piv50_c,
            "piv_global_cm3": piv_g, "piv50_global_cm3": piv50_g,
            "ci_paddick_component": inter * inter / (tv * piv_c),
            "ci_paddick_global": inter * inter / (tv * piv_g),
            "gi_component": piv50_c / piv_c, "gi_global": piv50_g / piv_g,
            "gm_cm_component": r_eq_cm(piv50_c) - r_eq_cm(piv_c),
            "coverage": inter / tv,
            "d2_gy": d2, "d50_gy": d50, "d95_gy": d95, "d98_gy": d98, "dmean_gy": dmean,
            "hi_icru83": (d2 - d98) / d50,
            "v100_pct": 100.0 * inter / tv,
            "v95_pct": 100.0 * self.cumulative_volume(0.95 * rx, model) / tv,
        }


# ---------------------------------------------------------------------------
# 4. Numerisches DVH (fuer die Eclipse-artige DVHSequence der RTDOSE)
# ---------------------------------------------------------------------------

def cumulative_dvh(contours: list, dose_fn: Callable[[np.ndarray], np.ndarray],
                   res_xy: float, bin_width: float, dz: float,
                   box: Optional[tuple] = None, total_volume: Optional[float] = None) -> dict:
    """
    Kumulatives DVH einer ROI: XOR-Rasterung der Konturen (``rasterize_contours``)
    auf Voxelmittelpunkte im Abstand ``res_xy``, z auf den Konturebenen (volle
    Schichtdicke ``dz``, Slab-Modell), Dosis aus ``dose_fn``.

    ``box = (x0, x1, y0, y1)`` begrenzt die Abtastung auf das Dosisgitter.  Ragt
    die ROI hinaus, bekommt der Teil ausserhalb Dosis 0 (wie ein TPS-DVH) und
    zaehlt nur im Bin 0, dessen Volumen dann ``total_volume`` ist (z.B.
    planimetrisch exakt).  Liegt sie ganz innerhalb, zaehlen nur die Voxel.

    Rueckgabe ``{'edges', 'volumes', 'dmin', 'dmax', 'dmean', 'total'}``:
    ``volumes[i]`` = Volumen mit D >= ``edges[i]`` (linke Binkante wie Eclipse).
    """
    pts = np.vstack([c for c in contours if len(c) >= 3])
    zc = np.unique(np.round(pts[:, 2], 6))
    lo, hi = pts[:, :2].min(axis=0) - res_xy, pts[:, :2].max(axis=0) + res_xy
    clipped = box is not None and bool(
        pts[:, 0].min() < box[0] or pts[:, 0].max() > box[1]
        or pts[:, 1].min() < box[2] or pts[:, 1].max() > box[3])
    if clipped:
        lo = np.maximum(lo, [box[0], box[2]])
        hi = np.minimum(hi, [box[1], box[3]])
    xc = np.arange(lo[0] + 0.5 * res_xy, hi[0], res_xy)
    yc = np.arange(lo[1] + 0.5 * res_xy, hi[1], res_xy)
    mask = rasterize_contours(contours, xc, yc, zc)
    kk, jj, ii = np.nonzero(mask)
    vox = np.column_stack([xc[ii], yc[jj], zc[kk]])
    dose = dose_fn(vox)
    v_vox = res_xy * res_xy * dz / 1000.0
    sampled = float(len(dose) * v_vox)
    total = float(total_volume) if (clipped and total_volume is not None) else sampled
    outside = max(total - sampled, 0.0)
    n_bins = int(math.floor(float(dose.max()) / bin_width)) + 2
    edges = np.arange(n_bins) * bin_width
    srt = np.sort(dose)
    # Anzahl Voxel mit D >= edge: len - Anzahl mit D < edge
    counts = len(srt) - np.searchsorted(srt, edges, side="left")
    volumes = counts * v_vox
    volumes[0] = total
    dmean = (float(dose.sum()) * v_vox) / total if total > 0 else 0.0
    return {
        "edges": edges, "volumes": volumes,
        "dmin": 0.0 if outside > 0 else float(dose.min()), "dmax": float(dose.max()),
        "dmean": dmean, "total": total,
    }
