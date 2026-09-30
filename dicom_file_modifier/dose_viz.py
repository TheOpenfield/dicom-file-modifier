"""
dose_viz.py - Validierungsansicht fuer die Dosisindex-Berechnung (dose_indices).

Erzeugt aus den ``DoseIndexArtifacts`` und den exportierten ROI-Konturen
(``rtstruct_writer.build_roi_specs``) zwei Dateien im Ausgabeordner:

  validation.html    Offline-Seite mit drei unabhaengigen Plotly-Figuren
                     (3D: Isodosen-/Zielflaechen + exportierte Konturen;
                     axialer Schichtbrowser mit CT-Hintergrund, Dosis-Wash und
                     Konturen je Ebene (Slider); kumulatives DVH je Ziel mit
                     Eclipse-DVH-Overlay) und einer Indextabelle.  plotly.js
                     liegt einmal inline (~4.9 MB), keine Netzverbindung noetig.
  dose_overview.png  statisch (matplotlib): axial / koronal / sagittal durch den
                     Schwerpunkt des ersten Ziels + DVH-Panel.

Datenbudget: Dosis-Wash und CT werden je Achse auf <= VIZ_MAX_PX Stuetzstellen
gestrided (echte Stichproben, keine Mittelung), Flaechen (Marching Cubes) auf
<= VIZ_MESH_MAX_PX in-plane, die DVH-Kurve auf <= VIZ_DVH_MAX_PTS Punkte.
Konturen (exportierte Isodosen/Hilfs-ROIs, Original-Ziel) werden nie
ausgeduennt.  Isodosen-Flaechen entstehen aus dem Dosisfeld (Marching Cubes
beim Level), die Zielflaeche aus der Zielmaske (Level 0.5).

Bekannte Einschraenkung: der Slider setzt die Sichtbarkeit aller Traces neu;
Legenden-Toggles gelten daher nur bis zum naechsten Schichtwechsel.

Bibliothek ohne CLI (Aufruf aus ``dose_indices.run_dose_indices``).
Konsolenausgabe ASCII, Dateien UTF-8.
"""

from __future__ import annotations

import html as _html
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pydicom

from . import dose as dm
from . import modifier as mod
from .dose_constants import LEVEL_COLORS

VIZ_MAX_PX = 128          # Wash/CT-Raster je Achse (HTML)
VIZ_MESH_MAX_PX = 96      # Marching-Cubes-Eingabe in-plane
VIZ_DVH_MAX_PTS = 2000    # Punkte je DVH-Kurve
CT_WINDOW = (40.0, 400.0)  # (Level, Width) in HU: Weichteilfenster
CT_MARGIN_MM = 15.0       # CT-Ausschnitt um die Feingitter-BBox
WASH_MIN_FRAC = 0.10      # Dosis-Wash unter 10 % Rx transparent
HTML_NAME = "validation.html"
PNG_NAME = "dose_overview.png"

_CSS = """
body{background:#14141e;color:#e8e8e8;font-family:'Segoe UI',Arial,sans-serif;margin:0;padding:16px 24px}
h1{font-size:20px;margin:0 0 4px}
h2{font-size:16px;margin:26px 0 8px;color:#cfd8dc}
.meta{color:#aab;font-size:13px;margin-bottom:12px;line-height:1.5}
table{border-collapse:collapse;font-size:13px;margin:8px 0}
th,td{border:1px solid #333;padding:3px 10px;text-align:right;white-space:nowrap}
th{background:#22223a}
td:first-child,th:first-child{text-align:left}
.flag{color:#ff6b6b;font-weight:bold}
.src{color:#889;font-size:11px}
.note{color:#aab;font-size:12px}
"""


# ---------------------------------------------------------------------------
# 1. Datenobjekte und Helfer
# ---------------------------------------------------------------------------

@dataclass
class CtWindow:
    """CT-Ausschnitt (HU) um das Feingitter: ``hu (nz, ny, nx)``, Achsen in mm."""
    hu: np.ndarray
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray

    @property
    def dx(self) -> float:
        return float(self.x[1] - self.x[0]) if len(self.x) > 1 else 1.0

    @property
    def dy(self) -> float:
        return float(self.y[1] - self.y[0]) if len(self.y) > 1 else 1.0

    @property
    def dz(self) -> float:
        return float(np.median(np.diff(self.z))) if len(self.z) > 1 else 1.0

    def nearest_k(self, z: float) -> Optional[int]:
        """Index der naechsten CT-Schicht (innerhalb dz/2), sonst None."""
        if len(self.z) == 0:
            return None
        k = int(np.argmin(np.abs(self.z - z)))
        return k if abs(float(self.z[k]) - z) <= 0.5 * self.dz + 1e-6 else None


@dataclass
class DvhCurve:
    """Kumulatives DVH (Dosis aufsteigend, Volumen cm3) plus D98/D50/D2."""
    dose_gy: np.ndarray
    volume_cm3: np.ndarray
    total_cm3: float
    d98_gy: float
    d50_gy: float
    d2_gy: float
    n_samples: int


def _rgb(color) -> str:
    r, g, b = (int(v) for v in color[:3])
    return f"rgb({r},{g},{b})"


def _mpl_color(color) -> tuple:
    return tuple(min(max(int(v), 0), 255) / 255.0 for v in color[:3])


def _stride(n: int, cap: int) -> int:
    return max(1, int(math.ceil(n / float(cap))))


def _fmt(v, nd: int = 3) -> str:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return "n/a"
    if math.isnan(f) or math.isinf(f):
        return "n/a"
    return f"{f:.{nd}f}"


def _esc(s) -> str:
    return _html.escape(str(s), quote=True)


def _mask_centroid(mask: np.ndarray, grid: dm.FineGrid) -> tuple:
    """((x, y, z) mm, (k, j, i)) des Maskenschwerpunkts (Gittermitte, wenn leer)."""
    idx = np.nonzero(mask)
    if len(idx[0]) == 0:
        k, j, i = len(grid.gz) // 2, len(grid.gy) // 2, len(grid.gx) // 2
        return (float(grid.gx[i]), float(grid.gy[j]), float(grid.gz[k])), (k, j, i)
    x, y, z = float(grid.gx[idx[2]].mean()), float(grid.gy[idx[1]].mean()), float(grid.gz[idx[0]].mean())
    k = int(np.argmin(np.abs(grid.gz - z)))
    j = int(np.argmin(np.abs(grid.gy - y)))
    i = int(np.argmin(np.abs(grid.gx - x)))
    return (x, y, z), (k, j, i)


def _rings_by_plane(contours: list, grid: dm.FineGrid) -> dict:
    """``{k: [(N,3) Ring, ...]}`` der Original-Konturen je Gitterebene."""
    out = {}
    for c in contours:
        c = np.asarray(c, dtype=float)
        if len(c) == 0:
            continue
        k = grid.plane_index(float(c[0, 2]))
        if k is not None:
            out.setdefault(k, []).append(c)
    return out


def _polyline_xyz(rings: list) -> tuple:
    """Geschlossene Ringe, durch None getrennt -> (xs, ys, zs) fuer Plotly."""
    xs, ys, zs = [], [], []
    for r in rings:
        r = np.asarray(r, dtype=float)
        if len(r) == 0:
            continue
        closed = np.round(np.vstack([r, r[:1]]), 3)
        xs.extend(closed[:, 0].tolist())
        ys.extend(closed[:, 1].tolist())
        zs.extend(closed[:, 2].tolist())
        xs.append(None)
        ys.append(None)
        zs.append(None)
    return xs, ys, zs


def _dose_anchors(rx: float, zmin: float, zmax: float) -> list:
    """``[(pos 0..1, (r, g, b)), ...]`` fuer Plotly-Colorscale und matplotlib-Colormap."""
    pts = [(pct / 100.0 * rx, LEVEL_COLORS[pct]) for pct in sorted(LEVEL_COLORS)]
    pts.append((max(zmax, 1.05 * rx), (139, 0, 0)))
    span = max(zmax - zmin, 1e-6)
    out = []
    for gy, col in pts:
        pos = min(max((gy - zmin) / span, 0.0), 1.0)
        if out and pos <= out[-1][0]:
            continue
        out.append((pos, col))
    if out[0][0] > 0.0:
        out.insert(0, (0.0, out[0][1]))
    if out[-1][0] < 1.0:
        out.append((1.0, out[-1][1]))
    return out


def _dose_range(art) -> tuple:
    rx = float(art.rx_gy)
    zmin = WASH_MIN_FRAC * rx
    dmax = float(np.nanmax(art.dose_fine)) if np.isfinite(np.nanmax(art.dose_fine)) else 0.0
    zmax = max(dmax, 1.1 * rx)
    return zmin, zmax


# ---------------------------------------------------------------------------
# 2. CT-Ausschnitt laden (nur Schichten im z-Bereich des Feingitters)
# ---------------------------------------------------------------------------

def load_ct_window(ct_index: Optional[dict], grid: dm.FineGrid, margin_mm: float = CT_MARGIN_MM,
                   notes: Optional[list] = None) -> Optional[CtWindow]:
    """
    Laedt die CT-Schichten mit z in ``[gz[0]-dz-margin, gz[-1]+dz+margin]``
    (Pfade aus ``ct_index['z_to_path']``), wandelt in HU (RescaleSlope/
    Intercept) und schneidet auf die Feingitter-BBox +- ``margin_mm`` zu.
    None (mit Hinweis in ``notes``), wenn keine Pfade vorliegen, das CT nicht
    axial ist oder keine Schicht im Bereich liegt.
    """
    def note(msg):
        if notes is not None:
            notes.append(msg)

    if not ct_index or not ct_index.get("z_to_path"):
        return None
    zlo = float(grid.gz[0]) - grid.dz - margin_mm
    zhi = float(grid.gz[-1]) + grid.dz + margin_mm
    zs = [z for z in sorted(ct_index["z_to_path"]) if zlo <= z <= zhi]
    if not zs:
        note("CT-Hintergrund: keine CT-Schicht im z-Bereich des Feingitters.")
        return None
    (xlo, ylo, _), (xhi, yhi, _) = grid.bbox
    dr, dc = (float(v) for v in ct_index["pixel_spacing"])   # Zeilenabstand (y), Spaltenabstand (x)
    x0, y0 = (float(v) for v in ct_index["ipp_xy"])
    hu_list, z_list = [], []
    i0 = i1 = j0 = j1 = None
    for z in zs:
        ds = pydicom.dcmread(ct_index["z_to_path"][z])
        iop = np.asarray(ds.ImageOrientationPatient, dtype=float)
        if not np.allclose(iop, [1, 0, 0, 0, 1, 0], atol=1e-3):
            note("CT-Hintergrund: CT nicht axial (IOP != 1,0,0,0,1,0); kein CT-Hintergrund.")
            return None
        arr = ds.pixel_array.astype(np.float32)
        arr = arr * float(getattr(ds, "RescaleSlope", 1.0)) + float(getattr(ds, "RescaleIntercept", 0.0))
        if i0 is None:
            ny, nx = arr.shape
            i0 = max(0, int(math.floor((xlo - margin_mm - x0) / dc)))
            i1 = min(nx, int(math.ceil((xhi + margin_mm - x0) / dc)) + 1)
            j0 = max(0, int(math.floor((ylo - margin_mm - y0) / dr)))
            j1 = min(ny, int(math.ceil((yhi + margin_mm - y0) / dr)) + 1)
            if i1 <= i0 or j1 <= j0:
                note("CT-Hintergrund: Feingitter liegt ausserhalb des CT-Bildbereichs.")
                return None
        hu_list.append(arr[j0:j1, i0:i1])
        z_list.append(float(z))
    return CtWindow(hu=np.stack(hu_list).astype(np.float32),
                    x=x0 + np.arange(i0, i1) * dc, y=y0 + np.arange(j0, j1) * dr,
                    z=np.asarray(z_list, dtype=float))


# ---------------------------------------------------------------------------
# 3. DVH-Kurve
# ---------------------------------------------------------------------------

def build_dvh_curve(dose_samples: np.ndarray, sample_weights: np.ndarray,
                    voxel_volume_mm3: float, max_pts: int = VIZ_DVH_MAX_PTS) -> DvhCurve:
    """
    Kumulatives DVH aus den gewichteten Feingitter-Stichproben eines Ziels
    (NaN verworfen, Sortierung/Cumsum wie ``dose.weighted_dose_statistics``,
    D98/D50/D2 identisch).  Kurve auf ``max_pts`` Punkte gestrided (Endpunkte
    bleiben) und um den Plateaupunkt (0 Gy, Gesamtvolumen) ergaenzt.
    """
    d = np.asarray(dose_samples, dtype=np.float64).ravel()
    w = np.asarray(sample_weights, dtype=np.float64).ravel()
    valid = ~np.isnan(d)
    d, w = d[valid], w[valid]
    nan = float("nan")
    if len(d) == 0 or w.sum() <= 0:
        return DvhCurve(np.zeros(0), np.zeros(0), 0.0, nan, nan, nan, 0)
    order = np.argsort(-d, kind="stable")
    ds, ws = d[order], w[order]
    cw = np.cumsum(ws) / ws.sum()
    total = float(ws.sum() * voxel_volume_mm3 / 1000.0)
    d98 = dm.dose_at_volume_fraction(ds, cw, 0.98)
    d50 = dm.dose_at_volume_fraction(ds, cw, 0.50)
    d2 = dm.dose_at_volume_fraction(ds, cw, 0.02)
    n = len(ds)
    idx = np.unique(np.linspace(0, n - 1, max_pts).astype(int)) if n > max_pts else np.arange(n)
    dose_curve = np.concatenate([[0.0], ds[idx][::-1]])
    vol_curve = np.concatenate([[total], (cw[idx] * total)[::-1]])
    return DvhCurve(dose_gy=dose_curve, volume_cm3=vol_curve, total_cm3=total,
                    d98_gy=d98, d50_gy=d50, d2_gy=d2, n_samples=int(n))


# ---------------------------------------------------------------------------
# 4. Plotly-Figuren
# ---------------------------------------------------------------------------

_DARK = dict(paper_bgcolor="rgb(20, 20, 30)", font=dict(color="white"))
_GRID = dict(gridcolor="rgba(255,255,255,0.15)", color="white", zeroline=False)


def _volume_to_mesh(vol: np.ndarray, grid: dm.FineGrid, level: float, stride: int = 1) -> tuple:
    """
    Marching Cubes bei ``level`` auf einem in-plane gestrideten, um ein Voxel
    gepaddeten Volumen (Randflaechen schliessen); Vertices in Patienten-mm.
    ``(None, None)``, wenn der Level ausserhalb des Wertebereichs liegt.
    """
    v = np.asarray(vol, dtype=np.float32)[:, ::stride, ::stride]
    v = np.pad(v, 1, mode="constant", constant_values=0.0)
    if not (float(np.nanmin(v)) < level < float(np.nanmax(v))):
        return None, None
    A = np.zeros((4, 4))
    A[3, 3] = 1.0
    A[2, 0] = grid.dz
    A[1, 1] = grid.res_xy * stride
    A[0, 2] = grid.res_xy * stride
    A[:3, 3] = [grid.gx[0] - grid.res_xy * stride, grid.gy[0] - grid.res_xy * stride,
                grid.gz[0] - grid.dz]
    verts, faces = mod._extract_surface(v, A, threshold=level, downsample=1)
    if verts is None:
        return None, None
    return np.round(verts, 3), faces


def _mesh_trace(go, verts, faces, name: str, color, opacity: float, visible=True):
    return go.Mesh3d(
        x=verts[:, 0], y=verts[:, 1], z=verts[:, 2], i=faces[:, 0], j=faces[:, 1], k=faces[:, 2],
        color=_rgb(color), opacity=opacity, name=name, showlegend=True, visible=visible,
        flatshading=True, lighting=dict(diffuse=0.8, specular=0.2, roughness=0.6),
        hoverinfo="name",
    )


def _polyline_trace3d(go, rings: list, name: str, color, width: float = 3.0, visible=True):
    xs, ys, zs = _polyline_xyz(rings)
    return go.Scatter3d(x=xs, y=ys, z=zs, mode="lines", line=dict(color=_rgb(color), width=width),
                        name=name, visible=visible, hoverinfo="name")


def _fig_3d(go, art, specs: list):
    grid = art.grid
    stride = _stride(max(len(grid.gx), len(grid.gy)), VIZ_MESH_MAX_PX)
    fig = go.Figure()
    field = np.nan_to_num(np.asarray(art.dose_fine, dtype=np.float32), nan=0.0)
    for lv in art.levels.values():
        verts, faces = _volume_to_mesh(field, grid, lv.gy, stride)
        if verts is None:
            continue
        is_rx = lv.key == "100"
        fig.add_trace(_mesh_trace(go, verts, faces, f"Isodose {lv.label} ({lv.gy:.2f} Gy) Flaeche",
                                  lv.color, 0.30 if is_rx else 0.12, True if is_rx else "legendonly"))
    for tm in art.targets.values():
        verts, faces = _volume_to_mesh(tm.structure.mask.astype(np.float32), grid, 0.5, stride)
        if verts is not None:
            fig.add_trace(_mesh_trace(go, verts, faces, f"Ziel {tm.name} (Maske)", tm.color, 0.35))
        fig.add_trace(_polyline_trace3d(go, tm.contours, f"Ziel {tm.name} (Original-Konturen)",
                                        tm.color, 3.0))
    for sp in specs:
        if not sp.contours_by_z:
            continue
        rings = dm.contours_flat(sp.contours_by_z)
        iso = sp.kind == "isodose"
        fig.add_trace(_polyline_trace3d(go, rings, sp.name, sp.color, 4.0 if iso else 2.0,
                                        True if iso else "legendonly"))
    fig.update_layout(
        title=dict(text="3D: Isodosen-/Zielflaechen (Marching Cubes) und exportierte Konturen "
                        "(Legende: Ebenen einzeln zuschaltbar)", font=dict(size=14)),
        scene=dict(xaxis_title="X [mm] (Links)", yaxis_title="Y [mm] (Posterior)",
                   zaxis_title="Z [mm] (Superior)", aspectmode="data", bgcolor="rgb(20, 20, 30)",
                   xaxis=dict(gridcolor="rgba(255,255,255,0.15)", color="white"),
                   yaxis=dict(gridcolor="rgba(255,255,255,0.15)", color="white"),
                   zaxis=dict(gridcolor="rgba(255,255,255,0.15)", color="white")),
        legend=dict(x=0.01, y=0.99, bgcolor="rgba(0,0,0,0.4)", font=dict(color="white")),
        margin=dict(l=0, r=0, b=0, t=50), height=750, **_DARK,
    )
    return fig


def _fig_slices(go, art, specs: list, ct: Optional[CtWindow], k_init: int, target_rings: dict):
    grid, rx = art.grid, float(art.rx_gy)
    gx, gy, gz = grid.gx, grid.gy, grid.gz
    s = _stride(max(len(gx), len(gy)), VIZ_MAX_PX)
    zmin, zmax = _dose_range(art)
    scale = [[p, _rgb(c)] for p, c in _dose_anchors(rx, zmin, zmax)]
    lvl, wid = CT_WINDOW
    vmin, vmax = lvl - wid / 2.0, lvl + wid / 2.0
    sc = _stride(max(len(ct.x), len(ct.y)), VIZ_MAX_PX) if ct is not None else 1
    xs_w, ys_w = np.round(gx[::s], 3), np.round(gy[::s], 3)
    traces, planes = [], []
    for k in range(len(gz)):
        start = len(traces)
        vis = bool(k == k_init)
        if ct is not None:
            ck = ct.nearest_k(float(gz[k]))
            z_ct = np.round(ct.hu[ck, ::sc, ::sc]).astype(int) if ck is not None else [[None]]
            traces.append(go.Heatmap(
                x=np.round(ct.x[::sc], 3) if ck is not None else [0.0],
                y=np.round(ct.y[::sc], 3) if ck is not None else [0.0],
                z=z_ct, colorscale="Gray", zmin=vmin, zmax=vmax, showscale=False,
                hoverinfo="skip", name="CT", visible=vis))
        wash = art.dose_fine[k, ::s, ::s].astype(np.float64)
        wash = np.where(np.isnan(wash) | (wash < zmin), np.nan, np.round(wash, 3))
        traces.append(go.Heatmap(
            x=xs_w, y=ys_w, z=wash, colorscale=scale, zmin=zmin, zmax=zmax, opacity=0.45,
            zsmooth=False, hoverongaps=False, colorbar=dict(title="Gy", x=1.02, len=0.8),
            name="Dosis", hovertemplate="x %{x:.2f}  y %{y:.2f}<br>%{z:.2f} Gy<extra>Dosis</extra>",
            visible=vis))
        for sp in specs:
            if sp.contours_by_z is None:
                continue
            xs, ys, _ = _polyline_xyz(sp.contours_by_z.get(float(gz[k]), []))
            iso = sp.kind == "isodose"
            traces.append(go.Scatter(
                x=xs, y=ys, mode="lines",
                line=dict(color=_rgb(sp.color), width=2.0 if iso else 1.5, dash="solid" if iso else "dot"),
                name=sp.name, legendgroup=sp.name, visible=vis, hoverinfo="name"))
        for tm in art.targets.values():
            xs, ys, _ = _polyline_xyz(target_rings.get(tm.name, {}).get(k, []))
            traces.append(go.Scatter(
                x=xs, y=ys, mode="lines", line=dict(color=_rgb(tm.color), width=2.5, dash="dash"),
                name=f"Ziel {tm.name} (Original)", legendgroup=f"ziel-{tm.name}", visible=vis,
                hoverinfo="name"))
        planes.append((start, len(traces)))
    fig = go.Figure(data=traces)
    n = len(traces)
    steps = []
    for k, (a, b) in enumerate(planes):
        vis = [False] * n
        for t in range(a, b):
            vis[t] = True
        steps.append(dict(method="restyle", args=[{"visible": vis}], label=f"{gz[k]:.1f}"))
    (xlo, ylo, _), (xhi, yhi, _) = grid.bbox
    if ct is not None:
        xlo, xhi = min(xlo, ct.x[0] - ct.dx / 2), max(xhi, ct.x[-1] + ct.dx / 2)
        ylo, yhi = min(ylo, ct.y[0] - ct.dy / 2), max(yhi, ct.y[-1] + ct.dy / 2)
    fig.update_layout(
        title=dict(text="Axialer Schichtbrowser: CT (W400/L40), Dosis-Wash (>= 10 % Rx), "
                        "exportierte Konturen und Original-Ziel je Ebene", font=dict(size=14)),
        xaxis=dict(title="X [mm] (Links)", range=[float(xlo), float(xhi)], constrain="domain", **_GRID),
        yaxis=dict(title="Y [mm] (Posterior)", range=[float(yhi), float(ylo)], scaleanchor="x",
                   scaleratio=1, **_GRID),
        sliders=[dict(active=int(k_init), currentvalue=dict(prefix="Ebene z = ", suffix=" mm"),
                      pad=dict(t=40), steps=steps)],
        plot_bgcolor="black", legend=dict(x=1.12, y=1.0, bgcolor="rgba(0,0,0,0.4)"),
        height=800, margin=dict(l=50, r=40, t=50, b=40), **_DARK,
    )
    return fig


def _fig_dvh(go, art, curves: dict):
    fig = go.Figure()
    edvh = getattr(art, "eclipse_dvh", {}) or {}
    for tm in art.targets.values():
        cv = curves[tm.name]
        col = _rgb(tm.color)
        fig.add_trace(go.Scatter(
            x=np.round(cv.dose_gy, 4), y=np.round(cv.volume_cm3, 5), mode="lines",
            name=f"{tm.name} (Tool)", line=dict(color=col, width=2.5),
            hovertemplate="%{x:.2f} Gy: %{y:.3f} cm3<extra>" + _esc(tm.name) + "</extra>"))
        if cv.n_samples:
            fig.add_trace(go.Scatter(
                x=[cv.d98_gy, cv.d50_gy, cv.d2_gy],
                y=[0.98 * cv.total_cm3, 0.5 * cv.total_cm3, 0.02 * cv.total_cm3],
                mode="markers+text",
                text=[f"D98 {cv.d98_gy:.2f}", f"D50 {cv.d50_gy:.2f}", f"D2 {cv.d2_gy:.2f}"],
                textposition="top right", textfont=dict(size=10),
                marker=dict(color=col, size=8, symbol="diamond"),
                name=f"{tm.name} D98/D50/D2", showlegend=False, hoverinfo="text"))
        e = edvh.get(tm.roi_number)
        if e is not None:
            fig.add_trace(go.Scatter(
                x=np.round(np.asarray(e.dose_gy, dtype=float), 4),
                y=np.round(np.asarray(e.volume_cm3, dtype=float), 5), mode="lines",
                name=f"{tm.name} (Eclipse-DVH)", line=dict(color=col, width=2, dash="dash"),
                hovertemplate="%{x:.2f} Gy: %{y:.3f} cm3<extra>Eclipse</extra>"))
    fig.add_vline(x=float(art.rx_gy), line_dash="dot", line_color="white",
                  annotation_text=f"Rx {art.rx_gy:.2f} Gy", annotation_position="top left",
                  annotation_font_color="white")
    fig.update_layout(
        title=dict(text="Kumulatives DVH der Ziele (gewichtete Feingitter-Stichproben; "
                        "Eclipse-DVH aus der RTDOSE gestrichelt)", font=dict(size=14)),
        xaxis=dict(title="Dosis [Gy]", rangemode="tozero", **_GRID),
        yaxis=dict(title="Volumen [cm3]", rangemode="tozero", **_GRID),
        plot_bgcolor="rgb(20, 20, 30)", legend=dict(x=0.01, y=0.01, bgcolor="rgba(0,0,0,0.4)"),
        height=500, margin=dict(l=50, r=30, t=50, b=40), **_DARK,
    )
    return fig


# ---------------------------------------------------------------------------
# 5. HTML-Seite (Kopf, Indextabelle, drei Figuren)
# ---------------------------------------------------------------------------

_TABLE_ROWS = (
    ("TV [cm3]", "components", "tv_cm3", 3), ("PIV [cm3]", "components", "piv_cm3", 3),
    ("TV&PIV [cm3]", "components", "tv_piv_cm3", 3), ("PIV50 [cm3]", "components", "piv50_cm3", 3),
    ("CI Paddick", "indices", "ci_paddick", 3), ("Coverage", "indices", "coverage", 3),
    ("Selektivitaet", "indices", "selectivity", 3), ("CI RTOG", "indices", "ci_rtog", 3),
    ("Dice", "indices", "dice", 3), ("GI", "indices", "gi", 2), ("GM [cm]", "indices", "gm_cm", 2),
    ("HI ICRU83", "indices", "hi_icru83", 3), ("D2 [Gy]", "dvh_stats", "d2_gy", 2),
    ("D50 [Gy]", "dvh_stats", "d50_gy", 2), ("D95 [Gy]", "dvh_stats", "d95_gy", 2),
    ("D98 [Gy]", "dvh_stats", "d98_gy", 2), ("Dmin [Gy]", "dvh_stats", "dmin_gy", 2),
    ("Dmax [Gy]", "dvh_stats", "dmax_gy", 2), ("Dmean [Gy]", "dvh_stats", "dmean_gy", 2),
    ("V95 [%]", "dvh_stats", "v95_pct", 1), ("V100 [%]", "dvh_stats", "v100_pct", 1),
)


def _header_html(art, title: str, ct_loaded: bool) -> str:
    s = art.settings or {}
    grid = art.grid
    bits = [f"Rx {art.rx_gy:.2f} Gy ({s.get('rx_source', art.rx_source)})",
            f"Gitter {float(s.get('grid_mm', grid.res_xy)):g} mm (z {grid.dz:.2f} mm)",
            f"Interpolation {s.get('dose_interp', '-')}", f"Volumenmodell {s.get('volume_model', '-')}",
            f"PIV-Scope {s.get('piv_scope', '-')}", f"Isodosen-Konturen {s.get('iso_contours', 'mask')}"]
    if s.get("eclipse_compat"):
        bits.append(f"Eclipse-kompatibel: {s['eclipse_compat']}")
    bits.append("CT-Hintergrund: " + ("ja" if ct_loaded else "nein"))
    levels = ", ".join(f"{lv.label} = {lv.gy:.2f} Gy ({lv.volume_cm3:.3f} cm3, {lv.n_components} Komp.)"
                       for lv in art.levels.values())
    return (f"<h1>{_esc(title)}</h1>\n<div class='meta'>{_esc(' | '.join(bits))}<br>"
            f"Isodosen (global): {_esc(levels)}<br>"
            f"<span class='note'>Feingitter {grid.shape[2]} x {grid.shape[1]} x {grid.shape[0]} Voxel; "
            f"Wash/CT auf max. {VIZ_MAX_PX} px je Achse gestrided, Konturen exakt.</span></div>\n")


def _index_table_html(art) -> str:
    ecl = getattr(art, "eclipse", {}) or {}
    names = list(art.targets)
    has_ecl = any(n in ecl for n in names)
    head = ["<tr><th>Kennwert</th>"]
    for n in names:
        head.append(f"<th>{_esc(n)}<br><span class='src'>Tool</span></th>")
        if has_ecl:
            head.append(f"<th>{_esc(n)}<br><span class='src'>Eclipse</span></th><th>Diff %</th>")
    head.append("</tr>")
    rows = ["<table>", "".join(head)]
    for label, block, key, nd in _TABLE_ROWS:
        cells = [f"<td>{_esc(label)}</td>"]
        for n in names:
            r = art.targets[n].result
            cells.append(f"<td>{_fmt(r[block].get(key), nd)}</td>")
            if has_ecl:
                rw = next((x for x in ecl.get(n, {}).get("rows", []) if x["key"] == key), None)
                if rw is None:
                    cells.append("<td></td><td></td>")
                else:
                    flag = " class='flag'" if rw.get("within_tol") is False else ""
                    src = f" <span class='src'>{_esc(rw['source'])}</span>" if rw.get("source") else ""
                    cells.append(f"<td{flag}>{_fmt(rw['eclipse'], nd)}{src}</td>"
                                 f"<td{flag}>{_fmt(rw['diff_pct'], 1)}</td>")
        rows.append("<tr>" + "".join(cells) + "</tr>")
    for label, kind in (("Schnitt [cm3]", "intersection"), ("Unterdosiert [cm3]", "underdosed"),
                        ("Spill [cm3]", "spill")):
        cells = [f"<td>{_esc(label)}</td>"]
        for n in names:
            tm = art.targets[n]
            v = tm.result["helper_rois"][kind]["volume_cm3"]
            nm = tm.helper_names.get(kind)
            cells.append(f"<td>{_fmt(v, 3)}" + (f" <span class='src'>{_esc(nm)}</span>" if nm else "") + "</td>")
            if has_ecl:
                cells.append("<td></td><td></td>")
        rows.append("<tr>" + "".join(cells) + "</tr>")
    rows.append("</table>")
    if has_ecl:
        notes = []
        for n in names:
            ec = ecl.get(n)
            if not ec:
                continue
            notes.append(f"{n}: Toleranz {ec['tol_pct']:g} %, {ec['n_flagged']} von {ec['n_compared']} "
                         f"Werten ausserhalb (Quellen: {', '.join(ec['sources'])})"
                         + (f"; {ec['piv_scope_note']}" if ec.get("piv_scope_note") else ""))
        rows.append("<div class='note'>Abgleich Eclipse: " + _esc("; ".join(notes)) + "</div>")
    return "\n".join(rows) + "\n"


def _assemble_html(title: str, header: str, table: str, figs: list, out_path: Path) -> None:
    import plotly.io as pio
    from plotly.offline import get_plotlyjs
    parts = [f"<!doctype html>\n<html lang=\"de\">\n<head>\n<meta charset=\"utf-8\">\n"
             f"<title>{_esc(title)}</title>\n<style>{_CSS}</style>\n"
             f"<script>{get_plotlyjs()}</script>\n</head>\n<body>\n", header, "<h2>Indizes</h2>\n", table]
    for div_id, fig, caption in figs:
        parts.append(f"<h2>{_esc(caption)}</h2>\n")
        parts.append(pio.to_html(fig, full_html=False, include_plotlyjs=False, div_id=div_id))
        parts.append("\n")
    parts.append("</body>\n</html>\n")
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("".join(parts))


# ---------------------------------------------------------------------------
# 6. Statisches PNG (axial / koronal / sagittal + DVH)
# ---------------------------------------------------------------------------

def _render_png(art, specs: list, ct: Optional[CtWindow], curves: dict, target_rings: dict,
                kji: tuple, title: str, out_path: Path) -> None:
    # OO-API statt pyplot (Plan P0.7): kein globaler Zustand, kein Backend-Wechsel
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from matplotlib.figure import Figure
    from matplotlib.lines import Line2D

    grid, rx = art.grid, float(art.rx_gy)
    gx, gy, gz = grid.gx, grid.gy, grid.gz
    res, dz = grid.res_xy, grid.dz
    k0, j0, i0 = kji
    zmin, zmax = _dose_range(art)
    cmap = LinearSegmentedColormap.from_list(
        "dose", [(p, _mpl_color(c)) for p, c in _dose_anchors(rx, zmin, zmax)])
    norm = Normalize(vmin=zmin, vmax=zmax)
    lvl, wid = CT_WINDOW
    vmin, vmax = lvl - wid / 2.0, lvl + wid / 2.0

    def ext(a, da, b, db):
        return (float(a[0] - da / 2), float(a[-1] + da / 2), float(b[0] - db / 2), float(b[-1] + db / 2))

    def wash(plane):
        p = np.asarray(plane, dtype=float)
        return np.ma.masked_where(~(p >= zmin), p)     # NaN und < 10 % Rx maskiert

    def iso_lines(ax, GA, GB, plane):
        field = np.nan_to_num(np.asarray(plane, dtype=float), nan=0.0)
        for lv in art.levels.values():
            if (field >= lv.gy).any() and (field < lv.gy).any():
                ax.contour(GA, GB, field, levels=[lv.gy], colors=[_mpl_color(lv.color)], linewidths=1.3)

    def target_lines(ax, GA, GB, plane_mask, color):
        m = np.asarray(plane_mask, dtype=float)
        if m.any() and not m.all():
            ax.contour(GA, GB, m, levels=[0.5], colors=[color], linewidths=1.5, linestyles="dashed")

    (xlo, ylo, zlo), (xhi, yhi, zhi) = grid.bbox
    if ct is not None:
        xlo, xhi = min(xlo, ct.x[0] - ct.dx / 2), max(xhi, ct.x[-1] + ct.dx / 2)
        ylo, yhi = min(ylo, ct.y[0] - ct.dy / 2), max(yhi, ct.y[-1] + ct.dy / 2)
        zlo, zhi = min(zlo, ct.z[0] - ct.dz / 2), max(zhi, ct.z[-1] + ct.dz / 2)

    fig = Figure(figsize=(24, 6.5))
    FigureCanvasAgg(fig)
    axes = fig.subplots(1, 4, gridspec_kw={"width_ratios": [1, 1, 1, 1.15]})

    # --- Axial (x-y bei z = gz[k0]) ---
    ax = axes[0]
    if ct is not None:
        ck = ct.nearest_k(float(gz[k0]))
        if ck is not None:
            ax.imshow(ct.hu[ck], cmap="gray", vmin=vmin, vmax=vmax, extent=ext(ct.x, ct.dx, ct.y, ct.dy),
                      origin="lower", interpolation="nearest", aspect="equal")
    ax.imshow(wash(art.dose_fine[k0]), cmap=cmap, norm=norm, alpha=0.45, extent=ext(gx, res, gy, res),
              origin="lower", interpolation="nearest", aspect="equal")
    for sp in specs:
        if not sp.contours_by_z:
            continue
        iso = sp.kind == "isodose"
        for r in sp.contours_by_z.get(float(gz[k0]), []):
            c = np.vstack([r, r[:1]])
            ax.plot(c[:, 0], c[:, 1], color=_mpl_color(sp.color), lw=1.3 if iso else 1.0, ls="-" if iso else ":")
    for tm in art.targets.values():
        for r in target_rings.get(tm.name, {}).get(k0, []):
            c = np.vstack([r, r[:1]])
            ax.plot(c[:, 0], c[:, 1], color=_mpl_color(tm.color), lw=1.5, ls="--")
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(ylo, yhi)
    ax.invert_yaxis()
    ax.set_title(f"Axial  z = {gz[k0]:.1f} mm")
    ax.set_xlabel("X [mm] (Links)")
    ax.set_ylabel("Y [mm] (Posterior)")

    # --- Koronal (x-z bei y = gy[j0]) ---
    ax = axes[1]
    if ct is not None:
        jc = int(np.argmin(np.abs(ct.y - gy[j0])))
        ax.imshow(ct.hu[:, jc, :], cmap="gray", vmin=vmin, vmax=vmax, extent=ext(ct.x, ct.dx, ct.z, ct.dz),
                  origin="lower", interpolation="nearest", aspect="equal")
    plane = art.dose_fine[:, j0, :]
    ax.imshow(wash(plane), cmap=cmap, norm=norm, alpha=0.45, extent=ext(gx, res, gz, dz),
              origin="lower", interpolation="nearest", aspect="equal")
    GX, GZ = np.meshgrid(gx, gz)
    iso_lines(ax, GX, GZ, plane)
    for tm in art.targets.values():
        target_lines(ax, GX, GZ, tm.structure.mask[:, j0, :], _mpl_color(tm.color))
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(zlo, zhi)
    ax.set_title(f"Koronal  y = {gy[j0]:.1f} mm")
    ax.set_xlabel("X [mm] (Links)")
    ax.set_ylabel("Z [mm] (Superior)")

    # --- Sagittal (y-z bei x = gx[i0]) ---
    ax = axes[2]
    if ct is not None:
        ic = int(np.argmin(np.abs(ct.x - gx[i0])))
        ax.imshow(ct.hu[:, :, ic], cmap="gray", vmin=vmin, vmax=vmax, extent=ext(ct.y, ct.dy, ct.z, ct.dz),
                  origin="lower", interpolation="nearest", aspect="equal")
    plane = art.dose_fine[:, :, i0]
    ax.imshow(wash(plane), cmap=cmap, norm=norm, alpha=0.45, extent=ext(gy, res, gz, dz),
              origin="lower", interpolation="nearest", aspect="equal")
    GY, GZ = np.meshgrid(gy, gz)
    iso_lines(ax, GY, GZ, plane)
    for tm in art.targets.values():
        target_lines(ax, GY, GZ, tm.structure.mask[:, :, i0], _mpl_color(tm.color))
    ax.set_xlim(ylo, yhi)
    ax.set_ylim(zlo, zhi)
    ax.set_title(f"Sagittal  x = {gx[i0]:.1f} mm")
    ax.set_xlabel("Y [mm] (Posterior)")
    ax.set_ylabel("Z [mm] (Superior)")

    # --- DVH ---
    ax = axes[3]
    edvh = getattr(art, "eclipse_dvh", {}) or {}
    for tm in art.targets.values():
        cv = curves[tm.name]
        col = _mpl_color(tm.color)
        ax.plot(cv.dose_gy, cv.volume_cm3, color=col, lw=2, label=f"{tm.name} (Tool)")
        if cv.n_samples:
            ax.plot([cv.d98_gy, cv.d50_gy, cv.d2_gy],
                    [0.98 * cv.total_cm3, 0.5 * cv.total_cm3, 0.02 * cv.total_cm3],
                    "D", color=col, ms=5)
            for dval, frac, lab in ((cv.d98_gy, 0.98, "D98"), (cv.d50_gy, 0.5, "D50"), (cv.d2_gy, 0.02, "D2")):
                ax.annotate(f"{lab} {dval:.2f}", (dval, frac * cv.total_cm3), textcoords="offset points",
                            xytext=(4, 4), fontsize=8, color=col)
        e = edvh.get(tm.roi_number)
        if e is not None:
            ax.plot(np.asarray(e.dose_gy, dtype=float), np.asarray(e.volume_cm3, dtype=float),
                    color=col, lw=1.5, ls="--", label=f"{tm.name} (Eclipse-DVH)")
    ax.axvline(rx, color="k", ls=":", lw=1)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Dosis [Gy]")
    ax.set_ylabel("Volumen [cm3]")
    ax.set_title(f"DVH  (Rx {rx:.2f} Gy)")
    ax.grid(alpha=0.3, ls="--")
    ax.legend(fontsize=8, loc="lower left")

    handles, labels = [], []
    for lv in art.levels.values():
        handles.append(Line2D([0], [0], color=_mpl_color(lv.color), lw=1.5))
        labels.append(f"Isodose {lv.label} ({lv.gy:.2f} Gy)")
    for sp in specs:
        if sp.contours_by_z and sp.kind != "isodose":
            handles.append(Line2D([0], [0], color=_mpl_color(sp.color), lw=1.0, ls=":"))
            labels.append(sp.name)
    for tm in art.targets.values():
        handles.append(Line2D([0], [0], color=_mpl_color(tm.color), lw=1.5, ls="--"))
        labels.append(f"Ziel {tm.name}")
    fig.legend(handles, labels, loc="lower center", ncol=min(6, max(1, len(handles))), fontsize=8, frameon=False)
    fig.suptitle(title, fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0, 0.07, 1, 0.95))
    fig.savefig(out_path, dpi=150, bbox_inches="tight")


# ---------------------------------------------------------------------------
# 7. Einstieg
# ---------------------------------------------------------------------------

def run_dose_visualization(art, specs: list, ct_index: Optional[dict], out_dir, *,
                           ct_background: bool = True, case_id: str = "", label: str = "",
                           verbose: bool = True) -> dict:
    """
    Schreibt ``validation.html`` und ``dose_overview.png`` nach ``out_dir`` und
    liefert ``{'viz_html_path', 'viz_png_path'}`` (None bei Fehlschlag).  HTML
    und PNG sind getrennt abgesichert; ein Fehler in einem verhindert das
    andere nicht.  ``ct_index`` (``rtstruct_writer.build_ct_slice_index``) ist
    optional; ohne ``z_to_path`` oder mit ``ct_background=False`` gibt es
    keinen CT-Hintergrund.
    """
    def say(msg):
        if verbose:
            print(msg)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out = {"viz_html_path": None, "viz_png_path": None}
    if not art.targets:
        say("  Visualisierung uebersprungen: keine Ziele.")
        return out
    grid = art.grid
    notes = []
    ct = None
    if ct_background and ct_index:
        try:
            ct = load_ct_window(ct_index, grid, notes=notes)
        except Exception as e:  # noqa: BLE001 - CT ist nur Hintergrund
            notes.append(f"CT-Hintergrund nicht geladen ({type(e).__name__}: {e}).")
    for n in notes:
        say(f"  Hinweis: {n}")
    first = next(iter(art.targets.values()))
    _centroid, kji = _mask_centroid(first.structure.mask, grid)
    target_rings = {tm.name: _rings_by_plane(tm.contours, grid) for tm in art.targets.values()}
    curves = {tm.name: build_dvh_curve(tm.dose_samples, tm.sample_weights, grid.voxel_volume_mm3)
              for tm in art.targets.values()}
    title = f"Dosisindex-Validierung {case_id}{label}".strip()

    try:
        import plotly.graph_objects as go
        figs = [
            ("viz-3d", _fig_3d(go, art, specs), "3D: Isodosen-/Zielflaechen und exportierte Konturen"),
            ("viz-slices", _fig_slices(go, art, specs, ct, kji[0], target_rings),
             "Axialer Schichtbrowser (CT, Dosis-Wash, Konturen)"),
            ("viz-dvh", _fig_dvh(go, art, curves), "Kumulatives DVH"),
        ]
        html_path = out_dir / HTML_NAME
        _assemble_html(title, _header_html(art, title, ct is not None), _index_table_html(art), figs, html_path)
        out["viz_html_path"] = str(html_path)
        say(f"  Gespeichert: {html_path} ({html_path.stat().st_size / 1e6:.1f} MB)")
    except ImportError:
        say("  Hinweis: plotly nicht installiert -- validation.html uebersprungen.")
    except Exception as e:  # noqa: BLE001
        say(f"  (!) validation.html uebersprungen: {type(e).__name__}: {e}")

    try:
        png_path = out_dir / PNG_NAME
        _render_png(art, specs, ct, curves, target_rings, kji, title, png_path)
        out["viz_png_path"] = str(png_path)
        say(f"  Gespeichert: {png_path}")
    except Exception as e:  # noqa: BLE001
        say(f"  (!) dose_overview.png uebersprungen: {type(e).__name__}: {e}")
    return out
