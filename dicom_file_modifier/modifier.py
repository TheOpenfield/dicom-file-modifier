#!/usr/bin/env python3
"""
CT DICOM Rigid Body Transformer
================================
Wendet eine starre Körper-Transformation (Translation + Rotation) auf eine
CT-DICOM-Serie an und speichert das Ergebnis als neue DICOM-Slices.

Zwei Methoden:
  resample  – Neuabtastung des Voxelgitters via inverses Mapping + trilineare
              Interpolation.  Minimale HU-Wert-Abweichung (~1–2 HU für Gewebe).
              Erzeugt achsiale Standard-Slices, kompatibel mit allen TPS.
  metadata  – Nur DICOM-Tags (ImagePositionPatient, ImageOrientationPatient)
              werden angepasst, Pixeldaten bleiben byte-identisch (HU exakt).
              Einschränkung: Manche Planungssysteme lehnen schräge IOP ab.

Rotationskonvention:
  Extrinsisch XYZ – erst rx um feste X-Achse, dann ry um feste Y-Achse,
  dann rz um feste Z-Achse (DICOM-Patientenkoordinaten: X=links, Y=posterior,
  Z=superior).  Rotationszentrum = geometrischer Mittelpunkt des Volumens.

Verwendung:
  python -m dicom_file_modifier.modifier <CT-Verzeichnis> [Optionen]

Beispiele:
  # 10 mm nach rechts verschieben + 15° um Z drehen, Ergebnis anzeigen
  python -m dicom_file_modifier.modifier data/0000000171/CT --tx 10 --rz 15

  # Nur Metadaten, Output in eigenem Verzeichnis
  python -m dicom_file_modifier.modifier data/0000000171/CT --ty -5 --rx 3 \\
      --method metadata --output output/ct_meta
"""

import os
import copy
import glob
import argparse
import numpy as np
import pydicom
from pydicom.uid import generate_uid
from scipy.ndimage import map_coordinates, spline_filter
from scipy.spatial.transform import Rotation


# ─────────────────────────────────────────────────────────────────────────────
#  Laden
# ─────────────────────────────────────────────────────────────────────────────

def load_ct_series(ct_dir: str) -> list:
    """Lädt alle DICOM-Dateien aus einem Verzeichnis, sortiert nach Z-Position."""
    files = sorted(glob.glob(os.path.join(ct_dir, "*.dcm")))
    if not files:
        raise FileNotFoundError(f"Keine DICOM-Dateien in {ct_dir!r}")

    slices = []
    for f in files:
        try:
            ds = pydicom.dcmread(f)
            if hasattr(ds, "ImagePositionPatient") and hasattr(ds, "PixelData"):
                slices.append(ds)
        except Exception:
            pass

    if not slices:
        raise ValueError(f"Keine gültigen CT-Slices in {ct_dir!r}")

    # Nach Z-Koordinate der Slice-Normalenkomponente sortieren
    slices.sort(key=lambda s: float(s.ImagePositionPatient[2]))
    return slices


def slices_to_hu(slices: list) -> np.ndarray:
    """Stapelt Pixel-Arrays und wandelt sie in Hounsfield-Einheiten (HU) um."""
    arrays = []
    for s in slices:
        arr = s.pixel_array.astype(np.float32)
        slope = float(getattr(s, "RescaleSlope", 1.0))
        intercept = float(getattr(s, "RescaleIntercept", 0.0))
        arrays.append(arr * slope + intercept)
    return np.stack(arrays, axis=0)  # (nz, ny, nx)


def extract_geometry(slices: list) -> dict:
    """
    Extrahiert räumliche Geometrie aus DICOM-Metadaten.

    Affine-Konvention:  P_patient = A @ [k, j, i, 1]^T
      k = Slice-Index, j = Zeilen-Index, i = Spalten-Index
    """
    iop = np.array([float(x) for x in slices[0].ImageOrientationPatient])
    # DICOM-Definition:
    #   iop[:3] = Richtungsvektor entlang steigendem Spalten-Index (= "row direction")
    #   iop[3:] = Richtungsvektor entlang steigendem Zeilen-Index  (= "col direction")
    row_dir = iop[:3]
    col_dir = iop[3:]
    normal  = np.cross(row_dir, col_dir)          # Slice-Normalenvektor

    ps = slices[0].PixelSpacing
    dr = float(ps[0])   # Zeilen-Abstand    [mm]  (PixelSpacing[0])
    dc = float(ps[1])   # Spalten-Abstand   [mm]  (PixelSpacing[1])

    ipp_list = [
        np.array([float(x) for x in s.ImagePositionPatient]) for s in slices
    ]

    if len(slices) > 1:
        dz = float(np.linalg.norm(ipp_list[1] - ipp_list[0]))
    else:
        dz = float(getattr(slices[0], "SliceThickness", 1.0))

    # Affine 4×4: Voxel (k, j, i) → Patient (x, y, z)
    # P = IPP[0]  +  k·dz·normal  +  j·dr·col_dir  +  i·dc·row_dir
    A = np.eye(4)
    A[:3, 0] = normal  * dz   # k-Achse (Slice-Richtung)
    A[:3, 1] = col_dir * dr   # j-Achse (Zeilen)
    A[:3, 2] = row_dir * dc   # i-Achse (Spalten)
    A[:3, 3] = ipp_list[0]    # Ursprung = IPP des ersten Slices

    return {
        "affine":   A,
        "row_dir":  row_dir,
        "col_dir":  col_dir,
        "normal":   normal,
        "dr": dr, "dc": dc, "dz": dz,
        "ipp":      ipp_list,
        "shape":    (len(slices), slices[0].Rows, slices[0].Columns),
    }


# ─────────────────────────────────────────────────────────────────────────────
#  Transformation
# ─────────────────────────────────────────────────────────────────────────────

def volume_center(geom: dict) -> np.ndarray:
    """Geometrischer Mittelpunkt des CT-Volumens in Patienten-Koordinaten [mm]."""
    nz, ny, nx = geom["shape"]
    c_vox = np.array([(nz - 1) / 2.0, (ny - 1) / 2.0, (nx - 1) / 2.0, 1.0])
    return (geom["affine"] @ c_vox)[:3]


def build_rigid_transform(
    rx: float, ry: float, rz: float,
    tx: float, ty: float, tz: float,
    center: np.ndarray,
) -> np.ndarray:
    """
    Erstellt eine 4×4 starre Transformationsmatrix.

    Reihenfolge:
      1. Rotation extrinsisch XYZ um das Volumenzentrum
      2. Translation

    Vorwärts-Transform (Original → Neu):
      P_neu = R · (P - C) + C + t
            = R · P  +  (−R·C + C + t)
    """
    R = Rotation.from_euler("XYZ", [rx, ry, rz], degrees=True).as_matrix()
    T = np.eye(4)
    T[:3, :3] = R
    T[:3, 3]  = -R @ center + center + np.array([tx, ty, tz])
    return T


def resample_volume(
    volume_hu: np.ndarray,
    affine: np.ndarray,
    T: np.ndarray,
    order: int = 1,
    chunk_slices: int = 20,
) -> np.ndarray:
    """
    Neuabtastung des CT-Volumens nach einer Starrkörper-Transformation.

    Inverses Mapping:  Für jeden Ausgabe-Voxel P_out wird die Quellposition
      P_in = T⁻¹(P_out)  berechnet und der HU-Wert interpoliert.

    Das Ausgabe-Volumen hat dieselbe Geometrie (IPP/IOP/Größe) wie die Eingabe.
    Bereiche außerhalb des Original-Volumens werden mit −1000 HU (Luft) gefüllt.

    Für order > 1 wird der Spline-Filter vorab auf das Gesamtvolumen angewandt,
    danach wird in Scheiben (chunk_slices) verarbeitet, um Speicher zu sparen.
    """
    T_inv   = np.linalg.inv(T)
    A_inv   = np.linalg.inv(affine)
    M       = A_inv @ T_inv @ affine   # kombinierte Voxel→Voxel-Abbildung

    nz, ny, nx = volume_hu.shape

    # Vorfilterung für kubische Interpolation (einmalig für das Gesamtvolumen)
    if order > 1:
        vol_proc = spline_filter(volume_hu.astype(np.float64), order=order,
                                 output=np.float64)
    else:
        vol_proc = volume_hu.astype(np.float64)

    output = np.full((nz, ny, nx), -1000.0, dtype=np.float32)

    for k0 in range(0, nz, chunk_slices):
        k1 = min(k0 + chunk_slices, nz)
        nk = k1 - k0
        n  = nk * ny * nx

        k_g, j_g, i_g = np.mgrid[k0:k1, 0:ny, 0:nx]

        coords_out    = np.empty((4, n), dtype=np.float64)
        coords_out[0] = k_g.ravel()
        coords_out[1] = j_g.ravel()
        coords_out[2] = i_g.ravel()
        coords_out[3] = 1.0

        coords_in = M @ coords_out   # (4, n)

        vals = map_coordinates(
            vol_proc,
            [coords_in[0], coords_in[1], coords_in[2]],
            order=order,
            mode="constant",
            cval=-1000.0,
            prefilter=False,   # Vorfilterung bereits erledigt
        )
        output[k0:k1] = vals.reshape(nk, ny, nx)

    return output


def apply_metadata_transform(slices: list, T: np.ndarray) -> list:
    """
    Aktualisiert ImagePositionPatient und ImageOrientationPatient (Vorwärts-Transform).
    Pixeldaten bleiben byte-identisch → HU-Werte exakt erhalten.
    """
    R = T[:3, :3]
    new_slices = []
    for ds in slices:
        nd = copy.deepcopy(ds)

        # IPP: neue Position = T · alte Position
        old_ipp = np.array([float(x) for x in ds.ImagePositionPatient])
        new_ipp = (T @ np.append(old_ipp, 1.0))[:3]
        nd.ImagePositionPatient = [f"{v:.6f}" for v in new_ipp]

        # IOP: Richtungsvektoren rotieren
        old_iop = np.array([float(x) for x in ds.ImageOrientationPatient])
        new_iop = np.concatenate([R @ old_iop[:3], R @ old_iop[3:]])
        nd.ImageOrientationPatient = [f"{v:.6f}" for v in new_iop]

        nd.SOPInstanceUID = generate_uid()
        new_slices.append(nd)
    return new_slices


# ─────────────────────────────────────────────────────────────────────────────
#  Speichern
# ─────────────────────────────────────────────────────────────────────────────

def save_ct_series(
    slices: list,
    output_dir: str,
    new_volume_hu: "np.ndarray | None" = None,
) -> None:
    """
    Speichert CT-Slices als DICOM-Dateien mit neuer SeriesInstanceUID.

    - new_volume_hu = None  →  Slices enthalten bereits aktualisierte Metadaten
                               (metadata-Methode).  Pixeldaten unverändert.
    - new_volume_hu ≠ None  →  Pixeldaten werden aus dem Volumen geschrieben;
                               Metadaten kommen aus den Original-Slices
                               (resample-Methode).
    """
    os.makedirs(output_dir, exist_ok=True)
    series_uid = generate_uid()

    for k, ds in enumerate(slices):
        nd = copy.deepcopy(ds)
        nd.SeriesInstanceUID = series_uid
        nd.SOPInstanceUID    = generate_uid()

        orig_desc = str(getattr(ds, "SeriesDescription", "CT"))
        nd.SeriesDescription = orig_desc + "_transformed"

        if new_volume_hu is not None:
            slope     = float(getattr(ds, "RescaleSlope",     1.0))
            intercept = float(getattr(ds, "RescaleIntercept", 0.0))

            stored = np.round((new_volume_hu[k].astype(np.float64) - intercept) / slope)

            # Auf gültigen int16-Wertebereich begrenzen
            stored = np.clip(stored, -32768, 32767).astype(np.int16)

            nd.PixelData         = stored.tobytes()
            nd.Rows, nd.Columns  = stored.shape
            nd.BitsAllocated     = 16
            nd.BitsStored        = 16
            nd.HighBit           = 15
            nd.PixelRepresentation = 1   # vorzeichenbehaftet (int16)

        out_path = os.path.join(output_dir, f"CT_{k:04d}.dcm")
        nd.save_as(out_path)

    print(f"  {len(slices)} Slices gespeichert -> {output_dir}")


# ─────────────────────────────────────────────────────────────────────────────
#  3D-Visualisierung
# ─────────────────────────────────────────────────────────────────────────────

def _extract_surface(
    volume_hu: np.ndarray,
    affine: np.ndarray,
    threshold: float = -300.0,
    downsample: int = 2,
):
    """
    Extrahiert eine Isofläche via Marching-Cubes-Algorithmus.

    Gibt (vertices_patient, faces) zurück; vertices in Patienten-Koordinaten [mm].
    Downsample reduziert die Rechenlast (Standard-Faktor 2).
    """
    try:
        from skimage.measure import marching_cubes
    except ImportError:
        raise ImportError(
            "scikit-image wird für die 3D-Visualisierung benötigt: "
            "pip install scikit-image"
        )

    vol_ds = volume_hu[::downsample, ::downsample, ::downsample]

    try:
        verts_vox, faces, _, _ = marching_cubes(vol_ds, level=threshold, step_size=1)
    except (ValueError, RuntimeError):
        return None, None

    # Voxel-Koordinaten auf Original-Auflösung skalieren
    verts_vox_full = verts_vox * downsample

    # Voxel (k, j, i) → Patienten (x, y, z)
    n    = verts_vox_full.shape[0]
    hom  = np.hstack([verts_vox_full, np.ones((n, 1))])   # (N, 4) [k, j, i, 1]
    verts_pat = (affine @ hom.T).T[:, :3]                  # (N, 3) [x, y, z]

    return verts_pat, faces


def visualize_3d(
    volume_original: np.ndarray,
    volume_transformed: np.ndarray,
    geom: dict,
    T: np.ndarray,
    method: str = "resample",
    output_html: "str | None" = None,
) -> None:
    """
    Interaktiver 3D-Plot: Original (blau, halbtransparent) vs. Transformiert (rot).
    Zeigt Körperoberfläche (HU > −300) und optional Knochen (HU > 400).

    Für die metadata-Methode: Vertices werden mit der Vorwärts-Transformation T
    in die neue Position verschoben (Pixeldaten identisch, daher kein Resampling
    sichtbar).
    """
    try:
        import plotly.graph_objects as go
    except ImportError:
        raise ImportError("plotly wird für die Visualisierung benötigt: pip install plotly")

    affine = geom["affine"]
    center = volume_center(geom)

    print("  Extrahiere Körperoberfläche (Original) …")
    verts_orig, faces_orig = _extract_surface(volume_original, affine, threshold=-300.0)

    if method == "metadata":
        # Pixeldaten identisch → transformierte Oberfläche durch Vorwärts-Transform
        if verts_orig is not None:
            n        = verts_orig.shape[0]
            hom      = np.hstack([verts_orig, np.ones((n, 1))])
            verts_trans = (T @ hom.T).T[:, :3]
            faces_trans = faces_orig
        else:
            verts_trans = faces_trans = None
    else:
        print("  Extrahiere Körperoberfläche (Transformiert) …")
        verts_trans, faces_trans = _extract_surface(
            volume_transformed, affine, threshold=-300.0
        )

    # ── Knochen ──────────────────────────────────────────────────────────────
    print("  Extrahiere Knochenoberfläche …")
    verts_bone_orig, faces_bone = _extract_surface(
        volume_original, affine, threshold=400.0, downsample=3
    )

    fig = go.Figure()

    # Original-Oberfläche (blau)
    if verts_orig is not None:
        fig.add_trace(go.Mesh3d(
            x=verts_orig[:, 0], y=verts_orig[:, 1], z=verts_orig[:, 2],
            i=faces_orig[:, 0], j=faces_orig[:, 1], k=faces_orig[:, 2],
            color="royalblue", opacity=0.20, name="Original (Körper)",
            showlegend=True,
            lighting=dict(diffuse=0.8, specular=0.3, roughness=0.5),
            lightposition=dict(x=1, y=1, z=2),
        ))

    # Transformierte Oberfläche (rot)
    if verts_trans is not None:
        fig.add_trace(go.Mesh3d(
            x=verts_trans[:, 0], y=verts_trans[:, 1], z=verts_trans[:, 2],
            i=faces_trans[:, 0], j=faces_trans[:, 1], k=faces_trans[:, 2],
            color="tomato", opacity=0.35, name="Transformiert (Körper)",
            showlegend=True,
            lighting=dict(diffuse=0.8, specular=0.3, roughness=0.5),
            lightposition=dict(x=1, y=1, z=2),
        ))

    # Knochen (gelb)
    if verts_bone_orig is not None:
        fig.add_trace(go.Mesh3d(
            x=verts_bone_orig[:, 0], y=verts_bone_orig[:, 1], z=verts_bone_orig[:, 2],
            i=faces_bone[:, 0], j=faces_bone[:, 1], k=faces_bone[:, 2],
            color="gold", opacity=0.70, name="Original (Knochen)",
            showlegend=True,
            lighting=dict(diffuse=0.7, specular=0.5),
        ))

    # Rotationszentrum
    fig.add_trace(go.Scatter3d(
        x=[center[0]], y=[center[1]], z=[center[2]],
        mode="markers",
        marker=dict(size=8, color="white", symbol="x", line=dict(color="black", width=2)),
        name="Rotationszentrum",
    ))

    # Koordinatenachsen
    axis_len = 60.0
    for axis_vec, color, label in [
        (np.array([1, 0, 0]), "red",   "X (L)"),
        (np.array([0, 1, 0]), "lime",  "Y (P)"),
        (np.array([0, 0, 1]), "cyan",  "Z (S)"),
    ]:
        end = center + axis_vec * axis_len
        fig.add_trace(go.Scatter3d(
            x=[center[0], end[0]],
            y=[center[1], end[1]],
            z=[center[2], end[2]],
            mode="lines+text",
            line=dict(color=color, width=5),
            text=["", label],
            textposition="top center",
            textfont=dict(color=color, size=12),
            name=f"Achse {label}",
            showlegend=False,
        ))

    # Transformationsparameter als Untertitel
    R = T[:3, :3]
    t = T[:3, 3]
    rx_eff, ry_eff, rz_eff = Rotation.from_matrix(R).as_euler("XYZ", degrees=True)

    title_text = (
        "CT Rigid Body Transformation<br>"
        f"<sub>Translation: ({t[0]:.1f}, {t[1]:.1f}, {t[2]:.1f}) mm  |  "
        f"Rotation: Rx={rx_eff:.1f}°  Ry={ry_eff:.1f}°  Rz={rz_eff:.1f}°  |  "
        f"Methode: {method}</sub>"
    )

    fig.update_layout(
        title=dict(text=title_text, font=dict(size=14)),
        scene=dict(
            xaxis_title="X [mm]  (Links)",
            yaxis_title="Y [mm]  (Posterior)",
            zaxis_title="Z [mm]  (Superior)",
            aspectmode="data",
            bgcolor="rgb(20, 20, 30)",
            xaxis=dict(gridcolor="rgba(255,255,255,0.15)", color="white"),
            yaxis=dict(gridcolor="rgba(255,255,255,0.15)", color="white"),
            zaxis=dict(gridcolor="rgba(255,255,255,0.15)", color="white"),
        ),
        legend=dict(x=0.01, y=0.99, bgcolor="rgba(0,0,0,0.4)",
                    font=dict(color="white")),
        paper_bgcolor="rgb(20, 20, 30)",
        font=dict(color="white"),
        margin=dict(l=0, r=0, b=0, t=80),
    )

    if output_html:
        fig.write_html(output_html)
        print(f"  Visualisierung gespeichert -> {output_html}")
    else:
        fig.show()


# ─────────────────────────────────────────────────────────────────────────────
#  Hauptprogramm
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "CT DICOM Rigid Body Transformer – "
            "verschiebt und dreht CT-Slices ohne Bildverzerrung."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Beispiel:\n"
            "  python -m dicom_file_modifier.modifier data/0000000171/CT "
            "--tx 10 --rz 15\n"
        ),
    )
    parser.add_argument("ct_dir",
                        help="Verzeichnis mit CT-DICOM-Dateien")
    parser.add_argument("--output", "-o", default="output/ct_transformed",
                        metavar="DIR",
                        help="Ausgabeverzeichnis  (Standard: output/ct_transformed)")

    grp_t = parser.add_argument_group("Translation [mm]")
    grp_t.add_argument("--tx", type=float, default=0.0, metavar="mm",
                       help="Verschiebung in X-Richtung (Links+)")
    grp_t.add_argument("--ty", type=float, default=0.0, metavar="mm",
                       help="Verschiebung in Y-Richtung (Posterior+)")
    grp_t.add_argument("--tz", type=float, default=0.0, metavar="mm",
                       help="Verschiebung in Z-Richtung (Superior+)")

    grp_r = parser.add_argument_group("Rotation [°]  –  extrinsisch XYZ um Volumenmitte")
    grp_r.add_argument("--rx", type=float, default=0.0, metavar="deg",
                       help="Rotation um X-Achse (Pitch)")
    grp_r.add_argument("--ry", type=float, default=0.0, metavar="deg",
                       help="Rotation um Y-Achse (Roll)")
    grp_r.add_argument("--rz", type=float, default=0.0, metavar="deg",
                       help="Rotation um Z-Achse (Yaw)")

    grp_m = parser.add_argument_group("Methode und Qualität")
    grp_m.add_argument("--method", choices=["resample", "metadata"],
                       default="resample",
                       help=(
                           "resample (Standard): Neuabtastung, Standard-Axial-Slices, "
                           "minimale HU-Änderung. "
                           "metadata: Pixeldaten unverändert, HU exakt erhalten, "
                           "schräge IOP."
                       ))
    grp_m.add_argument("--order", type=int, choices=[0, 1, 3], default=1,
                       metavar="N",
                       help=(
                           "Interpolationsordnung (nur resample): "
                           "0=Nächster Nachbar (keine HU-Änderung, Treppeneffekte), "
                           "1=linear (Standard, ~1–2 HU Abweichung), "
                           "3=kubisch (beste Qualität, ~0.5 HU Abweichung)"
                       ))

    grp_v = parser.add_argument_group("Visualisierung")
    grp_v.add_argument("--save-viz", metavar="FILE.html",
                       help="Visualisierung als HTML-Datei speichern "
                            "(Standard: <output>/visualization_3d.html)")
    grp_v.add_argument("--no-viz", action="store_true",
                       help="Visualisierung überspringen")

    args = parser.parse_args()

    # ── Laden ──────────────────────────────────────────────────────────────
    print(f"\nLade CT-Serie aus {args.ct_dir!r} …")
    slices   = load_ct_series(args.ct_dir)
    print(f"  {len(slices)} Slices geladen")

    volume_hu = slices_to_hu(slices)
    geom      = extract_geometry(slices)
    center    = volume_center(geom)

    nz, ny, nx = geom["shape"]
    print(f"  Volumengröße : {nz} × {ny} × {nx}  Voxel")
    print(f"  Voxelabstand : dz={geom['dz']:.2f} mm  |  "
          f"dr={geom['dr']:.2f} mm  |  dc={geom['dc']:.2f} mm")
    print(f"  Volumen-Mitte: ({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f}) mm")

    # ── Transformationsmatrix ───────────────────────────────────────────────
    print(f"\nTransformation:")
    print(f"  Translation  : tx={args.tx} mm,  ty={args.ty} mm,  tz={args.tz} mm")
    print(f"  Rotation     : rx={args.rx}°,  ry={args.ry}°,  rz={args.rz}°  [extrinsisch XYZ]")
    print(f"  Methode      : {args.method}")

    T = build_rigid_transform(
        args.rx, args.ry, args.rz,
        args.tx, args.ty, args.tz,
        center,
    )
    print(f"  T =\n{np.array2string(T, precision=4, suppress_small=True)}")

    # ── Transformation anwenden ─────────────────────────────────────────────
    if args.method == "metadata":
        print("\nAktualisiere DICOM-Metadaten (HU-Werte exakt erhalten) …")
        transformed_slices = apply_metadata_transform(slices, T)
        save_ct_series(transformed_slices, args.output)
        volume_transformed = volume_hu   # identische Pixeldaten für Viz

    else:
        print(f"\nNeuabtastung (Interpolationsordnung {args.order}) …")
        volume_transformed = resample_volume(
            volume_hu, geom["affine"], T, order=args.order
        )

        # Sanity-Check: HU-Wertebereich muss im Original und Ergebnis gleich sein
        # (Vergleich voxelweise ist sinnlos, da der Körper seine Position gewechselt hat)
        print(f"  HU-Bereich Original    : [{volume_hu.min():.0f}, {volume_hu.max():.0f}] HU")
        print(f"  HU-Bereich Transformiert: [{volume_transformed.min():.0f}, {volume_transformed.max():.0f}] HU")

        save_ct_series(slices, args.output, new_volume_hu=volume_transformed)

    # ── Visualisierung ──────────────────────────────────────────────────────
    if not args.no_viz:
        print("\nErstelle 3D-Visualisierung …")
        html_out = args.save_viz or os.path.join(args.output, "visualization_3d.html")
        visualize_3d(
            volume_original=volume_hu,
            volume_transformed=volume_transformed,
            geom=geom,
            T=T,
            method=args.method,
            output_html=html_out,
        )

    print("\nFertig!")


if __name__ == "__main__":
    main()
