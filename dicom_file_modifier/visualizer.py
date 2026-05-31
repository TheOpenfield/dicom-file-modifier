#!/usr/bin/env python3
"""
RTSTRUCT Visualizer
===================
Erzeugt aussagekraeftige Plots aus RTSTRUCT-Analyseergebnissen.

Plots (Analyse eines einzelnen RTSTRUCT):
  1. Volumen-Balkendiagramm  – alle Strukturen, farbcodiert nach Typ (Target/OAR)
  2. Formmetriken-Vergleich  – Sphaerizitaet, Kompaktheit, Elongation nebeneinander
  3. Abstands-Uebersicht     – Min/Hausdorff/Zentroid pro Strukturpaar als Balken
  4. Schwerpunkt-3D-Karte    – raeumliche Positionen aller Strukturzentren im Patientenraum
  5. statistics.txt          – Zahlenzusammenfassung aller Metriken

Case-Transform-Plots (Vorher/Nachher einer rigiden Transformation, von
``case_modifier`` aufgerufen via :func:`run_case_visualization`):
  A. transform_3d.html       – interaktiver 3D-Vergleich Original vs. Transformiert
                               (Konturen, Rotationszentrum, POINT-Marker und
                               optionale CT-Oberflaeche einzeln zuschaltbar)
  B. transform_overview.png  – statische tri-planare Projektion (axial/koronal/sagittal)
  C. displacement.png        – Schwerpunkt-Verschiebung pro ROI als Balkendiagramm

Verwendung:
  python -m dicom_file_modifier.visualizer <rtstruct.dcm> [Optionen]

Beispiele:
  python -m dicom_file_modifier.visualizer data/0000000171/test/1.dcm --output output/
  python -m dicom_file_modifier.visualizer data/0000000171/test/1.dcm \\
      --targets PTV,CTV --oars Parotis,Blase --output output/
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")          # kein interaktives Fenster notwenig
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

from dicom_file_modifier.analyzer import run_analysis, get_structure_names


# Farben
COLOR_TARGET = "#2471A3"   # Blau fuer Zielgebiete
COLOR_OAR    = "#CB4335"   # Rot fuer Risikoorgane
COLOR_MIN    = "#27AE60"   # Gruen fuer Minimalabstand
COLOR_HAUS   = "#E67E22"   # Orange fuer Hausdorff
COLOR_CENT   = "#8E44AD"   # Lila fuer Schwerpunktabstand

# Case-Transform-Vergleich (Vorher/Nachher)
COLOR_ORIG   = "#5D6D7E"   # Grau-Blau fuer Original-Geometrie
COLOR_TRANS  = "#C0392B"   # Rot fuer transformierte Geometrie
COLOR_CENTER = "#27AE60"   # Gruen fuer das Rotationszentrum (Drehpunkt)
COLOR_MARKER = "#F39C12"   # Orange fuer POINT-Marker


# ---------------------------------------------------------------------------
# Plot 1: Volumen-Balkendiagramm
# ---------------------------------------------------------------------------

def plot_volumes(results: dict, output_dir: Path) -> None:
    """
    Horizontales Balkendiagramm aller Strukturvolumen.
    Targets blau, OARs rot. Sortiert nach Volumen (absteigend).
    Sinnvoller als ein Histogramm, wenn < 20 Strukturen vorliegen.
    """
    names, volumes, colors = [], [], []

    for name, s in results["targets"].items():
        names.append(name)
        volumes.append(s["volume_cm3"])
        colors.append(COLOR_TARGET)

    for name, s in results["oars"].items():
        names.append(name)
        volumes.append(s["volume_cm3"])
        colors.append(COLOR_OAR)

    if not names:
        return

    # Absteigend nach Volumen sortieren
    order = np.argsort(volumes)[::-1]
    names   = [names[i]   for i in order]
    volumes = [volumes[i] for i in order]
    colors  = [colors[i]  for i in order]

    fig, ax = plt.subplots(figsize=(10, max(4, len(names) * 0.45 + 1)))
    bars = ax.barh(names, volumes, color=colors, edgecolor="white", linewidth=0.5)

    # Werte ans Ende der Balken schreiben
    for bar, vol in zip(bars, volumes):
        ax.text(bar.get_width() + 0.01 * max(volumes),
                bar.get_y() + bar.get_height() / 2,
                f"{vol:.1f}", va="center", ha="left", fontsize=8)

    ax.set_xlabel("Volumen (cm³)", fontsize=11)
    ax.set_title("Strukturvolumen (Targets vs. OARs)", fontsize=13, fontweight="bold")
    ax.set_xlim(0, max(volumes) * 1.15)
    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.3, linestyle="--")

    legend = [
        mpatches.Patch(color=COLOR_TARGET, label="Zielgebiet (Target)"),
        mpatches.Patch(color=COLOR_OAR,    label="Risikoorgan (OAR)"),
    ]
    ax.legend(handles=legend, loc="lower right", fontsize=9)

    plt.tight_layout()
    path = output_dir / "volumes.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Gespeichert: {path}")


# ---------------------------------------------------------------------------
# Plot 2: Formmetriken-Vergleich
# ---------------------------------------------------------------------------

def plot_shape_metrics(results: dict, output_dir: Path) -> None:
    """
    Gruppiertes Balkendiagramm der drei Formmetriken pro Struktur.
    Ermoeglicht schnellen visuellen Vergleich zwischen Targets und OARs.
    """
    all_structs = {**results["targets"], **results["oars"]}
    if not all_structs:
        return

    names        = list(all_structs.keys())
    sphericities = [s["shape"]["sphericity"]  for s in all_structs.values()]
    compactness  = [s["shape"]["compactness"] for s in all_structs.values()]
    elongations  = [s["shape"]["elongation"]  for s in all_structs.values()]

    x         = np.arange(len(names))
    n_targets = len(results["targets"])

    fig, axes = plt.subplots(1, 3, figsize=(16, max(4, len(names) * 0.4 + 2)))

    metric_data = [
        (axes[0], sphericities, "Sphaerizitaet",
         "1.0 = perfekte Kugel\n< 1 = unregelmaessig/elongiert"),
        (axes[1], compactness,  "Kompaktheit",
         "1.0 = konvex ausgefuellt\n< 1 = konkav/lueckenhaft"),
        (axes[2], elongations,  "Elongation",
         "> 1 = gestreckt\n= 1 = isotropisch"),
    ]

    for ax, values, title, subtitle in metric_data:
        bar_colors = [COLOR_TARGET if i < n_targets else COLOR_OAR
                      for i in range(len(names))]
        bars = ax.bar(x, values, color=bar_colors, edgecolor="white",
                      linewidth=0.5, width=0.6)

        for bar, val in zip(bars, values):
            if val > 0:
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 0.005 * max(values + [1]),
                        f"{val:.3f}", ha="center", va="bottom", fontsize=7)

        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=40, ha="right", fontsize=8)
        ax.set_title(f"{title}\n{subtitle}", fontsize=10, fontweight="bold")
        ax.set_ylim(0, max(values + [1.05]) * 1.12)
        ax.grid(axis="y", alpha=0.3, linestyle="--")

        # Referenzlinie bei 1.0 (Sphaerizitaet, Kompaktheit)
        if title != "Elongation":
            ax.axhline(1.0, color="gray", linestyle=":", linewidth=1, alpha=0.6)

    legend = [
        mpatches.Patch(color=COLOR_TARGET, label="Zielgebiet"),
        mpatches.Patch(color=COLOR_OAR,    label="Risikoorgan"),
    ]
    axes[0].legend(handles=legend, fontsize=8)

    fig.suptitle("Formmetriken aller Strukturen", fontsize=13, fontweight="bold", y=1.01)
    plt.tight_layout()
    path = output_dir / "shape_metrics.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Gespeichert: {path}")


# ---------------------------------------------------------------------------
# Plot 3: Abstands-Uebersicht
# ---------------------------------------------------------------------------

def plot_distances(results: dict, output_dir: Path, max_pairs: int = 25) -> None:
    """
    Gruppiertes Balkendiagramm: Min./Hausdorff-/Schwerpunktabstand pro Strukturpaar.

    Zeigt nur Target-OAR-Paare (klinisch relevant), sortiert nach minimalem Abstand
    (aufsteigend = kritischste zuerst). Auf max_pairs Eintraege begrenzt.
    Klinisch relevante Schwelle (5 mm) als gestrichelte rote Linie eingezeichnet.
    """
    if not results["distances"]:
        return

    target_names = set(results["targets"].keys())
    oar_names    = set(results["oars"].keys())

    # Nur Target-OAR-Paare behalten (klinisch relevant)
    relevant = [
        d for d in results["distances"]
        if (d["structure_a"] in target_names and d["structure_b"] in oar_names)
        or (d["structure_b"] in target_names and d["structure_a"] in oar_names)
    ]

    # Fallback: alle Paare wenn keine Target-OAR-Kombination vorhanden
    if not relevant:
        relevant = results["distances"]

    # Aufsteigend nach Minimalabstand sortieren (kritischste zuerst)
    relevant = sorted(relevant, key=lambda d: d["min_distance_mm"])[:max_pairs]

    pairs = [f"{d['structure_a']}\nvs\n{d['structure_b']}" for d in relevant]
    d_min  = [d["min_distance_mm"]       for d in relevant]
    d_haus = [d["hausdorff_distance_mm"] for d in relevant]
    d_cent = [d["centroid_distance_mm"]  for d in relevant]

    x     = np.arange(len(pairs))
    width = 0.25

    fig, ax = plt.subplots(figsize=(max(8, min(len(pairs) * 2.2, 28)), 6))

    ax.bar(x - width, d_min,  width, label="Min. Abstand",        color=COLOR_MIN,  alpha=0.85)
    ax.bar(x,         d_haus, width, label="Hausdorff-Abstand",   color=COLOR_HAUS, alpha=0.85)
    ax.bar(x + width, d_cent, width, label="Schwerpunkt-Abstand", color=COLOR_CENT, alpha=0.85)

    # Klinische Schwelle: 5 mm fuer Mindestabstand PTV-OAR
    ax.axhline(5.0, color="red", linestyle="--", linewidth=1.2, alpha=0.7,
               label="Klinische Schwelle (5 mm)")

    ax.set_xticks(x)
    ax.set_xticklabels(pairs, fontsize=8)
    ax.set_ylabel("Abstand (mm)", fontsize=11)
    ax.set_title("Abstands-Metriken zwischen Strukturpaaren", fontsize=13,
                 fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(axis="y", alpha=0.3, linestyle="--")
    ax.set_ylim(0, max(d_cent + d_haus + [10]) * 1.1)

    plt.tight_layout()
    path = output_dir / "distances.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Gespeichert: {path}")


# ---------------------------------------------------------------------------
# Plot 4: Schwerpunkt-3D-Karte
# ---------------------------------------------------------------------------

def plot_centroids_3d(results: dict, output_dir: Path) -> None:
    """
    3D-Streudiagramm der Strukturzentren im DICOM-Patientenkoordinatensystem
    (X=links, Y=posterior, Z=superior). Setzt Groesse proportional zum Volumen.
    Gibt raeumliche Lage und Abstands-Verhaeltnisse auf einen Blick wieder.
    """
    fig = plt.figure(figsize=(10, 8))
    ax  = fig.add_subplot(111, projection="3d")

    coords = []   # zum Setzen eines verzerrungsfreien Box-Aspektverhaeltnisses

    def _add_group(struct_dict, color, marker):
        for name, s in struct_dict.items():
            cx, cy, cz = s["centroid_mm"]
            vol = s["volume_cm3"]
            # Marker-Flaeche skaliert mit der Wurzel des Volumens, damit grosse
            # Strukturen sichtbar dominanter, aber kleine nicht unsichtbar sind
            # (lineare Skalierung liesse 100-cm3-Strukturen alles erschlagen).
            size = max(30, min(600, np.sqrt(max(vol, 0.0)) * 40))
            ax.scatter(cx, cy, cz, s=size, c=color, marker=marker,
                       alpha=0.75, edgecolors="white", linewidths=0.5)
            ax.text(cx, cy, cz + 3, name, fontsize=7, ha="center", color=color)
            coords.append((cx, cy, cz))

    _add_group(results["targets"], COLOR_TARGET, "o")
    _add_group(results["oars"],    COLOR_OAR,    "^")

    # Gleiches Laengen-zu-Pixel-Verhaeltnis fuer alle Achsen, sonst werden
    # raeumliche Abstaende im 3D-Plot verzerrt dargestellt (irrefuehrend).
    if coords:
        c_arr = np.array(coords)
        spans = np.ptp(c_arr, axis=0)
        spans[spans == 0] = 1.0
        ax.set_box_aspect(tuple(spans))

    ax.set_xlabel("X  [mm]  (Links)", fontsize=9)
    ax.set_ylabel("Y  [mm]  (Posterior)", fontsize=9)
    ax.set_zlabel("Z  [mm]  (Superior)", fontsize=9)
    ax.set_title("Raeumliche Verteilung der Strukturschwerpunkte",
                 fontsize=12, fontweight="bold")

    legend = [
        mpatches.Patch(color=COLOR_TARGET, label="Zielgebiet (Kreis)"),
        mpatches.Patch(color=COLOR_OAR,    label="Risikoorgan (Dreieck)"),
    ]
    ax.legend(handles=legend, fontsize=9, loc="upper left")

    plt.tight_layout()
    path = output_dir / "centroids_3d.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Gespeichert: {path}")


# ---------------------------------------------------------------------------
# Plot 5: Statistik-Textdatei
# ---------------------------------------------------------------------------

def write_statistics(results: dict, output_dir: Path) -> None:
    """Schreibt eine strukturierte Zusammenfassung aller Metriken als TXT."""

    def _stats(values):
        if not values:
            return {"n": 0, "mean": 0, "std": 0, "min": 0, "max": 0}
        return {
            "n":    len(values),
            "mean": float(np.mean(values)),
            "std":  float(np.std(values)),
            "min":  float(np.min(values)),
            "max":  float(np.max(values)),
        }

    path = output_dir / "statistics.txt"
    with open(path, "w", encoding="utf-8") as f:
        f.write("RTSTRUCT Analyse – Statistiken\n")
        f.write("=" * 50 + "\n\n")

        for section, struct_dict in [("ZIELGEBIETE", results["targets"]),
                                      ("RISIKOORGANE", results["oars"])]:
            f.write(f"{section}\n")
            f.write("-" * 50 + "\n")
            for name, s in struct_dict.items():
                sh = s["shape"]
                f.write(f"\n  {name}  (ROI #{s['roi_number']})\n")
                f.write(f"    Konturen  : {s['num_contours']} Schichten, "
                        f"{s['num_points']} Punkte\n")
                f.write(f"    Volumen   : {s['volume_cm3']:.3f} cm3\n")
                cx, cy, cz = s["centroid_mm"]
                f.write(f"    Schwerpunkt: ({cx:.1f}, {cy:.1f}, {cz:.1f}) mm\n")
                bx, by, bz = sh["bbox_size_mm"]
                f.write(f"    BBox Groesse: {bx:.1f} x {by:.1f} x {bz:.1f} mm\n")
                f.write(f"    Sphaerizitaet : {sh['sphericity']:.4f}\n")
                f.write(f"    Kompaktheit   : {sh['compactness']:.4f}\n")
                f.write(f"    Elongation    : {sh['elongation']:.4f}\n")

        if results["distances"]:
            f.write("\n\nABSTANDS-STATISTIK\n")
            f.write("-" * 50 + "\n")
            st_min  = _stats([d["min_distance_mm"]       for d in results["distances"]])
            st_haus = _stats([d["hausdorff_distance_mm"] for d in results["distances"]])
            st_cent = _stats([d["centroid_distance_mm"]  for d in results["distances"]])

            f.write(f"Anzahl Paare : {st_min['n']}\n\n")
            for label, st in [("Min. Abstand (mm)",        st_min),
                               ("Hausdorff-Abstand (mm)",  st_haus),
                               ("Schwerpunkt-Abstand (mm)", st_cent)]:
                f.write(f"  {label}\n")
                f.write(f"    Mittelwert : {st['mean']:.2f}\n")
                f.write(f"    Std.abw.   : {st['std']:.2f}\n")
                f.write(f"    Min        : {st['min']:.2f}\n")
                f.write(f"    Max        : {st['max']:.2f}\n\n")

            f.write("\nEINZELPAARE\n")
            f.write(f"{'Paar':<45} {'Min':>8} {'Hausdorff':>12} {'Zentroid':>10}\n")
            f.write("-" * 78 + "\n")
            for d in results["distances"]:
                pair = f"{d['structure_a']} vs {d['structure_b']}"
                f.write(f"{pair:<45} {d['min_distance_mm']:>8.2f} "
                        f"{d['hausdorff_distance_mm']:>12.2f} "
                        f"{d['centroid_distance_mm']:>10.2f}\n")

    print(f"  Gespeichert: {path}")


# ===========================================================================
# Case-Transform-Visualisierung  (Vorher/Nachher einer rigiden Transformation)
# ===========================================================================
#
# Diese Funktionen werden von ``case_modifier`` nach einem Transform-Lauf
# aufgerufen.  Sie vergleichen die *Original*-RTSTRUCT-Geometrie mit der
# *transformierten* und blenden das gewaehlte Rotationszentrum sowie die
# POINT-Marker ein.  Es werden nur die Konturpunkte benoetigt -- keine
# CT-Pixeldaten -- daher sind die Plots auch im ``metadata``-Modus guenstig.
# Die (teure) CT-Koerperoberflaeche ist optional zuschaltbar.


def _structure_pointclouds(
    ds, subsample_per_roi: "int | None" = 600
) -> dict:
    """
    Sammelt pro nicht-POINT-ROI die Konturpunkte als (N,3)-Array in LPS-mm.

    POINT-Typ-Konturen (Marker, Drehpunkt) werden uebersprungen -- diese werden
    separat behandelt.  Bei ``subsample_per_roi`` wird je ROI gleichmaessig
    (deterministisch) ausgeduennt, damit die Plots auch bei vielen ROIs
    reaktionsschnell bleiben.
    """
    names = get_structure_names(ds)
    clouds: dict = {}
    if not hasattr(ds, "ROIContourSequence"):
        return clouds

    for rc in ds.ROIContourSequence:
        if not hasattr(rc, "ContourSequence"):
            continue
        roi_num = int(getattr(rc, "ReferencedROINumber", -1))
        pts_list = []
        is_point = False
        for c in rc.ContourSequence:
            if str(getattr(c, "ContourGeometricType", "")) == "POINT":
                is_point = True
                break
            if hasattr(c, "ContourData"):
                arr = np.array(c.ContourData, dtype=np.float64).reshape(-1, 3)
                if arr.size:
                    pts_list.append(arr)
        if is_point or not pts_list:
            continue

        pts = np.vstack(pts_list)
        if subsample_per_roi and len(pts) > subsample_per_roi:
            idx = np.linspace(0, len(pts) - 1, subsample_per_roi).astype(int)
            pts = pts[idx]
        clouds[names.get(roi_num, f"ROI#{roi_num}")] = pts

    return clouds


def _stack_clouds(clouds: dict) -> np.ndarray:
    """Fasst alle ROI-Punktwolken zu einem (N,3)-Array zusammen (ggf. leer)."""
    if not clouds:
        return np.empty((0, 3))
    return np.vstack(list(clouds.values()))


# ---------------------------------------------------------------------------
# A. Statische tri-planare Vorher/Nachher-Projektion
# ---------------------------------------------------------------------------

def plot_transform_overview(
    orig_ds, new_ds,
    center: np.ndarray,
    drehpunkt_pos: np.ndarray,
    markers: "list | None",
    output_dir: Path,
) -> None:
    """
    Drei orthografische Projektionen (axial X-Y, koronal X-Z, sagittal Y-Z) der
    Strukturpunkte: Original (grau) ueberlagert mit Transformiert (rot).
    Rotationszentrum als Stern, Translationsvektor als Pfeil, POINT-Marker als x.

    Statisches PNG fuer Reports/CI -- zeigt auf einen Blick Richtung und Betrag
    der Verschiebung sowie eine evtl. Rotation (verkippte Punktwolke).
    """
    orig = _stack_clouds(_structure_pointclouds(orig_ds))
    new  = _stack_clouds(_structure_pointclouds(new_ds))
    if orig.size == 0 and new.size == 0:
        return

    center        = np.asarray(center, dtype=np.float64).reshape(3)
    drehpunkt_pos = np.asarray(drehpunkt_pos, dtype=np.float64).reshape(3)

    # (Index a, Index b, Titel, X-Achsen-Label, Y-Achsen-Label, Y invertieren?)
    planes = [
        (0, 1, "Axial (X-Y)",    "X [mm] (Links)",     "Y [mm] (Posterior)", True),
        (0, 2, "Koronal (X-Z)",  "X [mm] (Links)",     "Z [mm] (Superior)",  False),
        (1, 2, "Sagittal (Y-Z)", "Y [mm] (Posterior)", "Z [mm] (Superior)",  False),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    for ax, (a, b, title, xl, yl, inv_y) in zip(axes, planes):
        if orig.size:
            ax.scatter(orig[:, a], orig[:, b], s=2, c=COLOR_ORIG,
                       alpha=0.25, linewidths=0, label="Original")
        if new.size:
            ax.scatter(new[:, a], new[:, b], s=2, c=COLOR_TRANS,
                       alpha=0.30, linewidths=0, label="Transformiert")

        # Translationsvektor: Rotationszentrum -> Drehpunkt (= Zentrum + t)
        ax.annotate(
            "", xy=(drehpunkt_pos[a], drehpunkt_pos[b]),
            xytext=(center[a], center[b]),
            arrowprops=dict(arrowstyle="->", color="black", lw=1.6),
        )
        ax.scatter([center[a]], [center[b]], marker="*", s=240, c=COLOR_CENTER,
                   edgecolors="black", linewidths=0.6, zorder=6,
                   label="Rotationszentrum")
        ax.scatter([drehpunkt_pos[a]], [drehpunkt_pos[b]], marker="P", s=90,
                   c=COLOR_CENTER, edgecolors="black", linewidths=0.6, zorder=6,
                   label="Drehpunkt (transformiert)")

        if markers:
            mx = [p[a] for _, p in markers]
            my = [p[b] for _, p in markers]
            ax.scatter(mx, my, marker="x", s=30, c=COLOR_MARKER, linewidths=1.0,
                       zorder=5, label="POINT-Marker")

        ax.set_xlabel(xl, fontsize=9)
        ax.set_ylabel(yl, fontsize=9)
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.set_aspect("equal", adjustable="datalim")
        ax.grid(alpha=0.25, linestyle="--")
        if inv_y:
            ax.invert_yaxis()

    # Eine gemeinsame Legende (Eintraege aus dem ersten Panel reichen)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=len(labels),
               fontsize=9, frameon=True)
    fig.suptitle("Rigide Transformation -- Vorher (grau) / Nachher (rot)",
                 fontsize=13, fontweight="bold")
    plt.tight_layout(rect=(0, 0.05, 1, 0.97))
    path = output_dir / "transform_overview.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Gespeichert: {path}")


# ---------------------------------------------------------------------------
# B. Schwerpunkt-Verschiebung pro ROI
# ---------------------------------------------------------------------------

def plot_displacement_summary(
    orig_ds,
    T: np.ndarray,
    translation: "tuple | np.ndarray",
    output_dir: Path,
    max_rois: int = 40,
) -> None:
    """
    Balkendiagramm der Schwerpunkt-Verschiebung ``|T(c) - c|`` pro ROI (mm),
    absteigend sortiert.  Eine gestrichelte Linie markiert den reinen
    Translationsbetrag ``|t|``: Strukturen genau im Rotationszentrum wandern um
    ``|t|``, weiter entfernte zusaetzlich durch den Rotations-Hebelarm.

    Reines QA-Diagramm -- macht sichtbar, wie stark die Rotation (nicht nur die
    Translation) einzelne Strukturen bewegt.
    """
    clouds = _structure_pointclouds(orig_ds, subsample_per_roi=None)
    if not clouds:
        return

    rows = []
    for name, pts in clouds.items():
        c = pts.mean(axis=0)
        c_new = (T @ np.append(c, 1.0))[:3]
        rows.append((name, float(np.linalg.norm(c_new - c))))

    rows.sort(key=lambda r: r[1], reverse=True)
    clipped = len(rows) > max_rois
    rows = rows[:max_rois]
    names = [r[0] for r in rows]
    disp  = [r[1] for r in rows]

    t_norm = float(np.linalg.norm(np.asarray(translation, dtype=np.float64)))

    fig, ax = plt.subplots(figsize=(10, max(4, len(names) * 0.32 + 1)))
    bars = ax.barh(names, disp, color=COLOR_TRANS, edgecolor="white", linewidth=0.5)
    for bar, d in zip(bars, disp):
        ax.text(bar.get_width() + 0.01 * max(disp + [1]),
                bar.get_y() + bar.get_height() / 2,
                f"{d:.1f}", va="center", ha="left", fontsize=8)

    if t_norm > 0:
        ax.axvline(t_norm, color="black", linestyle="--", linewidth=1.2,
                   alpha=0.7, label=f"|Translation| = {t_norm:.1f} mm")
        ax.legend(loc="lower right", fontsize=9)

    title = "Schwerpunkt-Verschiebung pro ROI"
    if clipped:
        title += f"  (Top {max_rois})"
    ax.set_xlabel("Verschiebung |T(c) - c|  [mm]", fontsize=11)
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.set_xlim(0, max(disp + [1]) * 1.15)
    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.3, linestyle="--")

    plt.tight_layout()
    path = output_dir / "displacement.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Gespeichert: {path}")


# ---------------------------------------------------------------------------
# C. Interaktiver 3D-Vergleich (Plotly)
# ---------------------------------------------------------------------------

def plot_transform_3d(
    orig_ds, new_ds,
    center: np.ndarray,
    drehpunkt_pos: np.ndarray,
    markers: "list | None",
    output_dir: Path,
    T: "np.ndarray | None" = None,
    geom: "dict | None" = None,
    volume_hu: "np.ndarray | None" = None,
    ct_surface: bool = False,
    max_points: int = 9000,
) -> None:
    """
    Interaktiver 3D-Vergleich (Plotly HTML) der Original- vs. transformierten
    Konturgeometrie.  Jede Ebene ist ueber die Legende einzeln zuschaltbar:

      - "Konturen (Original)"        -- ausgeduennte Punktwolke, grau-blau
      - "Konturen (Transformiert)"   -- ausgeduennte Punktwolke, rot
      - "POINT-Marker (Original)"    -- alle POINT-ROIs mit Namen
      - "POINT-Marker (Transformiert)"
      - "Drehpunkt / Rotationszentrum" -- Stern + Translationsvektor
      - "DICOM-Achsen"               -- L/P/S-Triade
      - "CT-Koerperoberflaeche"      -- nur wenn ``ct_surface`` und Volumen
                                        vorhanden (per Default ausgeblendet)

    Die Toggles fuer Rotationszentrum und Marker erfuellen die Anforderung,
    diese Referenzen vor/nach der Transformation ein- und ausblenden zu koennen.
    """
    try:
        import plotly.graph_objects as go
    except ImportError:
        print("  Hinweis: plotly nicht installiert -- 3D-HTML uebersprungen.")
        return

    center        = np.asarray(center, dtype=np.float64).reshape(3)
    drehpunkt_pos = np.asarray(drehpunkt_pos, dtype=np.float64).reshape(3)

    orig = _stack_clouds(_structure_pointclouds(orig_ds))
    new  = _stack_clouds(_structure_pointclouds(new_ds))

    def _thin(pts):
        if len(pts) > max_points:
            idx = np.linspace(0, len(pts) - 1, max_points).astype(int)
            return pts[idx]
        return pts

    orig = _thin(orig)
    new  = _thin(new)

    fig = go.Figure()

    # optionale CT-Koerperoberflaeche (teuer -> nur auf Wunsch)
    if ct_surface and volume_hu is not None and geom is not None:
        try:
            from . import modifier as mod
            print("  Extrahiere CT-Koerperoberflaeche (Original) …")
            verts, faces = mod._extract_surface(
                volume_hu, geom["affine"], threshold=-300.0, downsample=3
            )
            if verts is not None:
                fig.add_trace(go.Mesh3d(
                    x=verts[:, 0], y=verts[:, 1], z=verts[:, 2],
                    i=faces[:, 0], j=faces[:, 1], k=faces[:, 2],
                    color="royalblue", opacity=0.12,
                    name="CT-Koerperoberflaeche (Original)",
                    showlegend=True, visible="legendonly",
                    lighting=dict(diffuse=0.8, roughness=0.6),
                ))
        except Exception as e:   # Visualisierung darf den Lauf nie kippen
            print(f"  Hinweis: CT-Oberflaeche uebersprungen ({e}).")

    if orig.size:
        fig.add_trace(go.Scatter3d(
            x=orig[:, 0], y=orig[:, 1], z=orig[:, 2], mode="markers",
            marker=dict(size=1.6, color=COLOR_ORIG, opacity=0.35),
            name="Konturen (Original)",
        ))
    if new.size:
        fig.add_trace(go.Scatter3d(
            x=new[:, 0], y=new[:, 1], z=new[:, 2], mode="markers",
            marker=dict(size=1.6, color=COLOR_TRANS, opacity=0.5),
            name="Konturen (Transformiert)",
        ))

    # POINT-Marker (Original + transformiert via T)
    if markers:
        m_orig = np.array([p for _, p in markers], dtype=np.float64)
        m_names = [n for n, _ in markers]
        fig.add_trace(go.Scatter3d(
            x=m_orig[:, 0], y=m_orig[:, 1], z=m_orig[:, 2],
            mode="markers+text", text=m_names, textfont=dict(size=8),
            marker=dict(size=4, color=COLOR_MARKER, symbol="diamond"),
            name="POINT-Marker (Original)", visible="legendonly",
        ))
        if T is not None:
            homog = np.hstack([m_orig, np.ones((len(m_orig), 1))])
            m_new = (np.asarray(T) @ homog.T).T[:, :3]
            fig.add_trace(go.Scatter3d(
                x=m_new[:, 0], y=m_new[:, 1], z=m_new[:, 2],
                mode="markers+text", text=m_names, textfont=dict(size=8),
                marker=dict(size=4, color=COLOR_TRANS, symbol="diamond"),
                name="POINT-Marker (Transformiert)", visible="legendonly",
            ))

    fig.add_trace(go.Scatter3d(
        x=[center[0]], y=[center[1]], z=[center[2]], mode="markers+text",
        text=["Rotationszentrum"], textfont=dict(size=10, color=COLOR_CENTER),
        marker=dict(size=7, color=COLOR_CENTER, symbol="x",
                    line=dict(color="black", width=2)),
        name="Drehpunkt / Rotationszentrum",
    ))
    # Translationsvektor c -> c + t
    fig.add_trace(go.Scatter3d(
        x=[center[0], drehpunkt_pos[0]],
        y=[center[1], drehpunkt_pos[1]],
        z=[center[2], drehpunkt_pos[2]],
        mode="lines", line=dict(color=COLOR_CENTER, width=5, dash="dash"),
        name="Translationsvektor",
    ))

    # DICOM-Achsentriade am Rotationszentrum
    axis_len = 60.0
    for vec, color, label in [
        (np.array([1, 0, 0]), "red",  "X (L)"),
        (np.array([0, 1, 0]), "lime", "Y (P)"),
        (np.array([0, 0, 1]), "cyan", "Z (S)"),
    ]:
        end = center + vec * axis_len
        fig.add_trace(go.Scatter3d(
            x=[center[0], end[0]], y=[center[1], end[1]], z=[center[2], end[2]],
            mode="lines+text", line=dict(color=color, width=4),
            text=["", label], textposition="top center",
            textfont=dict(color=color, size=11),
            name=f"Achse {label}", showlegend=False, visible="legendonly",
        ))

    fig.update_layout(
        title=dict(
            text="RTSTRUCT Rigid Body Transform -- Original vs. Transformiert",
            font=dict(size=15)),
        scene=dict(
            xaxis_title="X [mm] (Links)",
            yaxis_title="Y [mm] (Posterior)",
            zaxis_title="Z [mm] (Superior)",
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
        margin=dict(l=0, r=0, b=0, t=60),
    )

    path = output_dir / "transform_3d.html"
    fig.write_html(str(path))
    print(f"  Gespeichert: {path}")


def run_case_visualization(
    orig_ds, new_ds,
    center: np.ndarray,
    drehpunkt_pos: np.ndarray,
    translation: "tuple | np.ndarray",
    T: np.ndarray,
    output_dir: Path,
    markers: "list | None" = None,
    geom: "dict | None" = None,
    volume_hu: "np.ndarray | None" = None,
    ct_surface: bool = False,
) -> None:
    """
    Erzeugt alle Case-Transform-Plots in ``output_dir``.

    ``orig_ds`` / ``new_ds`` sind das RTSTRUCT vor bzw. nach der Transformation,
    ``markers`` ist die ``[(name, position_lps), ...]``-Liste der POINT-Marker
    aus dem Original-RS (transformierte Positionen werden via ``T`` berechnet).
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"\nErstelle Case-Transform-Visualisierungen in {output_dir} …")
    plot_transform_overview(orig_ds, new_ds, center, drehpunkt_pos,
                            markers, output_dir)
    plot_displacement_summary(orig_ds, T, translation, output_dir)
    plot_transform_3d(orig_ds, new_ds, center, drehpunkt_pos, markers,
                      output_dir, T=T, geom=geom, volume_hu=volume_hu,
                      ct_surface=ct_surface)
    print("Case-Visualisierung abgeschlossen.")


# ---------------------------------------------------------------------------
# Haupt-Workflow
# ---------------------------------------------------------------------------

def run_visualization(results: dict, output_dir: Path) -> None:
    """Erstellt alle Plots fuer ein bereits analysiertes results-Dict."""
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"\nErstelle Visualisierungen in {output_dir} ...")
    plot_volumes(results, output_dir)
    plot_shape_metrics(results, output_dir)
    plot_distances(results, output_dir)
    plot_centroids_3d(results, output_dir)
    write_statistics(results, output_dir)
    print("Visualisierung abgeschlossen.")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="RTSTRUCT Visualizer – Plots aus DICOM RT Structure Sets",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Beispiele:\n"
            "  python -m dicom_file_modifier.visualizer data/0000000171/test/1.dcm\n"
            "  python -m dicom_file_modifier.visualizer data/rtstruct.dcm "
            "--targets PTV,CTV --oars Parotis --output output/plots\n"
        ),
    )
    parser.add_argument("file", help="Pfad zur RTSTRUCT DICOM Datei")
    parser.add_argument("--output", "-o", default="output",
                        help="Ausgabeverzeichnis (Standard: output)")
    parser.add_argument("--targets", type=str, default=None,
                        help="Komma-getrennte Zielgebiet-Namen (z.B. PTV,CTV)")
    parser.add_argument("--oars", type=str, default=None,
                        help="Komma-getrennte Risikoorgan-Namen")
    args = parser.parse_args()

    target_list = args.targets.split(",") if args.targets else None
    oar_list    = args.oars.split(",")    if args.oars    else None

    results = run_analysis(
        filepath=args.file,
        target_names=target_list,
        oar_names=oar_list,
    )

    run_visualization(results, Path(args.output))


if __name__ == "__main__":
    main()
