#!/usr/bin/env python3
"""
RTSTRUCT Visualizer
===================
Erzeugt aussagekraeftige Plots aus RTSTRUCT-Analyseergebnissen.

Plots:
  1. Volumen-Balkendiagramm  – alle Strukturen, farbcodiert nach Typ (Target/OAR)
  2. Formmetriken-Vergleich  – Sphaerizitaet, Kompaktheit, Elongation nebeneinander
  3. Abstands-Uebersicht     – Min/Hausdorff/Zentroid pro Strukturpaar als Balken
  4. Schwerpunkt-3D-Karte    – raeumliche Positionen aller Strukturzentren im Patientenraum
  5. statistics.txt          – Zahlenzusammenfassung aller Metriken

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

from dicom_file_modifier.analyzer import run_analysis


# Farben
COLOR_TARGET = "#2471A3"   # Blau fuer Zielgebiete
COLOR_OAR    = "#CB4335"   # Rot fuer Risikoorgane
COLOR_MIN    = "#27AE60"   # Gruen fuer Minimalabstand
COLOR_HAUS   = "#E67E22"   # Orange fuer Hausdorff
COLOR_CENT   = "#8E44AD"   # Lila fuer Schwerpunktabstand


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

    def _add_group(struct_dict, color, marker):
        for name, s in struct_dict.items():
            cx, cy, cz = s["centroid_mm"]
            vol = s["volume_cm3"]
            size = max(30, min(600, vol * 8))     # Groesse skaliert mit Volumen
            ax.scatter(cx, cy, cz, s=size, c=color, marker=marker,
                       alpha=0.75, edgecolors="white", linewidths=0.5)
            ax.text(cx, cy, cz + 3, name, fontsize=7, ha="center", color=color)

    _add_group(results["targets"], COLOR_TARGET, "o")
    _add_group(results["oars"],    COLOR_OAR,    "^")

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
