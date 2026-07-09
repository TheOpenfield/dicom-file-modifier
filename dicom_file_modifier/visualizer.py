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

from dicom_file_modifier.analyzer import (
    run_analysis, get_structure_names,
    CAT_TARGET, CAT_OAR_SERIAL, CAT_OAR_PARALLEL, CAT_HELPER,
    CATEGORY_ORDER,
)


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

# Einheitliche Kategorie-Farben/-Beschriftungen fuer ALLE Einzel-RTSTRUCT-Plots
CATEGORY_COLORS = {
    CAT_TARGET:       "#2471A3",   # Blau   – Zielvolumina (GTV/PTV)
    CAT_OAR_SERIAL:   "#C0392B",   # Rot    – serielle OARs (max-dosis-kritisch)
    CAT_OAR_PARALLEL: "#E67E22",   # Orange – parallele OARs (Volumeneffekt)
    CAT_HELPER:       "#7F8C8D",   # Grau   – Hilfs-/Planungsstrukturen
}
CATEGORY_LABELS = {
    CAT_TARGET:       "Zielvolumen (GTV/PTV)",
    CAT_OAR_SERIAL:   "OAR seriell",
    CAT_OAR_PARALLEL: "OAR parallel",
    CAT_HELPER:       "Hilfs-/Planungsstruktur",
}


def _iter_structures(results: dict):
    """Alle ausgewerteten Strukturen (Targets, OARs, Helpers) als (name, s)."""
    for bucket in ("targets", "oars", "helpers"):
        for name, s in results.get(bucket, {}).items():
            yield name, s


def _category_sort_key(s: dict):
    """Sortierschlüssel: Kategorie-Reihenfolge, dann (für Targets) Läsion."""
    cat = s.get("category", CAT_HELPER)
    ci = CATEGORY_ORDER.index(cat) if cat in CATEGORY_ORDER else len(CATEGORY_ORDER)
    if cat == CAT_TARGET:
        # GTV vor PTV innerhalb derselben Läsion
        key = s.get("lesion_key") or s["name"]
        prefix = 0 if s["name"].upper().startswith("GTV") else 1
        return (ci, key, prefix)
    return (ci, -s.get("volume_cm3", 0.0), s["name"])


# ---------------------------------------------------------------------------
# Plot 1: Volumen-Balkendiagramm
# ---------------------------------------------------------------------------

def plot_volumes(results: dict, output_dir: Path) -> None:
    """
    Horizontales Balkendiagramm aller Strukturvolumen auf LOG-Skala,
    farbcodiert nach Kategorie und nach Kategorie gruppiert (Targets als
    GTV/PTV-Paare benachbart, dann serielle/parallele OARs, dann – schraffiert
    und abgesetzt – Hilfsstrukturen). Eine lineare Skala würde die kleinen
    SRS-Targets (0.03 cm³) neben dem Ganzhirn (>1000 cm³) unsichtbar machen.
    """
    rows = [(name, s) for name, s in _iter_structures(results)
            if s.get("volume_cm3", 0) > 0]
    if not rows:
        return
    rows.sort(key=lambda r: _category_sort_key(r[1]))

    names   = [r[0] for r in rows]
    volumes = [r[1]["volume_cm3"] for r in rows]
    cats    = [r[1].get("category", CAT_HELPER) for r in rows]
    colors  = [CATEGORY_COLORS.get(c, "#95A5A6") for c in cats]
    hatches = ["///" if c == CAT_HELPER else "" for c in cats]

    y = np.arange(len(names))
    fig, ax = plt.subplots(figsize=(10, max(4, len(names) * 0.34 + 1.5)))
    vmin = max(0.01, min(volumes) * 0.5)
    # Balken vom Log-Achsenboden (vmin) bis zum echten Volumen. NICHT left=vmin
    # + width=volume verwenden: dort wäre die rechte Kante vmin+volume und würde
    # kleine Targets auf der Log-Skala überzeichnen. Breite = volume - vmin
    # lässt die Spitze exakt auf dem Volumen landen (vmin < jedes Volumen).
    bars = ax.barh(y, np.array(volumes) - vmin, color=colors, edgecolor="white",
                   linewidth=0.5, left=vmin)
    for bar, hatch in zip(bars, hatches):
        if hatch:
            bar.set_hatch(hatch)

    for yi, vol in zip(y, volumes):
        ax.text(vol * 1.05, yi, f"{vol:.2f}" if vol < 10 else f"{vol:.0f}",
                va="center", ha="left", fontsize=7)

    ax.set_xscale("log")
    ax.set_xlim(vmin, max(volumes) * 2.2)
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Volumen (cm³, log-Skala)", fontsize=11)
    ax.set_title("Strukturvolumen nach Kategorie", fontsize=13, fontweight="bold")
    ax.grid(axis="x", alpha=0.3, linestyle="--", which="both")

    present = [c for c in CATEGORY_ORDER if c in cats]
    legend = [mpatches.Patch(facecolor=CATEGORY_COLORS.get(c, "#95A5A6"),
                             hatch="///" if c == CAT_HELPER else None,
                             edgecolor="white", label=CATEGORY_LABELS.get(c, c))
              for c in present]
    ext = results.get("meta", {}).get("external_names", [])
    if ext:
        legend.append(mpatches.Patch(facecolor="none", edgecolor="none",
                                     label=f"(External ausgeblendet: {', '.join(ext)})"))
    ax.legend(handles=legend, loc="lower right", fontsize=8)

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
    Heatmap-Tabelle der Formmetriken (Sphärizität, Solidität, Elongation) –
    NUR für anatomische Einzelkomponenten-Strukturen (Targets + OARs), nach
    Kategorie gruppiert. Hilfs-/Vereinigungsstrukturen werden ausgelassen, da
    Hüllen-/Voxel-Formmetriken für Mehrkomponenten-Strukturen bedeutungslos
    sind. Eine Heatmap ersetzt die unleserlichen 40 gedrehten x-Labels und
    lässt Ausreißer (z.B. Rückenmark-Elongation) sofort erkennen.
    """
    import matplotlib.colors as mcolors

    # Nur anatomische Einzelkomponenten-Strukturen mit gültiger Formmetrik
    # (gleiche Auswahl wie plot_sphericity_vs_elongation -> konsistente Plots;
    # ungültige/Mehrkomponenten-Strukturen erscheinen weiter in statistics.txt).
    rows = [(name, s) for name, s in _iter_structures(results)
            if s.get("category") in (CAT_TARGET, CAT_OAR_SERIAL, CAT_OAR_PARALLEL)
            and s["shape"].get("shape_valid", False)]
    if not rows:
        return
    rows.sort(key=lambda r: _category_sort_key(r[1]))

    names = [r[0] for r in rows]
    cats  = [r[1].get("category") for r in rows]
    sph = np.array([r[1]["shape"]["sphericity"]  for r in rows])
    sol = np.array([r[1]["shape"]["solidity"]    for r in rows])
    elo = np.array([r[1]["shape"]["elongation"]  for r in rows])
    n = len(names)

    fig, axes = plt.subplots(
        1, 4, figsize=(11, max(4, n * 0.32 + 1.5)),
        gridspec_kw={"width_ratios": [0.12, 1, 1, 1], "wspace": 0.08})
    y = np.arange(n)

    # Spalte 0: Kategorie-Farbstreifen
    cax = axes[0]
    cat_rgb = [mcolors.to_rgb(CATEGORY_COLORS.get(c, "#95A5A6")) for c in cats]
    cax.imshow(np.array(cat_rgb).reshape(n, 1, 3), aspect="auto")
    cax.set_xticks([])
    cax.set_yticks(y)
    cax.set_yticklabels(names, fontsize=7)
    cax.set_title("Kat.", fontsize=8)

    def _strip(ax, vals, title, cmap, vmin, vmax, fmt="{:.2f}"):
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
        ax.imshow(vals.reshape(n, 1), aspect="auto", cmap=cmap, norm=norm)
        for yi, v in zip(y, vals):
            ax.text(0, yi, fmt.format(v), ha="center", va="center", fontsize=7,
                    color="black")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(title, fontsize=9, fontweight="bold")

    # Sphärizität & Solidität: 0..1, grün = rund/konvex (gut)
    _strip(axes[1], sph, "Sphärizität\n(1=Kugel)", "RdYlGn", 0.0, 1.0)
    _strip(axes[2], sol, "Solidität\n(1=konvex)", "RdYlGn", 0.0, 1.0)
    # Elongation: ab 1; höher = gestreckter -> Orangetöne
    _strip(axes[3], elo, "Elongation\n(1=isotrop)", "Oranges",
           1.0, float(max(elo.max(), 2.0)), fmt="{:.1f}")

    # Kategorie-Trennlinien
    bounds = [i for i in range(1, n) if cats[i] != cats[i - 1]]
    for ax in axes:
        for b in bounds:
            ax.axhline(b - 0.5, color="black", linewidth=1.0, alpha=0.6)

    present = [c for c in CATEGORY_ORDER if c in cats]
    legend = [mpatches.Patch(color=CATEGORY_COLORS.get(c, "#95A5A6"),
                             label=CATEGORY_LABELS.get(c, c)) for c in present]
    fig.legend(handles=legend, loc="lower center", ncol=len(legend),
               fontsize=8, frameon=True, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Formmetriken (anatomische Strukturen; Hilfsstrukturen ausgelassen)",
                 fontsize=12, fontweight="bold", y=1.0)
    plt.tight_layout(rect=(0, 0.03, 1, 0.98))
    path = output_dir / "shape_metrics.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Gespeichert: {path}")


# ---------------------------------------------------------------------------
# Plot 3: Abstands-Uebersicht
# ---------------------------------------------------------------------------

def plot_distances(results: dict, output_dir: Path, max_pairs: int = 25) -> None:
    """
    Lollipop-Diagramm der klinisch kritischen Abstände: Target ↔ serielle OARs
    (Hirnstamm, Rückenmark, Sehnerven, Chiasma, Hypophyse), aufsteigend nach
    Minimalabstand. Punkt = Min-Abstand, offener Marker = HD95; Linie zwischen
    beiden. Schwellen bei 3 mm / 5 mm. Ersetzt das alte Balkendiagramm, das von
    Containment-Artefakten (PTV in eigener Vereinigung / im Ganzhirn) dominiert
    wurde – diese sind hier durch die Kategorie-Filterung ausgeschlossen.
    """
    dist = results.get("distances", [])
    if not dist:
        return

    serial = [d for d in dist
              if d.get("pair_type") == "target-oar" and d.get("oar_subtype") == "serial"]
    subtitle = "Target ↔ serielle OARs"
    if not serial:   # Fallback: alle Target-OAR-Paare
        serial = [d for d in dist if d.get("pair_type") == "target-oar"]
        subtitle = "Target ↔ OAR"
    if not serial:
        return

    serial = sorted(serial, key=lambda d: d["min_distance_mm"])[:max_pairs]
    labels = [f"{d['structure_a']} → {d['structure_b']}" for d in serial]
    dmin   = [d["min_distance_mm"] for d in serial]
    dh95   = [d.get("hd95_mm", d["hausdorff_distance_mm"]) for d in serial]

    y = np.arange(len(serial))[::-1]   # kleinster Abstand oben
    fig, ax = plt.subplots(figsize=(10, max(4, len(serial) * 0.34 + 1.5)))

    for yi, mn, h95 in zip(y, dmin, dh95):
        ax.plot([mn, h95], [yi, yi], color="#BDC3C7", linewidth=1.5, zorder=1)
    flagged = [mn < 5.0 for mn in dmin]
    ax.scatter(dmin, y, s=55, zorder=3, label="Min. Abstand",
               color=[("#C0392B" if f else COLOR_MIN) for f in flagged],
               edgecolors="black", linewidths=0.5)
    ax.scatter(dh95, y, s=45, zorder=2, facecolors="none",
               edgecolors=COLOR_HAUS, linewidths=1.4, label="HD95")

    ax.axvline(5.0, color="red", linestyle="--", linewidth=1.2, alpha=0.7,
               label="Schwelle 5 mm")
    ax.axvline(3.0, color="darkred", linestyle=":", linewidth=1.2, alpha=0.7,
               label="Schwelle 3 mm")

    for yi, mn in zip(y, dmin):
        ax.text(mn, yi + 0.28, f"{mn:.1f}", ha="center", va="bottom", fontsize=7,
                color="black")

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Abstand (mm)", fontsize=11)
    ax.set_title(f"Kritische Abstände: {subtitle}\n(aufsteigend, < 5 mm rot markiert)",
                 fontsize=12, fontweight="bold")
    ax.set_xlim(left=0)
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(axis="x", alpha=0.3, linestyle="--")

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
    3D-Streudiagramm der Strukturzentren (LPS: X=links, Y=posterior, Z=superior),
    farbcodiert nach Kategorie, Markergröße ∝ √Volumen. Statt aller ~40 Labels
    (vorher unleserlicher Label-Brei) werden nur serielle OARs beschriftet und
    die Targets durchnummeriert (Zuordnung Nummer→Läsion in der Seitenlegende).
    Marker und External sind ausgeschlossen.
    """
    fig = plt.figure(figsize=(11, 8))
    ax  = fig.add_subplot(111, projection="3d")
    coords = []

    # Targets durchnummerieren (nach Läsion, GTV/PTV teilen sich eine Nummer)
    target_rows = sorted(results.get("targets", {}).items(),
                         key=lambda kv: _category_sort_key(kv[1]))
    lesion_to_num: dict = {}
    num_legend = []
    for name, s in target_rows:
        key = s.get("lesion_key") or name
        if key not in lesion_to_num:
            lesion_to_num[key] = len(lesion_to_num) + 1
            num_legend.append(f"{lesion_to_num[key]}: {key}")
        cx, cy, cz = s["centroid_mm"]
        size = max(25, min(500, np.sqrt(max(s["volume_cm3"], 0.0)) * 40))
        ax.scatter(cx, cy, cz, s=size, c=CATEGORY_COLORS[CAT_TARGET],
                   marker="o", alpha=0.8, edgecolors="white", linewidths=0.5)
        ax.text(cx, cy, cz + 2, str(lesion_to_num[key]), fontsize=7,
                ha="center", color="black", fontweight="bold")
        coords.append((cx, cy, cz))

    # OARs: seriell (Dreieck, beschriftet) + parallel (Quadrat, unbeschriftet)
    for name, s in results.get("oars", {}).items():
        sub = s.get("oar_subtype")
        col = CATEGORY_COLORS[CAT_OAR_SERIAL if sub == "serial" else CAT_OAR_PARALLEL]
        marker = "^" if sub == "serial" else "s"
        cx, cy, cz = s["centroid_mm"]
        size = max(25, min(500, np.sqrt(max(s["volume_cm3"], 0.0)) * 40))
        ax.scatter(cx, cy, cz, s=size, c=col, marker=marker,
                   alpha=0.75, edgecolors="white", linewidths=0.5)
        if sub == "serial":
            ax.text(cx, cy, cz + 3, name, fontsize=7, ha="center", color=col)
        coords.append((cx, cy, cz))

    if coords:
        spans = np.ptp(np.array(coords), axis=0)
        spans[spans == 0] = 1.0
        ax.set_box_aspect(tuple(spans))

    ax.set_xlabel("X [mm] (Links)", fontsize=9)
    ax.set_ylabel("Y [mm] (Posterior)", fontsize=9)
    ax.set_zlabel("Z [mm] (Superior)", fontsize=9)
    ax.set_title("Räumliche Verteilung der Strukturschwerpunkte",
                 fontsize=12, fontweight="bold")

    legend = [
        mpatches.Patch(color=CATEGORY_COLORS[CAT_TARGET], label="Target (Kreis, nummeriert)"),
        mpatches.Patch(color=CATEGORY_COLORS[CAT_OAR_SERIAL], label="OAR seriell (Dreieck)"),
        mpatches.Patch(color=CATEGORY_COLORS[CAT_OAR_PARALLEL], label="OAR parallel (Quadrat)"),
    ]
    ax.legend(handles=legend, fontsize=8, loc="upper left")
    if num_legend:
        fig.text(0.015, 0.5, "Targets:\n" + "\n".join(num_legend),
                 fontsize=7, va="center", ha="left",
                 bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

    plt.tight_layout(rect=(0.16, 0, 1, 1))
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
        return {"n": len(values), "mean": float(np.mean(values)),
                "std": float(np.std(values)), "min": float(np.min(values)),
                "max": float(np.max(values))}

    def _write_struct(f, name, s):
        sh = s["shape"]
        flag = "" if sh.get("shape_valid", False) else "   [!] Mehrkomponenten/ungueltige Formmetrik"
        f.write(f"\n  {name}  (ROI #{s['roi_number']}, {s.get('category','?')})\n")
        f.write(f"    Konturen   : {s['num_contours']} Schichten, {s['num_points']} Punkte\n")
        f.write(f"    Volumen    : {s['volume_cm3']:.3f} cm3"
                f"   (Voxel-Quervergleich: {sh.get('volume_voxel_cm3', 0):.3f} cm3)\n")
        f.write(f"    Aequiv.-Durchmesser: {sh.get('equivalent_diameter_mm', 0):.1f} mm"
                f"   max. 3D-Durchmesser: {sh.get('max_diameter_mm', 0):.1f} mm\n")
        cx, cy, cz = s["centroid_mm"]
        f.write(f"    Schwerpunkt: ({cx:.1f}, {cy:.1f}, {cz:.1f}) mm\n")
        bx, by, bz = sh["bbox_size_mm"]
        f.write(f"    BBox       : {bx:.1f} x {by:.1f} x {bz:.1f} mm\n")
        f.write(f"    Sphaerizitaet: {sh['sphericity']:.4f}   "
                f"Soliditaet: {sh['solidity']:.4f}   Elongation: {sh['elongation']:.4f}"
                f"{flag}\n")

    path = output_dir / "statistics.txt"
    with open(path, "w", encoding="utf-8") as f:
        f.write("RTSTRUCT Analyse - Statistiken\n")
        f.write("=" * 60 + "\n\n")

        meta = results.get("meta", {})
        if meta:
            f.write("KATEGORIEN: " + ", ".join(
                f"{k}={v}" for k, v in meta.get("category_counts", {}).items()) + "\n")
            if meta.get("external_names"):
                f.write(f"External (ausgeschlossen): {', '.join(meta['external_names'])}\n")
            f.write(f"POINT-Marker (ausgeschlossen): {meta.get('marker_count', 0)}\n")
            f.write("Hinweis: Sphaerizitaet & Soliditaet sind auf (0,1] beschraenkt "
                    "(konsistente Voxelmaske); fuer Mehrkomponenten-Hilfsstrukturen "
                    "sind Huellen-Formmetriken nicht aussagekraeftig.\n\n")

        # OARs nach Subtyp aufteilen
        oar_serial = {n: s for n, s in results["oars"].items() if s.get("oar_subtype") == "serial"}
        oar_parallel = {n: s for n, s in results["oars"].items() if s.get("oar_subtype") != "serial"}
        for section, struct_dict in [
            ("ZIELGEBIETE (Targets)", results["targets"]),
            ("RISIKOORGANE - seriell", oar_serial),
            ("RISIKOORGANE - parallel", oar_parallel),
            ("HILFS-/PLANUNGSSTRUKTUREN", results.get("helpers", {})),
        ]:
            if not struct_dict:
                continue
            f.write(f"\n{section}\n")
            f.write("-" * 60 + "\n")
            for name, s in struct_dict.items():
                _write_struct(f, name, s)

        dist = results.get("distances", [])
        toar = [d for d in dist if d.get("pair_type") == "target-oar"]
        gtvptv = [d for d in dist if d.get("pair_type") == "gtv-ptv"]

        if toar:
            f.write("\n\nABSTANDS-STATISTIK (nur Target<->OAR Paare)\n")
            f.write("-" * 60 + "\n")
            f.write(f"Anzahl Paare : {len(toar)}\n\n")
            for label, key in [("Min. Abstand (mm)", "min_distance_mm"),
                               ("HD95 (mm)", "hd95_mm"),
                               ("Hausdorff (mm)", "hausdorff_distance_mm"),
                               ("ASSD (mm)", "assd_mm")]:
                st = _stats([d[key] for d in toar])
                f.write(f"  {label}: Mittel {st['mean']:.2f}, Std {st['std']:.2f}, "
                        f"Min {st['min']:.2f}, Max {st['max']:.2f}\n")

            f.write("\nKRITISCHSTE PAARE (aufsteigend nach Min-Abstand, Top 30)\n")
            f.write(f"{'Paar':<46}{'Min':>7}{'HD95':>7}{'Haus':>7}{'ASSD':>7}{'Zentr':>7}\n")
            f.write("-" * 81 + "\n")
            for d in sorted(toar, key=lambda e: e["min_distance_mm"])[:30]:
                pair = f"{d['structure_a']} -> {d['structure_b']}"
                f.write(f"{pair[:45]:<46}{d['min_distance_mm']:>7.2f}{d['hd95_mm']:>7.2f}"
                        f"{d['hausdorff_distance_mm']:>7.2f}{d['assd_mm']:>7.2f}"
                        f"{d['centroid_distance_mm']:>7.2f}\n")

        if gtvptv:
            f.write("\n\nGTV<->PTV MARGIN-CHECK (pro Laesion)\n")
            f.write(f"{'GTV -> PTV':<46}{'Min':>7}{'HD95':>7}{'ASSD':>7}\n")
            f.write("-" * 67 + "\n")
            for d in gtvptv:
                pair = f"{d['structure_a']} -> {d['structure_b']}"
                f.write(f"{pair[:45]:<46}{d['min_distance_mm']:>7.2f}"
                        f"{d['hd95_mm']:>7.2f}{d['assd_mm']:>7.2f}\n")

    print(f"  Gespeichert: {path}")


# ---------------------------------------------------------------------------
# Plot 6: Naehe-Matrix  PTV x kritische OARs   (NEU)
# ---------------------------------------------------------------------------

def _ptv_rows(results: dict):
    """PTVs (oder ersatzweise alle Targets), sortiert nach Läsion."""
    ptvs = [(n, s) for n, s in results.get("targets", {}).items()
            if n.upper().startswith("PTV")]
    if not ptvs:
        ptvs = list(results.get("targets", {}).items())
    return sorted(ptvs, key=lambda kv: _category_sort_key(kv[1]))


def _critical_oar_cols(results: dict):
    """Spalten der Nähe-Matrix: serielle OARs + ausgewählte kleine parallele.

    Parallele werden per case-insensitivem Teilstring ausgewählt (Auge/eye,
    Linse/lens, Hippocampus), nicht per exaktem ROI-Namen – sonst blieben die
    Spalten bei abweichender Namenskonvention anderer Datensätze leer.
    """
    oars = results.get("oars", {})
    serial = [n for n, s in oars.items() if s.get("oar_subtype") == "serial"]
    par_keywords = ("auge", "eye", "linse", "lens", "hippocamp")
    key_par = [n for n, s in oars.items()
               if s.get("oar_subtype") != "serial"
               and any(k in n.lower() for k in par_keywords)]
    return serial + key_par


def plot_proximity_matrix(results: dict, output_dir: Path) -> None:
    """
    Heatmap der Minimalabstände PTV × kritische OARs (serielle + kleine
    parallele). Ein Blick zeigt, welche Metastase welchem Risikoorgan
    gefährlich nahe kommt. Farbschwellen: rot ≤2, orange ≤5, gelb ≤10,
    hellgrün ≤20, grün >20 mm. Containment-Strukturen (Ganzhirn, Vereinigungen)
    sind hier bewusst NICHT enthalten.
    """
    from matplotlib.colors import ListedColormap, BoundaryNorm

    toar = [d for d in results.get("distances", []) if d.get("pair_type") == "target-oar"]
    ptvs = _ptv_rows(results)
    cols = _critical_oar_cols(results)
    if not toar or not ptvs or not cols:
        return

    lut = {}
    for d in toar:
        lut[(d["structure_a"], d["structure_b"])] = d["min_distance_mm"]
    M = np.full((len(ptvs), len(cols)), np.nan)
    for i, (pn, _) in enumerate(ptvs):
        for j, on in enumerate(cols):
            v = lut.get((pn, on), lut.get((on, pn)))
            if v is not None:
                M[i, j] = v

    bounds = [0, 2, 5, 10, 20, 1e9]
    cmap = ListedColormap(["#C0392B", "#E67E22", "#F1C40F", "#A9DFBF", "#27AE60"])
    norm = BoundaryNorm(bounds, cmap.N)

    fig, ax = plt.subplots(figsize=(max(6, len(cols) * 0.95 + 2),
                                    max(4, len(ptvs) * 0.5 + 2)))
    ax.imshow(np.ma.masked_invalid(M), aspect="auto", cmap=cmap, norm=norm)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(len(ptvs)))
    ax.set_yticklabels([s.get("lesion_key") or n for n, s in ptvs], fontsize=8)
    for i in range(len(ptvs)):
        for j in range(len(cols)):
            if not np.isnan(M[i, j]):
                ax.text(j, i, f"{M[i, j]:.1f}", ha="center", va="center",
                        fontsize=7, color="white" if M[i, j] <= 5 else "black")
    ax.set_title("Nähe-Matrix: PTV × kritische OARs  –  Min-Abstand (mm)\n"
                 "rot ≤2  orange ≤5  gelb ≤10  hellgrün ≤20  grün >20",
                 fontsize=11, fontweight="bold")
    plt.tight_layout()
    path = output_dir / "proximity_matrix.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Gespeichert: {path}")


# ---------------------------------------------------------------------------
# Plot 7: Naechstes kritisches OAR pro PTV   (NEU)
# ---------------------------------------------------------------------------

def plot_nearest_critical_oar(results: dict, output_dir: Path) -> None:
    """Pro PTV das EINE nächstgelegene serielle OAR (Triage-Liste).

    Balken = Minimalabstand, beschriftet mit dem betreffenden OAR. Aufsteigend
    sortiert (gefährlichste Läsion oben), Schwellen bei 3 und 5 mm.
    """
    toar = [d for d in results.get("distances", [])
            if d.get("pair_type") == "target-oar" and d.get("oar_subtype") == "serial"]
    if not toar:
        return
    ptvs = [n for n, _ in _ptv_rows(results)]

    rows = []
    for p in ptvs:
        cands = [(d["structure_b"] if d["structure_a"] == p else d["structure_a"],
                  d["min_distance_mm"])
                 for d in toar if p in (d["structure_a"], d["structure_b"])]
        if cands:
            oar, mn = min(cands, key=lambda c: c[1])
            rows.append((p, oar, mn))
    if not rows:
        return
    rows.sort(key=lambda r: r[2], reverse=True)   # größter unten, kleinster oben
    labels = [results["targets"][p].get("lesion_key") or p for p, _, _ in rows]
    vals = [r[2] for r in rows]
    colors = ["#C0392B" if v < 3 else ("#E67E22" if v < 5 else COLOR_MIN) for v in vals]

    fig, ax = plt.subplots(figsize=(9, max(3.5, len(rows) * 0.4 + 1.5)))
    y = np.arange(len(rows))
    ax.barh(y, vals, color=colors, edgecolor="white", linewidth=0.5)
    for yi, (_, oar, mn) in zip(y, rows):
        ax.text(mn + 0.5, yi, f"{mn:.1f} mm  → {oar}", va="center", ha="left", fontsize=8)
    ax.axvline(5.0, color="red", linestyle="--", linewidth=1.2, alpha=0.7, label="5 mm")
    ax.axvline(3.0, color="darkred", linestyle=":", linewidth=1.2, alpha=0.7, label="3 mm")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Min-Abstand zum nächsten seriellen OAR (mm)", fontsize=11)
    ax.set_title("Nächstes kritisches OAR pro PTV", fontsize=12, fontweight="bold")
    ax.set_xlim(0, max(vals) * 1.35 + 2)
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(axis="x", alpha=0.3, linestyle="--")
    plt.tight_layout()
    path = output_dir / "nearest_critical_oar.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Gespeichert: {path}")


# ---------------------------------------------------------------------------
# Plot 8: Sphaerizitaet vs. Elongation   (NEU)
# ---------------------------------------------------------------------------

def plot_sphericity_vs_elongation(results: dict, output_dir: Path) -> None:
    """Streudiagramm des Formcharakters (nur anatomische Strukturen).

    x = Elongation (≥1), y = Sphärizität (0..1), Größe ∝ √Volumen, Farbe =
    Kategorie. Reale Organe trennen sich natürlich (Augen ~rund; Rückenmark /
    Sehnerven stark elongiert). Hilfsstrukturen sind ausgeschlossen.
    """
    rows = [(n, s) for n, s in _iter_structures(results)
            if s.get("category") in (CAT_TARGET, CAT_OAR_SERIAL, CAT_OAR_PARALLEL)
            and s["shape"].get("shape_valid", False)]
    if not rows:
        return
    fig, ax = plt.subplots(figsize=(9, 7))
    for cat in (CAT_TARGET, CAT_OAR_SERIAL, CAT_OAR_PARALLEL):
        grp = [s for _, s in rows if s.get("category") == cat]
        if not grp:
            continue
        xs = [s["shape"]["elongation"] for s in grp]
        ys = [s["shape"]["sphericity"] for s in grp]
        sz = [max(30, min(700, np.sqrt(max(s["volume_cm3"], 0.0)) * 45)) for s in grp]
        ax.scatter(xs, ys, s=sz, color=CATEGORY_COLORS[cat], alpha=0.7,
                   edgecolors="white", linewidths=0.6, label=CATEGORY_LABELS[cat])
    # serielle OARs beschriften (klinisch wichtig + meist Ausreißer)
    for n, s in rows:
        if s.get("oar_subtype") == "serial":
            ax.annotate(n, (s["shape"]["elongation"], s["shape"]["sphericity"]),
                        fontsize=7, xytext=(4, 4), textcoords="offset points")
    ax.set_xlabel("Elongation (Hauptachsen-Verhältnis, ≥1)", fontsize=11)
    ax.set_ylabel("Sphärizität (1 = Kugel)", fontsize=11)
    ax.set_title("Formcharakter: Sphärizität vs. Elongation\n"
                 "(Markergröße ∝ √Volumen; Hilfsstrukturen ausgeschlossen)",
                 fontsize=12, fontweight="bold")
    ax.set_ylim(0, 1.05)
    ax.grid(alpha=0.3, linestyle="--")
    ax.legend(fontsize=9)
    plt.tight_layout()
    path = output_dir / "sphericity_vs_elongation.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Gespeichert: {path}")


# ---------------------------------------------------------------------------
# Plot 9: GTV->PTV Margin pro Laesion   (NEU)
# ---------------------------------------------------------------------------

def plot_gtv_ptv_margin(results: dict, output_dir: Path) -> None:
    """Pro Läsion GTV- vs. PTV-Volumen (log) plus impliziter isotroper Margin.

    Macht sichtbar, ob jede Metastase einen konsistenten GTV→PTV-Margin
    erhalten hat (Margin ~ (d_PTV − d_GTV)/2 aus Äquivalent-Durchmessern).
    """
    targets = results.get("targets", {})
    by_key: dict = {}
    for n, s in targets.items():
        key = s.get("lesion_key") or n
        prefix = "GTV" if n.upper().startswith("GTV") else ("PTV" if n.upper().startswith("PTV") else None)
        if prefix:
            by_key.setdefault(key, {})[prefix] = s
    pairs = [(k, d["GTV"], d["PTV"]) for k, d in by_key.items() if "GTV" in d and "PTV" in d]
    if not pairs:
        return
    pairs.sort(key=lambda t: t[2]["volume_cm3"])

    keys = [p[0] for p in pairs]
    gtv_v = [p[1]["volume_cm3"] for p in pairs]
    ptv_v = [p[2]["volume_cm3"] for p in pairs]
    margins = [max(0.0, (p[2]["shape"]["equivalent_diameter_mm"]
                         - p[1]["shape"]["equivalent_diameter_mm"]) / 2.0) for p in pairs]

    y = np.arange(len(keys))
    h = 0.38
    fig, ax = plt.subplots(figsize=(9, max(3.5, len(keys) * 0.5 + 1.5)))
    vmin = max(0.001, min(gtv_v) * 0.5)   # stets < jedes Volumen -> Breite > 0
    # Breite = Volumen - vmin, damit die Balkenspitze auf der Log-Skala exakt
    # auf dem echten Volumen liegt (left=vmin + width=volume würde überzeichnen).
    ax.barh(y + h / 2, np.array(gtv_v) - vmin, height=h, color="#2471A3",
            label="GTV", left=vmin)
    ax.barh(y - h / 2, np.array(ptv_v) - vmin, height=h, color="#85C1E9",
            label="PTV", left=vmin)
    for yi, g, p, m in zip(y, gtv_v, ptv_v, margins):
        ax.text(p * 1.05, yi - h / 2, f"{p:.2f} cm³", va="center", ha="left", fontsize=7)
        ax.text(g * 1.05, yi + h / 2, f"{g:.2f}", va="center", ha="left", fontsize=7)
        ax.text(vmin, yi, f"  Margin≈{m:.1f} mm", va="center", ha="left",
                fontsize=7, color="black", fontweight="bold")
    ax.set_xscale("log")
    ax.set_xlim(vmin, max(ptv_v) * 2.5)
    ax.set_yticks(y)
    ax.set_yticklabels(keys, fontsize=8)
    ax.set_xlabel("Volumen (cm³, log)", fontsize=11)
    ax.set_title("GTV→PTV pro Läsion (Volumen + impliziter Margin)",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(axis="x", alpha=0.3, linestyle="--", which="both")
    plt.tight_layout()
    path = output_dir / "gtv_ptv_margin.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
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
    # Pro-Plot-Isolation: ein fehlschlagender Plot darf die anderen nicht kippen.
    plots = (
        lambda: plot_transform_overview(orig_ds, new_ds, center, drehpunkt_pos,
                                        markers, output_dir),
        lambda: plot_displacement_summary(orig_ds, T, translation, output_dir),
        lambda: plot_transform_3d(orig_ds, new_ds, center, drehpunkt_pos, markers,
                                  output_dir, T=T, geom=geom, volume_hu=volume_hu,
                                  ct_surface=ct_surface),
    )
    for fn in plots:
        try:
            fn()
        except Exception as e:
            print(f"  (!) Case-Plot uebersprungen: {e}")
    print("Case-Visualisierung abgeschlossen.")


# ---------------------------------------------------------------------------
# Haupt-Workflow
# ---------------------------------------------------------------------------

def run_visualization(results: dict, output_dir: Path) -> None:
    """Erstellt alle Plots fuer ein bereits analysiertes results-Dict."""
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"\nErstelle Visualisierungen in {output_dir} ...")
    expected = {
        "plot_volumes": "volumes.png",
        "plot_shape_metrics": "shape_metrics.png",
        "plot_distances": "distances.png",
        "plot_centroids_3d": "centroids_3d.png",
        "plot_proximity_matrix": "proximity_matrix.png",
        "plot_nearest_critical_oar": "nearest_critical_oar.png",
        "plot_sphericity_vs_elongation": "sphericity_vs_elongation.png",
        "plot_gtv_ptv_margin": "gtv_ptv_margin.png",
    }
    for fn in (plot_volumes, plot_shape_metrics, plot_distances,
               plot_centroids_3d, plot_proximity_matrix,
               plot_nearest_critical_oar, plot_sphericity_vs_elongation,
               plot_gtv_ptv_margin):
        try:
            fn(results, output_dir)
            png = expected.get(fn.__name__)
            if png and not (output_dir / png).exists():
                # Plot hat sich wegen fehlender Daten still beendet -> sichtbar machen
                print(f"  (-) {fn.__name__} uebersprungen (keine passenden Daten)")
        except Exception as e:   # ein fehlgeschlagener Plot darf den Rest nicht stoppen
            print(f"  (!) {fn.__name__} uebersprungen: {e}")
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
