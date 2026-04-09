# DICOM File Visualizer
# Visualisiert die Analyseergebnisse aus der JSON-Datei

import json
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

def load_analysis_results(json_path):
    """Lädt die Analyseergebnisse aus der JSON-Datei."""
    with open(json_path, 'r', encoding='utf-8') as f:
        return json.load(f)

def plot_volume_histogram(data, output_dir):
    """Erstellt Histogramme der Volumen für targets und oars."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Targets
    target_volumes = [struct['volume_cm3'] for struct in data['targets'].values()]
    ax1.hist(target_volumes, bins=10, alpha=0.7, color='blue', edgecolor='black')
    ax1.set_title('Volumen-Verteilung der Zielgebiete (Targets)')
    ax1.set_xlabel('Volumen (cm³)')
    ax1.set_ylabel('Anzahl')
    ax1.grid(True, alpha=0.3)

    # OARs
    oar_volumes = [struct['volume_cm3'] for struct in data['oars'].values()]
    ax2.hist(oar_volumes, bins=10, alpha=0.7, color='red', edgecolor='black')
    ax2.set_title('Volumen-Verteilung der Risikoorgane (OARs)')
    ax2.set_xlabel('Volumen (cm³)')
    ax2.set_ylabel('Anzahl')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / 'volume_histograms.png', dpi=300, bbox_inches='tight')
    plt.close()

def plot_distance_scatter(data, output_dir):
    """Erstellt Scatterplots der Abstände."""
    if not data['distances']:
        return

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))

    min_distances = [d['min_distance_mm'] for d in data['distances']]
    hausdorff_distances = [d['hausdorff_distance_mm'] for d in data['distances']]
    centroid_distances = [d['centroid_distance_mm'] for d in data['distances']]

    # Min Distance
    ax1.scatter(range(len(min_distances)), min_distances, alpha=0.7, color='green')
    ax1.set_title('Minimale Abstände')
    ax1.set_xlabel('Paar-Index')
    ax1.set_ylabel('Abstand (mm)')
    ax1.grid(True, alpha=0.3)

    # Hausdorff Distance
    ax2.scatter(range(len(hausdorff_distances)), hausdorff_distances, alpha=0.7, color='orange')
    ax2.set_title('Hausdorff-Abstände')
    ax2.set_xlabel('Paar-Index')
    ax2.set_ylabel('Abstand (mm)')
    ax2.grid(True, alpha=0.3)

    # Centroid Distance
    ax3.scatter(range(len(centroid_distances)), centroid_distances, alpha=0.7, color='purple')
    ax3.set_title('Schwerpunkt-Abstände')
    ax3.set_xlabel('Paar-Index')
    ax3.set_ylabel('Abstand (mm)')
    ax3.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / 'distance_scatter.png', dpi=300, bbox_inches='tight')
    plt.close()

def plot_sphericity_vs_volume(data, output_dir):
    """Erstellt Scatterplot Sphärizität vs Volumen."""
    fig, ax = plt.subplots(figsize=(8, 6))

    all_volumes = []
    all_sphericities = []

    # Targets
    for struct in data['targets'].values():
        all_volumes.append(struct['volume_cm3'])
        all_sphericities.append(struct['shape']['sphericity'])

    # OARs
    for struct in data['oars'].values():
        all_volumes.append(struct['volume_cm3'])
        all_sphericities.append(struct['shape']['sphericity'])

    ax.scatter(all_volumes, all_sphericities, alpha=0.7, color='teal')
    ax.set_title('Sphärizität vs Volumen')
    ax.set_xlabel('Volumen (cm³)')
    ax.set_ylabel('Sphärizität')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / 'sphericity_vs_volume.png', dpi=300, bbox_inches='tight')
    plt.close()

def calculate_statistics(data, output_dir):
    """Berechnet und speichert Statistiken als TXT."""
    stats = {}

    # Targets
    target_volumes = [struct['volume_cm3'] for struct in data['targets'].values()]
    stats['targets'] = {
        'count': len(target_volumes),
        'volume_mean': np.mean(target_volumes),
        'volume_std': np.std(target_volumes),
        'volume_min': np.min(target_volumes),
        'volume_max': np.max(target_volumes)
    }

    # OARs
    oar_volumes = [struct['volume_cm3'] for struct in data['oars'].values()]
    stats['oars'] = {
        'count': len(oar_volumes),
        'volume_mean': np.mean(oar_volumes),
        'volume_std': np.std(oar_volumes),
        'volume_min': np.min(oar_volumes),
        'volume_max': np.max(oar_volumes)
    }

    # Distances
    if data['distances']:
        min_distances = [d['min_distance_mm'] for d in data['distances']]
        hausdorff_distances = [d['hausdorff_distance_mm'] for d in data['distances']]
        centroid_distances = [d['centroid_distance_mm'] for d in data['distances']]

        stats['distances'] = {
            'count': len(min_distances),
            'min_distance_mean': np.mean(min_distances),
            'min_distance_std': np.std(min_distances),
            'hausdorff_mean': np.mean(hausdorff_distances),
            'hausdorff_std': np.std(hausdorff_distances),
            'centroid_mean': np.mean(centroid_distances),
            'centroid_std': np.std(centroid_distances)
        }

    # Speichere als TXT
    with open(output_dir / 'statistics.txt', 'w', encoding='utf-8') as f:
        f.write("DICOM Analyse Statistiken\n")
        f.write("=" * 40 + "\n\n")

        f.write("ZIELGEBIETE (Targets):\n")
        f.write(f"Anzahl: {stats['targets']['count']}\n")
        f.write(f"Volumen Mittelwert: {stats['targets']['volume_mean']:.3f} cm³\n")
        f.write(f"Volumen Standardabweichung: {stats['targets']['volume_std']:.3f} cm³\n")
        f.write(f"Volumen Minimum: {stats['targets']['volume_min']:.3f} cm³\n")
        f.write(f"Volumen Maximum: {stats['targets']['volume_max']:.3f} cm³\n\n")

        f.write("RISIKOORGANE (OARs):\n")
        f.write(f"Anzahl: {stats['oars']['count']}\n")
        f.write(f"Volumen Mittelwert: {stats['oars']['volume_mean']:.3f} cm³\n")
        f.write(f"Volumen Standardabweichung: {stats['oars']['volume_std']:.3f} cm³\n")
        f.write(f"Volumen Minimum: {stats['oars']['volume_min']:.3f} cm³\n")
        f.write(f"Volumen Maximum: {stats['oars']['volume_max']:.3f} cm³\n\n")

        if 'distances' in stats:
            f.write("ABSTÄNDE:\n")
            f.write(f"Anzahl Paare: {stats['distances']['count']}\n")
            f.write(f"Min. Abstand Mittelwert: {stats['distances']['min_distance_mean']:.3f} mm\n")
            f.write(f"Min. Abstand Std: {stats['distances']['min_distance_std']:.3f} mm\n")
            f.write(f"Hausdorff Mittelwert: {stats['distances']['hausdorff_mean']:.3f} mm\n")
            f.write(f"Hausdorff Std: {stats['distances']['hausdorff_std']:.3f} mm\n")
            f.write(f"Schwerpunkt Mittelwert: {stats['distances']['centroid_mean']:.3f} mm\n")
            f.write(f"Schwerpunkt Std: {stats['distances']['centroid_std']:.3f} mm\n")

def main():
    # Pfad zur JSON-Datei (angenommen im output/ Ordner)
    output_dir = Path('output')
    json_file = output_dir / 'RS.1.3.46.670589.13.8605032.20260409113651.374505_analysis.json'

    if not json_file.exists():
        print(f"JSON-Datei nicht gefunden: {json_file}")
        return

    # Daten laden
    data = load_analysis_results(json_file)

    # Plots erstellen
    plot_volume_histogram(data, output_dir)
    plot_distance_scatter(data, output_dir)
    plot_sphericity_vs_volume(data, output_dir)

    # Statistiken berechnen
    calculate_statistics(data, output_dir)

    print("Visualisierung abgeschlossen. Plots und Statistiken in output/ gespeichert.")

if __name__ == "__main__":
    main()