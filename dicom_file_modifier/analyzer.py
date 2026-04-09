#!/usr/bin/env python3
"""
Analyzer Analyse von Zielgebieten und Risikoorganen
aus DICOM RT Structure Set Dateien.

Berechnet pro Struktur:
  - Volumen (cm³)
  - Schwerpunkt (x, y, z in mm)
  - Bounding Box
  - Formmetriken: Sphärizität, Kompaktheit, Elongation

Berechnet zwischen Strukturen:
  - Minimaler Abstand (mm)
  - Hausdorff-Abstand (mm)
  - Schwerpunkt-Abstand (mm)

Verwendung:
  python -m dicom_file_modifier.analyzer rtstruct.dcm [--targets PTV,CTV,GTV] [--oars Parotis,Rueckenmark]
  python -m dicom_file_modifier.analyzer rtstruct.dcm --list    # Nur Strukturnamen auflisten

Benötigte Packages:
  pip install pydicom numpy scipy shapely matplotlib
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional

# Für Standalone-Ausführung
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pydicom
from scipy.spatial import ConvexHull
from scipy.spatial.distance import directed_hausdorff
from shapely.geometry import Polygon


# ---------------------------------------------------------------------------
# 1. DICOM RTSTRUCT einlesen
# ---------------------------------------------------------------------------

def load_rtstruct(filepath: str) -> pydicom.Dataset:
    """Lädt eine DICOM RTSTRUCT Datei und prüft die Modalität."""
    ds = pydicom.dcmread(filepath)
    if ds.Modality != "RTSTRUCT":
        raise ValueError(f"Datei ist keine RTSTRUCT (Modalität: {ds.Modality})")
    return ds


def get_structure_names(ds: pydicom.Dataset) -> Dict[int, str]:
    """Gibt ein Dictionary {ROI-Nummer: Name} zurück."""
    return {
        roi.ROINumber: roi.ROIName
        for roi in ds.StructureSetROISequence
    }


def get_structure_type(ds: pydicom.Dataset) -> Dict[int, str]:
    """Gibt ein Dictionary {ROI-Nummer: RT ROI Interpreted Type} zurück."""
    type_map = {}
    if hasattr(ds, "RTROIObservationsSequence"):
        for obs in ds.RTROIObservationsSequence:
            roi_num = obs.ReferencedROINumber
            rt_type = getattr(obs, "RTROIInterpretedType", "UNKNOWN")
            type_map[roi_num] = rt_type
    return type_map


# ---------------------------------------------------------------------------
# 2. Konturen extrahieren
# ---------------------------------------------------------------------------

def extract_contours(ds: pydicom.Dataset, roi_number: int) -> List[np.ndarray]:
    """
    Extrahiert die Konturen einer Struktur als Liste von Nx3 Arrays.
    Jedes Array enthält die (x, y, z) Koordinaten einer Kontur-Schicht.
    """
    contours = []
    for roi_contour in ds.ROIContourSequence:
        if roi_contour.ReferencedROINumber != roi_number:
            continue
        if not hasattr(roi_contour, "ContourSequence"):
            continue
        for contour in roi_contour.ContourSequence:
            pts = np.array(contour.ContourData).reshape(-1, 3)
            contours.append(pts)
    return contours


def contours_to_points(contours: List[np.ndarray]) -> np.ndarray:
    """Fasst alle Konturpunkte zu einem einzigen Nx3 Array zusammen."""
    if not contours:
        return np.empty((0, 3))
    return np.vstack(contours)


# ---------------------------------------------------------------------------
# 3. Volumenberechnung
# ---------------------------------------------------------------------------

def polygon_area(pts_2d: np.ndarray) -> float:
    """Berechnet die Fläche eines 2D-Polygons mit der Shoelace-Formel."""
    try:
        poly = Polygon(pts_2d)
        if not poly.is_valid:
            poly = poly.buffer(0)
        return abs(poly.area)
    except Exception:
        return 0.0


def compute_volume(contours: List[np.ndarray]) -> float:
    """
    Berechnet das Volumen in cm³.
    Summiert Konturflächen × Schichtabstand.
    """
    if len(contours) < 2:
        return 0.0

    # Schicht-Z-Werte und Flächen sammeln
    slice_data = {}
    for pts in contours:
        z = round(pts[0, 2], 3)
        area = polygon_area(pts[:, :2])
        slice_data.setdefault(z, 0.0)
        slice_data[z] += area  # Mehrere Konturen pro Schicht addieren

    z_values = sorted(slice_data.keys())
    if len(z_values) < 2:
        return 0.0

    # Mittleren Schichtabstand bestimmen
    dz = np.mean(np.diff(z_values))

    total_area = sum(slice_data.values())
    volume_mm3 = total_area * abs(dz)
    return volume_mm3 / 1000.0  # mm³ -> cm³


# ---------------------------------------------------------------------------
# 4. Schwerpunkt
# ---------------------------------------------------------------------------

def compute_centroid(contours: List[np.ndarray]) -> np.ndarray:
    """
    Berechnet den flächengewichteten Schwerpunkt (x, y, z) in mm.
    Jede Schicht wird mit ihrer Konturfläche gewichtet.
    """
    if not contours:
        return np.array([0.0, 0.0, 0.0])

    weighted_sum = np.zeros(3)
    total_weight = 0.0

    for pts in contours:
        area = polygon_area(pts[:, :2])
        centroid_xy = pts[:, :2].mean(axis=0)
        z = pts[0, 2]
        weighted_sum += area * np.array([centroid_xy[0], centroid_xy[1], z])
        total_weight += area

    if total_weight == 0:
        all_pts = contours_to_points(contours)
        return all_pts.mean(axis=0)

    return weighted_sum / total_weight


# ---------------------------------------------------------------------------
# 5. Formanalyse
# ---------------------------------------------------------------------------

def compute_shape_metrics(contours: List[np.ndarray], volume_cm3: float) -> Dict:
    """
    Berechnet Formmetriken:
    - Sphärizität: Verhältnis zur Kugel gleichen Volumens (1.0 = perfekte Kugel)
    - Kompaktheit: Volumen / Volumen der konvexen Hülle
    - Elongation: Verhältnis der Hauptachsenlängen (PCA)
    - Bounding Box: (x_min, x_max, y_min, y_max, z_min, z_max) in mm
    """
    all_pts = contours_to_points(contours)
    metrics = {
        "sphericity": 0.0,
        "compactness": 0.0,
        "elongation": 0.0,
        "bbox_mm": (0, 0, 0, 0, 0, 0),
        "bbox_size_mm": (0, 0, 0),
    }

    if len(all_pts) < 4 or volume_cm3 <= 0:
        return metrics

    # Bounding Box
    mins = all_pts.min(axis=0)
    maxs = all_pts.max(axis=0)
    metrics["bbox_mm"] = tuple(np.round(np.concatenate([mins, maxs]), 2))
    metrics["bbox_size_mm"] = tuple(np.round(maxs - mins, 2))

    # Sphärizität: (π^(1/3) * (6V)^(2/3)) / A
    # Vereinfachte Version mit konvexer Hülle als Oberflächenschätzung
    volume_mm3 = volume_cm3 * 1000.0
    try:
        hull = ConvexHull(all_pts)
        surface_area = hull.area
        if surface_area > 0:
            metrics["sphericity"] = round(
                (np.pi ** (1 / 3) * (6 * volume_mm3) ** (2 / 3)) / surface_area, 4
            )
        # Kompaktheit
        hull_vol = hull.volume
        if hull_vol > 0:
            metrics["compactness"] = round(volume_mm3 / hull_vol, 4)
    except Exception:
        pass

    # Elongation via PCA
    try:
        centered = all_pts - all_pts.mean(axis=0)
        cov = np.cov(centered.T)
        eigenvalues = np.sort(np.linalg.eigvalsh(cov))[::-1]
        if eigenvalues[-1] > 0:
            metrics["elongation"] = round(
                np.sqrt(eigenvalues[0] / eigenvalues[-1]), 4
            )
    except Exception:
        pass

    return metrics


# ---------------------------------------------------------------------------
# 6. Abstandsberechnungen
# ---------------------------------------------------------------------------

def min_distance(pts_a: np.ndarray, pts_b: np.ndarray,
                 sample_size: int = 5000) -> float:
    """Minimaler Abstand zwischen zwei Punktwolken in mm."""
    if len(pts_a) == 0 or len(pts_b) == 0:
        return float("inf")

    # Bei großen Punktwolken: zufällig samplen für Performance
    if len(pts_a) > sample_size:
        idx = np.random.choice(len(pts_a), sample_size, replace=False)
        pts_a = pts_a[idx]
    if len(pts_b) > sample_size:
        idx = np.random.choice(len(pts_b), sample_size, replace=False)
        pts_b = pts_b[idx]

    from scipy.spatial import cKDTree
    tree = cKDTree(pts_b)
    dists, _ = tree.query(pts_a, k=1)
    return float(np.min(dists))


def hausdorff_distance(pts_a: np.ndarray, pts_b: np.ndarray,
                       sample_size: int = 3000) -> float:
    """Hausdorff-Abstand zwischen zwei Punktwolken in mm."""
    if len(pts_a) == 0 or len(pts_b) == 0:
        return float("inf")

    if len(pts_a) > sample_size:
        idx = np.random.choice(len(pts_a), sample_size, replace=False)
        pts_a = pts_a[idx]
    if len(pts_b) > sample_size:
        idx = np.random.choice(len(pts_b), sample_size, replace=False)
        pts_b = pts_b[idx]

    d1 = directed_hausdorff(pts_a, pts_b)[0]
    d2 = directed_hausdorff(pts_b, pts_a)[0]
    return float(max(d1, d2))


def centroid_distance(c1: np.ndarray, c2: np.ndarray) -> float:
    """Euklidischer Abstand zwischen zwei Schwerpunkten in mm."""
    return float(np.linalg.norm(c1 - c2))


# ---------------------------------------------------------------------------
# 7. Gesamtanalyse
# ---------------------------------------------------------------------------

def analyze_structure(ds: pydicom.Dataset, roi_number: int,
                      roi_name: str) -> dict:
    """Vollständige Analyse einer einzelnen Struktur."""
    contours = extract_contours(ds, roi_number)
    all_pts = contours_to_points(contours)
    volume = compute_volume(contours)
    centroid = compute_centroid(contours)
    shape = compute_shape_metrics(contours, volume)

    return {
        "roi_number": roi_number,
        "name": roi_name,
        "num_contours": len(contours),
        "num_points": len(all_pts),
        "volume_cm3": round(volume, 3),
        "centroid_mm": tuple(np.round(centroid, 2)),
        "shape": shape,
        "contours": contours,       # für Abstandsberechnung
        "all_points": all_pts,
    }


def run_analysis(filepath: str,
                 target_names: Optional[List[str]] = None,
                 oar_names: Optional[List[str]] = None,
                 list_only: bool = False,
                 output_path: Optional[str] = None) -> Dict:
    """
    Hauptfunktion: Lädt RTSTRUCT, analysiert Strukturen, berechnet Abstände.

    Parameters
    ----------
    filepath : Pfad zur RTSTRUCT DICOM Datei
    target_names : Liste der Zielgebiet-Namen (z.B. ["PTV", "CTV", "GTV"])
                   Wenn None, werden alle als "TV" typisierten Strukturen verwendet.
    oar_names : Liste der Risikoorgan-Namen
                Wenn None, werden alle als "OAR" typisierten verwendet.
    list_only : Nur Strukturnamen auflisten
    output_path : Pfad zum Output-Ordner, um Ergebnisse zu speichern

    Returns
    -------
    Dictionary mit allen Analyseergebnissen
    """
    ds = load_rtstruct(filepath)
    names = get_structure_names(ds)
    types = get_structure_type(ds)

    print(f"\n{'=' * 60}")
    print(f"RTSTRUCT Analyse: {Path(filepath).name}")
    print(f"Patient: {getattr(ds, 'PatientName', 'N/A')}")
    print(f"Studie:  {getattr(ds, 'StudyDescription', 'N/A')}")
    print(f"Anzahl Strukturen: {len(names)}")
    print(f"{'=' * 60}")

    # Alle Strukturen auflisten
    print("\nVerfügbare Strukturen:")
    print(f"{'Nr':<6} {'Name':<30} {'Typ':<15}")
    print("-" * 51)
    for roi_num, roi_name in sorted(names.items()):
        rt_type = types.get(roi_num, "—")
        print(f"{roi_num:<6} {roi_name:<30} {rt_type:<15}")

    if list_only:
        return {"structures": names, "types": types}

    # Strukturen klassifizieren
    def match_names(name_list, roi_names_dict, type_filter=None):
        """Findet ROI-Nummern anhand von Namen oder Typ."""
        matched = {}
        if name_list:
            for roi_num, roi_name in roi_names_dict.items():
                for pattern in name_list:
                    if pattern.lower() in roi_name.lower():
                        matched[roi_num] = roi_name
                        break
        elif type_filter:
            for roi_num, roi_name in roi_names_dict.items():
                rt_type = types.get(roi_num, "")
                if rt_type.upper() in type_filter:
                    matched[roi_num] = roi_name
        return matched

    # Targets: PTV, CTV, GTV oder alles mit Typ "PTV"/"CTV"/"GTV"
    target_filter = ["PTV", "CTV", "GTV", "TV"]
    targets = match_names(target_names, names, target_filter)

    # OARs: explizit benannt oder Typ "OAR"/"ORGAN"/"AVOIDANCE"
    oar_filter = ["OAR", "ORGAN", "AVOIDANCE"]
    oars = match_names(oar_names, names, oar_filter)

    if not targets:
        print("\n⚠ Keine Zielgebiete gefunden. Verwende --targets um Namen anzugeben.")
    if not oars:
        print("⚠ Keine Risikoorgane gefunden. Verwende --oars um Namen anzugeben.")

    # Strukturen analysieren
    results = {"targets": {}, "oars": {}, "distances": []}

    print(f"\n{'=' * 60}")
    print("ZIELGEBIETE")
    print(f"{'=' * 60}")
    for roi_num, roi_name in targets.items():
        r = analyze_structure(ds, roi_num, roi_name)
        results["targets"][roi_name] = r
        _print_structure(r)

    print(f"\n{'=' * 60}")
    print("RISIKOORGANE")
    print(f"{'=' * 60}")
    for roi_num, roi_name in oars.items():
        r = analyze_structure(ds, roi_num, roi_name)
        results["oars"][roi_name] = r
        _print_structure(r)

    # Abstände berechnen
    all_structures = {**results["targets"], **results["oars"]}
    struct_list = list(all_structures.keys())

    if len(struct_list) >= 2:
        print(f"\n{'=' * 60}")
        print("ABSTÄNDE")
        print(f"{'=' * 60}")
        print(f"{'Struktur A':<20} {'Struktur B':<20} {'Min (mm)':<12} "
              f"{'Hausdorff (mm)':<16} {'Zentroid (mm)':<14}")
        print("-" * 82)

        for i in range(len(struct_list)):
            for j in range(i + 1, len(struct_list)):
                name_a = struct_list[i]
                name_b = struct_list[j]
                sa = all_structures[name_a]
                sb = all_structures[name_b]

                d_min = min_distance(sa["all_points"], sb["all_points"])
                d_haus = hausdorff_distance(sa["all_points"], sb["all_points"])
                d_cent = centroid_distance(
                    np.array(sa["centroid_mm"]),
                    np.array(sb["centroid_mm"])
                )

                dist_entry = {
                    "structure_a": name_a,
                    "structure_b": name_b,
                    "min_distance_mm": round(d_min, 2),
                    "hausdorff_distance_mm": round(d_haus, 2),
                    "centroid_distance_mm": round(d_cent, 2),
                }
                results["distances"].append(dist_entry)

                print(f"{name_a:<20} {name_b:<20} {d_min:<12.2f} "
                      f"{d_haus:<16.2f} {d_cent:<14.2f}")

    # Speichere Ergebnisse als JSON, wenn output_path angegeben
    if output_path:
        output_dir = Path(output_path)
        output_dir.mkdir(parents=True, exist_ok=True)
        filename = Path(filepath).stem + "_analysis.json"
        output_file = output_dir / filename

        # Entferne große Daten für JSON
        clean_results = {
            "targets": {},
            "oars": {},
            "distances": results["distances"]
        }
        for key, structures in [("targets", results["targets"]), ("oars", results["oars"])]:
            for name, struct in structures.items():
                clean_struct = {k: v for k, v in struct.items() if k not in ["contours", "all_points"]}
                clean_results[key][name] = clean_struct

        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(clean_results, f, indent=2, ensure_ascii=False)
        print(f"\nErgebnisse gespeichert in: {output_file}")

    return results


def _print_structure(r: dict):
    """Gibt die Analyseergebnisse einer Struktur formatiert aus."""
    s = r["shape"]
    print(f"\n  - {r['name']} (ROI #{r['roi_number']})")
    print(f"    Konturen: {r['num_contours']} Schichten, "
          f"{r['num_points']} Punkte")
    print(f"    Volumen:        {r['volume_cm3']:.3f} cm³")
    print(f"    Schwerpunkt:    x={r['centroid_mm'][0]:.1f}, "
          f"y={r['centroid_mm'][1]:.1f}, z={r['centroid_mm'][2]:.1f} mm")
    print(f"    Bounding Box:   {s['bbox_size_mm'][0]:.1f} × "
          f"{s['bbox_size_mm'][1]:.1f} × {s['bbox_size_mm'][2]:.1f} mm")
    print(f"    Sphärizität:    {s['sphericity']:.4f}  (1.0 = Kugel)")
    print(f"    Kompaktheit:    {s['compactness']:.4f}  "
          f"(Vol / konvexe Hülle)")
    print(f"    Elongation:     {s['elongation']:.4f}  "
          f"(Hauptachsen-Verhältnis)")


# ---------------------------------------------------------------------------
# 8. CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Analyse von DICOM RT Structure Sets",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Beispiele:
  %(prog)s rtstruct.dcm --list
  %(prog)s rtstruct.dcm --targets PTV,CTV,GTV --oars Parotis,Rueckenmark,Blase
  %(prog)s rtstruct.dcm   # Auto-Erkennung über DICOM RT ROI Type
        """,
    )
    parser.add_argument("file", help="Pfad zur RTSTRUCT DICOM Datei")
    parser.add_argument("--list", action="store_true",
                        help="Nur Strukturnamen auflisten")
    parser.add_argument("--targets", type=str, default=None,
                        help="Komma-getrennte Zielgebiet-Namen (z.B. PTV,CTV)")
    parser.add_argument("--oars", type=str, default=None,
                        help="Komma-getrennte Risikoorgan-Namen")
    parser.add_argument("--output", type=str, default="output",
                        help="Ausgabeordner für Ergebnisse (Standard: output)")

    args = parser.parse_args()

    target_list = args.targets.split(",") if args.targets else None
    oar_list = args.oars.split(",") if args.oars else None

    results = run_analysis(
        filepath=args.file,
        target_names=target_list,
        oar_names=oar_list,
        list_only=args.list,
        output_path=args.output,
    )

    print(f"\n{'=' * 60}")
    print("Analyse abgeschlossen.")
    print(f"{'=' * 60}\n")


if __name__ == "__main__":
    main()