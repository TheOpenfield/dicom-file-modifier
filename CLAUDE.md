# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Python toolkit for analyzing, modifying, and visualizing DICOM RT Structure Set (RTSTRUCT) files and CT DICOM series used in radiotherapy planning. The package is `dicom_file_modifier/` and ships three runnable modules.

## Common Commands

Install dependencies:
```bash
pip install -r requirements.txt
```

Each module is invoked as `python -m dicom_file_modifier.<module>`:

```bash
# Geometric analysis of an RTSTRUCT file → JSON + console summary
python -m dicom_file_modifier.analyzer data/0000000171/RS.dcm --output output/

# Plots + statistics.txt from an RTSTRUCT file
python -m dicom_file_modifier.visualizer data/0000000171/RS.dcm --output output/
# Optional structure selection:
python -m dicom_file_modifier.visualizer data/0000000171/RS.dcm \
    --targets GTV,PTV --oars Hirnstamm,Rueckenmark --output output/

# Rigid body transform of a CT series (translation mm, rotation deg)
python -m dicom_file_modifier.modifier data/0000000171/CT \
    --tx 10 --ty 0 --tz -5 --rx 0 --ry 0 --rz 15 \
    --output output/ct_transformed
```

Modifier-specific flags worth knowing: `--method {resample,metadata}` (default `resample`; `metadata` keeps pixel data byte-identical and only rewrites IPP/IOP), `--order {0,1,3}` (interpolation order; 0 preserves exact discrete HU values), `--no-viz` to skip the Plotly HTML output.

There is no test suite, lint config, or build step in this repo.

## Architecture

Three independent modules, each importable as a library and runnable as a CLI. They share no internal state — coupling is via files on disk (`data/` for inputs, `output/` for results).

### `analyzer.py` — RTSTRUCT geometric analysis
Pipeline: `load_rtstruct` → `extract_contours` per ROI → metric functions → `run_analysis` aggregates everything into a single results dict that is also written as `<rtstruct-stem>_analysis.json`.

Key conventions:
- Contours are kept as a `list[np.ndarray]`, one (N,3) array per slice. `contours_to_points` flattens them when a unified point cloud is needed.
- Volume is planimetric (Shoelace area × mean slice spacing). Surface area for sphericity uses the convex hull of the point cloud, not the true contour surface.
- Distance computations subsample point clouds (defaults: 5000 for min-distance via `cKDTree`, 3000 for Hausdorff) for tractability — see the README "Performance Notes" section.
- Structure classification (`PTV`/`CTV`/`GTV`/`OAR`) comes from `RTROIObservationsSequence.RTROIInterpretedType`.

### `modifier.py` — CT rigid body transform
Two transform paths sharing the same affine math:
- `--method resample` (default): builds the voxel-to-patient affine `A` from IOP/IPP/PixelSpacing/slice-spacing, composes `M = A⁻¹ T⁻¹ A`, and **inverse-maps** each output voxel back into the source volume using `scipy.ndimage.map_coordinates`. Processed in **20-slice chunks** to keep peak memory ~80 MB on typical 512×512×320 volumes. Out-of-bounds voxels are filled with −1000 HU.
- `--method metadata`: pixel bytes untouched; only `ImagePositionPatient` and `ImageOrientationPatient` are rewritten via the forward transform `T`. Guarantees exact HU preservation but produces non-axial slices that some TPS may not accept.

Rotations use **extrinsic XYZ Euler angles** (`Rotation.from_euler("XYZ", ...)` — uppercase = extrinsic) about the volume's geometric centre; the offset is folded into `T` so a single 4×4 matrix represents the whole transform. All output series get fresh `SeriesInstanceUID` and per-slice `SOPInstanceUID`s. Optional Plotly HTML viz extracts surfaces with marching cubes (`skimage.measure.marching_cubes`).

### `visualizer.py` — plots and statistics
Takes either an analysis-results dict (from `analyzer.run_analysis`) or, via its CLI, runs the analyzer first and then renders. Outputs are written into `output/`: `volumes.png`, `shape_metrics.png`, `distances.png` (top-25 most-critical Target↔OAR pairs), `centroids_3d.png`, and `statistics.txt`. The CLI accepts the same `--targets`/`--oars` filters as the analyzer for selecting which ROIs to include.

## Notes for editing

- DICOM patient coordinate system is **LPS** (X=left, Y=posterior, Z=superior); all distances/translations are in mm, rotations in degrees.
- `data/` and `output/` are gitignored — don't commit DICOM files or generated artifacts.
- Some user-facing strings and argparse help text are in German; keep that consistent within each module rather than mixing languages.
- The README contains substantial mathematical documentation (volume, sphericity, Hausdorff, affine math, interpolation orders) — consult it before changing the geometric formulas, since the implementations are derived from those exact definitions.
