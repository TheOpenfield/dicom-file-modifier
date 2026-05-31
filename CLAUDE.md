# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Python toolkit for analyzing, modifying, and visualizing DICOM RT Structure Set (RTSTRUCT) files and CT DICOM series used in radiotherapy planning. The package is `dicom_file_modifier/` and ships four runnable modules.

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

# Lockstep transform: CT + RTSTRUCT in one go (case-folder workflow)
python -m dicom_file_modifier.case_modifier data/0000000171 \
    --tx 10 --ty 0 --tz -5 --rx 0 --ry 0 --rz 15 \
    --center marker:HS1 --output output/run1
# List POINT-type markers (potential rotation centres) and exit
python -m dicom_file_modifier.case_modifier data/0000000171 --list-markers
# Identity-transform self-test (exit 0 = pass)
python -m dicom_file_modifier.case_modifier data/0000000171 --self-test
```

Modifier-specific flags worth knowing: `--method {resample,metadata}` (default `resample`; `metadata` keeps pixel data byte-identical and only rewrites IPP/IOP), `--order {0,1,3}` (interpolation order; 0 preserves exact discrete HU values), `--no-viz` to skip the Plotly HTML output.

`case_modifier` adds: `--center {volume,marker:NAME,x,y,z}` (rotation centre; default = interactive prompt with marker list, or volume centre if `--non-interactive`), `--label TEXT` (suffix for output dir / RS filename / `StructureSetLabel` / `SeriesDescription`; default `_RB`), `--new-frame-of-reference` (mints a new `FrameOfReferenceUID` for the transformed pair; default behaviour keeps the original FoR for legacy plan/dose linkage), `--dry-run`, `--verify`, `--no-viz` (skip the before/after plots), `--viz-ct-surface` (also extract the CT body surface into the 3D HTML).

There is no test suite, lint config, or build step in this repo.

## Architecture

Four runnable modules. `analyzer`, `modifier`, and `visualizer` are independent — they share no internal state and couple only via files on disk (`data/` inputs, `output/` results). `case_modifier` is an orchestrator: it imports building blocks from `modifier` and `analyzer` to transform a CT and its companion RTSTRUCT in lockstep.

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

Rotations use **intrinsic XYZ Euler angles** (`Rotation.from_euler("XYZ", ...)` — in SciPy uppercase = intrinsic; the same matrix as an extrinsic ZYX rotation) about the volume's geometric centre; the offset is folded into `T` so a single 4×4 matrix represents the whole transform. All output series get fresh `SeriesInstanceUID` and per-slice `SOPInstanceUID`s. Optional Plotly HTML viz extracts surfaces with marching cubes (`skimage.measure.marching_cubes`).

### `case_modifier.py` — case-level lockstep transform of CT + RTSTRUCT
Orchestrator that takes a case folder of the form `data/<id>/CT/*.dcm` + `data/<id>/RS*.dcm` and applies the same rigid `T` to both. Pipeline:

1. `discover_case` auto-finds `CT/` subdir and the unique `RS*.dcm` (override with `--rs`); warns about sibling `RP*`/`RD*` that do not get transformed.
2. `validate_ct_geometry` enforces uniform `ImageOrientationPatient`, `PixelSpacing`, and slice spacing within 1%. `validate_for_consistency` checks that the RTSTRUCT references the CT's `FrameOfReferenceUID`.
3. `find_point_markers` enumerates all ROIs whose `ContourGeometricType == "POINT"` — these become valid rotation centres alongside `volume` and explicit `x,y,z`. `parse_center_spec` handles all three forms; `interactive_center_prompt` is used when stdin is a TTY and no `--center` is given.
4. CT transform runs through `modifier`'s `build_rigid_transform`, `resample_volume` / `apply_metadata_transform`, and `save_ct_series`. `save_ct_series` was extended (Stage 1) to return a `{series_uid, frame_of_reference_uid_used, sop_map}` dict so the RS rewrite can map old→new SOP UIDs.
5. `transform_rtstruct` applies `T` to every `ContourData` triple (rigid → no point-cloud distortion; volume preserved), rewrites every `ReferencedSOPInstanceUID` via `sop_map`, updates `RTReferencedSeriesSequence.SeriesInstanceUID` to point at the new CT, optionally mints a new `FrameOfReferenceUID` (--new-frame-of-reference), and inserts a synthetic POINT-type ROI named `Drehpunkt` at `centre + (tx,ty,tz)` so the planner sees the rotation centre at a glance.
6. Aria-visible metadata: `StructureSetLabel` truncated to DICOM SH (16 chars) with the suffix preserved; `StructureSetDescription` carries the full transform string (LO, 64 chars) via `build_transform_description`; `SeriesDescription` mirrors the CT's; `SeriesNumber += 1000` so the transformed series sorts adjacent to but distinct from the original.
7. `--self-test` runs an identity transform end-to-end and asserts that all `ContourData` values match the input within 1e-4 mm. `--verify` after a real run computes per-ROI centroids and compares with `T @ centroid_orig` (centroids are linear under rigid transforms). `--dry-run` validates inputs and prints `T` + planned output paths without writing.

After a real (non-`--dry-run`) write, unless `--no-viz` is given, it calls `visualizer.run_case_visualization` with the in-memory original + transformed RS, `T`, the rotation centre, the `Drehpunkt` position, and the POINT-marker list, writing the three before/after plots into the case output dir. The call is wrapped in a try/except so a visualisation failure never invalidates the already-written transform.

Default FoR behaviour is **keep** (original FoR is preserved on transformed CT+RS), with a runtime warning that legacy plans/doses sharing the FoR will be auto-overlaid by Aria. Use `--new-frame-of-reference` to mint a fresh FoR when this is undesired.

### `visualizer.py` — plots and statistics
Takes either an analysis-results dict (from `analyzer.run_analysis`) or, via its CLI, runs the analyzer first and then renders. Outputs are written into `output/`: `volumes.png`, `shape_metrics.png`, `distances.png` (top-25 most-critical Target↔OAR pairs), `centroids_3d.png`, and `statistics.txt`. The CLI accepts the same `--targets`/`--oars` filters as the analyzer for selecting which ROIs to include.

This module also hosts the **case-transform visualisation** used by `case_modifier` (it has no CLI of its own). `run_case_visualization(orig_ds, new_ds, center, drehpunkt_pos, translation, T, output_dir, markers=…, geom=…, volume_hu=…, ct_surface=…)` compares the original RTSTRUCT to the transformed one and emits: `transform_3d.html` (interactive Plotly — original vs transformed contour point clouds, with the rotation centre/`Drehpunkt`, translation vector, POINT markers, axis triad, and an opt-in CT body surface all as independently legend-toggleable traces), `transform_overview.png` (static tri-planar axial/coronal/sagittal before/after), and `displacement.png` (per-ROI centroid displacement `‖T(c)−c‖` vs the `‖t‖` reference line). It works from contour points alone — no CT pixels needed unless `ct_surface=True`, which calls `modifier._extract_surface`. `_structure_pointclouds` skips POINT-type contours so markers/`Drehpunkt` don't pollute the structure clouds.

## Notes for editing

- DICOM patient coordinate system is **LPS** (X=left, Y=posterior, Z=superior); all distances/translations are in mm, rotations in degrees.
- `data/` and `output/` are gitignored — don't commit DICOM files or generated artifacts.
- Some user-facing strings and argparse help text are in German; keep that consistent within each module rather than mixing languages.
- The README contains substantial mathematical documentation (volume, sphericity, Hausdorff, affine math, interpolation orders) — consult it before changing the geometric formulas, since the implementations are derived from those exact definitions.
