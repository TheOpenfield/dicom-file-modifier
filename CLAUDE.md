# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Python toolkit for analyzing, modifying, and visualizing DICOM RT Structure Set (RTSTRUCT) files and CT DICOM series used in radiotherapy planning. The package is `dicom_file_modifier/` and ships four runnable modules.

## Common Commands

Install dependencies (Python ≥ 3.12; `uv.lock` pins the tested versions on Python 3.14):
```bash
uv sync                  # .venv with Python 3.14 and the locked dependencies
pip install -e .         # alternative: minimum versions from pyproject.toml (requirements.txt does the same)
```

Each module is invoked as `python -m dicom_file_modifier.<module>`:

```bash
# Geometric analysis of an RTSTRUCT file → JSON + console summary
python -m dicom_file_modifier.analyzer data/<case-id>/RS.dcm --output output/
# Synthetic geometry self-test of the analyzer (XOR holes, z-gaps; exit 0 = pass)
python -m dicom_file_modifier.analyzer --self-test

# Plots + statistics.txt from an RTSTRUCT file
python -m dicom_file_modifier.visualizer data/<case-id>/RS.dcm --output output/
# Optional structure selection:
python -m dicom_file_modifier.visualizer data/<case-id>/RS.dcm \
    --targets GTV,PTV --oars Hirnstamm,Rueckenmark --output output/

# Rigid body transform of a CT series (translation mm, rotation deg)
python -m dicom_file_modifier.modifier data/<case-id>/CT \
    --tx 10 --ty 0 --tz -5 --rx 0 --ry 0 --rz 15 \
    --output output/ct_transformed

# Lockstep transform: CT + RTSTRUCT in one go (case-folder workflow)
python -m dicom_file_modifier.case_modifier data/<case-id> \
    --tx 10 --ty 0 --tz -5 --rx 0 --ry 0 --rz 15 \
    --center marker:HS1 --output output/run1
# List POINT-type markers (potential rotation centres) and exit
python -m dicom_file_modifier.case_modifier data/<case-id> --list-markers
# Identity-transform self-test (exit 0 = pass)
python -m dicom_file_modifier.case_modifier data/<case-id> --self-test

# Dose indices (Paddick CI, GI, ICRU 83 HI, ...) from RS + RD (+ RP) → console/JSON/TXT/CSV
# plus a separate RTSTRUCT with isodose ROIs and the intersection/underdosed/spill helper contours
python -m dicom_file_modifier.dose_indices data/<case-id> --output output/
# Every setting on Eclipse conventions (CT-pixel-aligned grid, half end slabs, global PIV, field isolines)
python -m dicom_file_modifier.dose_indices data/<case-id> --eclipse-compat high --label _ECL
# Compare with Eclipse: the RTDOSE DVHSequence is read automatically; further TPS values via
# data/<case-id>/eclipse_ref.json (--eclipse-ref) or the command line; --no-viz skips validation.html/dose_overview.png
python -m dicom_file_modifier.dose_indices data/<case-id> --eclipse-compat high --eclipse-values "PIV=1.45,V10Gy=6.1"
# ROI table + prescriptions, analytic phantom self-test, writer self-test (exit 0 = pass)
python -m dicom_file_modifier.dose_indices data/<case-id> --list
python -m dicom_file_modifier.dose_indices --self-test
python -m dicom_file_modifier.rtstruct_writer --self-test

# Synthetic demo case (CT + RS + RP + RD, no patient data; deterministic UIDs, demo_expected.json)
python -m dicom_file_modifier.demo output/demo [--layout flat] [--rd-set plan+2beams|2plans] \
    [--no-dvh] [--explicit-vr] [--rle] [--eclipse-ref ptv1|both] [--no-empty-roi] [--size medium]
# Golden-output harness (CLI parity); work dir %USERPROFILE%\dfm-golden (outside the repo)
python tools/golden.py make-inputs
python tools/golden.py run --label G0 [--only PATTERN] [--real data/<case-id> ...]
python tools/golden.py compare G0 G1 [--mode exact|migration] [--expected] [--stats]
```

`dose_indices` flags: `--target NAME[,NAME]` (default: PTVs classified as targets; a RTPLAN `DoseReferenceDescription` prefix narrows to one), `--rx GY` / `--rx-pct-of-max PCT` (default: RTPLAN `TargetPrescriptionDose`), `--isodose 100,50[,80,12Gy]`, `--grid {1.0,0.5,0.25,0.1}` (in-plane mm, default 0.25; z stays on the dose planes), `--dose-interp {linear,cubic}`, `--volume-model {slab,eclipse}` (eclipse = end slabs count half), `--piv-scope {global,component}` (default component: only the Rx-isodose component(s) overlapping the target), `--iso-contours {mask,field}`, `--eclipse-compat {high,default}` (sets grid = 1 or 2 CT pixels aligned to the CT pixel raster, eclipse volume model, global PIV, linear, field isolines, no simplification), `--label` (default `_IDX`), `--no-rs`, `--include-target`, `--simplify-mm`, `--transfer-syntax`, `--max-name-len`, `--append-csv PATH` (cross-case collection table; refuses a file with an older header), `--eclipse-ref PATH` (default `<case>/eclipse_ref.json` if present), `--eclipse-values "[NAME:]KEY=VAL,..."` (aliases TV/VRX/PIV/V50/CI/GI/HI/D98/...; `V<n>Gy` resolved via Rx), `--eclipse-tol-pct` (default 5; rows outside are marked `!` and listed as warnings, exit code unchanged), `--no-eclipse-dvh`, `--no-viz`, `--no-viz-ct`.

Modifier-specific flags worth knowing: `--method {resample,metadata}` (default `resample`; `metadata` keeps pixel data byte-identical and only rewrites IPP/IOP), `--order {0,1,3}` (interpolation order; 0 preserves exact discrete HU values), `--no-viz` to skip the Plotly HTML output.

`case_modifier` adds: `--center {volume,marker:NAME,x,y,z}` (rotation centre; default = interactive prompt with marker list, or volume centre if `--non-interactive`), `--label TEXT` (suffix for output dir / RS filename / `StructureSetLabel` / `SeriesDescription`; default `_RB`), `--new-frame-of-reference` (mints a new `FrameOfReferenceUID` for the transformed pair; default behaviour keeps the original FoR for legacy plan/dose linkage), `--dry-run`, `--verify`, `--no-viz` (skip the before/after plots), `--viz-ct-surface` (also extract the CT body surface into the 3D HTML).

There is no pytest suite or lint config yet. `pyproject.toml` declares pytest and ruff as dependency groups and builds the package with hatchling. The only automated checks are:
- `python -m dicom_file_modifier.analyzer --self-test`: synthetic geometry (XOR holes, keyhole contours, z-gaps).
- `python -m dicom_file_modifier.case_modifier <case> --self-test`: rigid-body round trip. It runs on the synthetic demo case too.
- `python -m dicom_file_modifier.dose_indices --self-test`: analytic sphere phantom with closed-form CI/GI/HI expectations.
- `python -m dicom_file_modifier.rtstruct_writer --self-test`: mask → contour → RTSTRUCT round trip.
- The golden-output harness `tools/golden.py`.

**Golden harness (`tools/golden.py`):**
- It runs 43 scenarios as subprocesses on the demo variants, plus optional local real cases (`real<n>_*`, case folder name replaced by `REAL<n>`).
- One scenario is `cli_surface`: `tools/cli_surface.py` captures every `main()` parser's options, dests, defaults, choices and types.
- It writes normalized snapshots. UIDs become `UID<n>` by first occurrence; paths, timestamps and `tempfile` names become placeholders; DICOM is stored as a tag dump plus a pixel hash, and CT pixel arrays go to `arrays.npz`.
- The DICOM dump is independent of the pydicom version:
  - `FileMetaInformationGroupLength` is masked, because pydicom 3 mints random UIDs of 62–64 characters.
  - Tags without a keyword are written as `(gggg,eeee)`.
  - Ambiguous VRs (`US or SS`) with integer values are stored as integers.
- Every step records its runtime and, on Windows, the peak working set and peak commit of the child process. Neither is compared.
- `compare --mode exact` is for refactoring PRs. `--mode migration` is for environment upgrades and uses the tolerances in `tools/golden_rules.json`:
  - Integer strings (IS, counts) stay exact.
  - Contour rings whose point count changed are compared geometrically (`contour`: same plane, Hausdorff ≤ simplify + extra; XOR area only when no ring is a vertex subset of the other).
- `--stats` appends to the report:
  - every numeric deviation, including those within tolerance, per rule;
  - pixel differences per scenario;
  - runtime and peak memory;
  - stderr warnings that are new or gone.
- Run it before and after every change to a CLI module. Scenarios tagged `known_bug` record current failures that later fixes are expected to change (`expected_changes` in the rules file).
- The work dir `%USERPROFILE%\dfm-golden` must stay outside the repo, because real-case runs contain patient data.
- No goldens are committed. The synthetic snapshots are about 50 MB, and exact parity is only reliable on the same CPU.
- The local reference is `G3b` (Python 3.14 from `uv.lock`).
- CI (Phase 1) is to run the base and the PR commit on the same runner and compare them exactly.
- `docs/migration-py314.md` documents the 3.8 → 3.14 drift. `tools/legacy/requirements-py38.txt` pins the old environment.

## Architecture

Five runnable modules plus three libraries, and the synthetic test-case generator `demo` (with `phantom` and `_compat`). `analyzer`, `modifier`, and `visualizer` are independent — they share no internal state and couple only via files on disk (`data/` inputs, `output/` results). `case_modifier` is an orchestrator: it imports building blocks from `modifier` and `analyzer` to transform a CT and its companion RTSTRUCT in lockstep. `dose_indices` is a second orchestrator for RS + RD (+ RP) built on the libraries `dose.py` (dose grid numerics, Eclipse DVH reader), `rtstruct_writer.py` (isodose RTSTRUCT export) and `dose_viz.py` (validation view).

### `demo.py` / `phantom.py` / `_compat.py` — synthetic test case
`demo.make_demo_case(out_dir, DemoSpec(...))` writes a head-first-supine phantom case.
- **Contents:**
  - CT: 128² px at 2 mm, 80 slices at 1 mm, unsigned 12-bit with intercept −1024.
  - RS: 14 ROIs incl. `Rückenmark` (ISO_IR 100), an XOR-hole OAR, a two-component OAR, `h_PTV_gesamt` typed PTV, a CONTROL ring, POINT markers `HS1`/`Iso` and an empty ROI.
  - RP: two DoseReferences, PTV_1 20 Gy and PTV_2 18 Gy.
  - RD: 1 mm grid, uint32 at `DoseGridScaling` 1e-6, relative GFOV, `DVHSequence` for PTV_1/PTV_2/BODY/Hirnstamm.
- **File names** follow Eclipse (`CT.<uid>.dcm`, `RS.<uid>.dcm`, ...).
- **UIDs** are deterministic via `generate_uid(PYDICOM_ROOT_UID, entropy_srcs=[seed, role, i])`.
- **DS values** are formatted by `demo._ds` (≤ 16 chars, identical on pydicom 2/3); files are written through `_compat.dcmwrite_file_format`.
- **`phantom.py`:**
  - `HillField` gives one smooth dose field per target.
  - `SphereTarget.expected` computes closed-form slab/eclipse values and PIV/PIV50 component vs global, neglecting cross-talk between the fields (< 0.02 Gy).
  - `cumulative_dvh` builds a numeric DVH on a fine grid. The BODY DVH is sampled only inside the dose-grid box; outside it counts as 0 Gy, like a TPS.
- **Validation:** the default case passes `validate_ct_geometry`, `validate_index_against_rs`, `dose_grid_from_dataset` (no warnings) and `read_dvh_sequence` (no notes). `dose_indices` reproduces `demo_expected.json` within ~1 %.
- **Scope:** head-first only. Feet-first is out of scope for the tool.

### `dose.py` / `dose_indices.py` / `rtstruct_writer.py` / `dose_viz.py` — dose indices
Pipeline: `discover_dose_case` (RS*/RD*/RP*.dcm, CT/ optional) → `dose.dose_grid_from_dataset` (`DoseGrid`: float32 Gy array `(k,j,i)`, affine with the same `P = A @ [k,j,i,1]` convention as `modifier.extract_geometry`, GFOV relative/absolute, validation of units/GFOV/FoR) → `select_targets` / `resolve_prescription` / `parse_isodose_levels` → `compute_dose_indices`: one `FineGrid` (in-plane `--grid`, z = native dose planes ∩ CT planes, bbox = targets ∪ lowest isodose level + margin; `--eclipse-compat` aligns the voxel centres to the CT pixel raster) → `sample_dose_on_grid` (`map_coordinates` per plane on a cropped native array; cubic prefilters once like `modifier.resample_volume`; NaN outside the grid) → `rasterize_structure` (wrapper around `analyzer.rasterize_contours`, XOR per ring, plus per-plane slab weights: `eclipse` halves the first/last slab of every contiguous z-run) → isodose masks, `ndimage.label` components (`--piv-scope component` keeps only components overlapping the target) → `evaluate_target` (volumes as weighted voxel sums, weighted DVH percentiles, all indices) → `build_eclipse_reference` / `compare_with_eclipse` (RTDOSE `DVHSequence` via `dose.read_dvh_sequence`, `eclipse_ref.json`, `--eclipse-values`; precedence cli > json > dvh, gaps filled by `derive_eclipse_values`) → RS export → `dose_viz.run_dose_visualization` (`validation.html` + `dose_overview.png`, wrapped in try/except like the RS export) → report/JSON/TXT/CSV.

Key conventions:
- All volumes come from the same fine grid; target-derived volumes (TV, TV∩PIV, underdosed) carry the slab weights, isodose-derived ones (PIV, PIV50, spill) do not, so TV∩PIV ≤ min(TV, PIV) always holds. Comparisons are inclusive (`>=`). `D_x` = dose received by x % of the weighted target volume (D98 = 2nd weighted percentile).
- Both volume models are always evaluated (only the weights differ); the report prints the chosen one and one line for the other. NaN → `null` in JSON, empty cell in CSV.
- `rtstruct_writer.write_isodose_rtstruct` builds a *fresh* Dataset (patient/study/FoR copied from the original RS, everything else new), one `CLOSED_PLANAR` contour per ring (holes are separate rings with negative area, Eclipse convention), `ContourImageSequence` referencing the CT slice at that z, explicit VR LE with a proper `FileMetaDataset`, written under `writing_validation_mode = RAISE`, then re-read by `verify_rtstruct`. Isodose contours come from `dose.mask_to_contours` (marching squares at 0.5) or, with `--iso-contours field` / `--eclipse-compat`, from `dose.field_to_contours` (isoline of the sampled dose; one vertex coordinate lies exactly on the grid lines, which is how Eclipse exports its "High"/"Default" resolution contours).
- ROI names: `ISO_100%_20.0Gy`, `<Ziel>_x_ISO100`, `<Ziel>_minus_ISO100`, `ISO100_minus_<Ziel>`; `analyzer.classify_structure` files CONTROL/DOSE_REGION types and these name patterns under HELPER so they stay out of the clinical plots.
- New-module console output is pure ASCII (`TV&PIV`, `>=`, `cm3`): the Windows console is cp1252 and `∩`/`≥` would raise `UnicodeEncodeError`. Files are UTF-8.
- Eclipse comparison: `dose.EclipseDVH` reads the DVH bins at the **left edge** (bin 0 = total volume, last non-empty bin at `DVHMaximumDose`); a target DVH gives TV/TV∩PIV/D-values, only a body (`EXTERNAL`) DVH gives PIV/PIV50. `compare_with_eclipse` always compares against `piv_global_cm3`/`piv50_global_cm3`, appends out-of-tolerance rows to `result['warnings']` and never changes the exit code. `CSV_COLUMNS` carries the `ecl_*`/`d_*_pct` columns, so `_check_csv_header` refuses `--append-csv` onto an old collection file before anything is computed. JSON stays additive: `targets[name].eclipse`, `meta.eclipse_reference`, `outputs.viz_html_path/viz_png_path`.
- `dose_viz` reads `tm.result`, `art.eclipse` and `art.eclipse_dvh` (not `art.results`, which is only filled after the hook). It needs `ct_index['z_to_path']` (added to `build_ct_slice_index`) to load only the CT slices inside the fine grid; wash/CT are strided to ≤128 px per axis, meshes to ≤96 px in-plane, contours are never thinned; isodose meshes come from the dose field (marching cubes at the level), the target mesh from its mask. With `--no-rs` the ROI specs are still built for the view, so the report shows the would-be ROI names. Plotly's slider restyles all traces, so legend toggles do not survive a slice change (documented limitation).

### `analyzer.py` — RTSTRUCT geometric analysis
Pipeline: `load_rtstruct` → `extract_contours` per ROI → metric functions → `run_analysis` aggregates everything into a single results dict that is also written as `<rtstruct-stem>_analysis.json`.

Key conventions:
- Contours are kept as a `list[np.ndarray]`, one (N,3) array per contour (a slice may hold several: islands and holes). `contours_to_points` flattens them when a unified point cloud is needed.
- Contours are grouped into slices (`_group_slices`, z within 0.05 mm) and all contours of a slice are XOR-combined (`_slice_geometries`: Shapely `symmetric_difference` after `make_valid`), so a nested contour is a hole (DICOM `CLOSED_PLANAR_XOR` / Eclipse semantics; duplicate contours cancel). Volume is planimetric: XOR area × nominal slice spacing (`_nominal_slice_spacing` = median of the smallest cluster of z-differences; larger differences are reported as `n_gaps` and are never bridged). Sphericity/solidity rasterise each structure onto its own local voxel grid via `rasterize_contours` (XOR per nearest z-plane, bbox-cropped `matplotlib.path`; reusable for any grid, e.g. a future dose grid) and take volume + surface from that one mask (surface via `skimage.measure.marching_cubes`); solidity divides that mask volume by the convex hull of the contour vertices, and the hull is also the sphericity fallback.
- Distance computations are exact and deterministic (one `cKDTree` query pair yields min/Hausdorff/HD95/ASSD); only clouds > 50,000 points are thinned to that cap with a fixed seed — see the README "Performance Notes" section.
- Structure classification into Target / serial-OAR / parallel-OAR / helper / external / marker is done by `classify_structure` (ROI name + `RTROIInterpretedType` + contour geometry; name rules override mistagged DICOM types).

### `modifier.py` — CT rigid body transform
Two transform paths sharing the same affine math:
- `--method resample` (default): builds the voxel-to-patient affine `A` from IOP/IPP/PixelSpacing/slice-spacing, composes `M = A⁻¹ T⁻¹ A`, and **inverse-maps** each output voxel back into the source volume using `scipy.ndimage.map_coordinates`. Processed in **20-slice chunks** to bound transient coordinate-array memory (two 4×n float64 arrays of ~160 MB each per chunk, vs. ~2.5 GB per whole-volume array) on typical 512×512×320 volumes. Out-of-bounds voxels are filled with −1000 HU.
- `--method metadata`: pixel bytes untouched; only `ImagePositionPatient` and `ImageOrientationPatient` are rewritten via the forward transform `T`. Guarantees exact HU preservation but produces non-axial slices that some TPS may not accept.

Rotations use **intrinsic XYZ Euler angles** (`Rotation.from_euler("XYZ", ...)` — in SciPy uppercase = intrinsic; the same matrix as an extrinsic ZYX rotation) about the volume's geometric centre; the offset is folded into `T` so a single 4×4 matrix represents the whole transform. All output series get fresh `SeriesInstanceUID` and per-slice `SOPInstanceUID`s. Optional Plotly HTML viz extracts surfaces with marching cubes (`skimage.measure.marching_cubes`).

### `case_modifier.py` — case-level lockstep transform of CT + RTSTRUCT
Orchestrator that takes a case folder of the form `data/<id>/CT/*.dcm` + `data/<id>/RS*.dcm` and applies the same rigid `T` to both. Pipeline:

1. `discover_case` auto-finds `CT/` subdir and the unique `RS*.dcm` (override with `--rs`); warns about sibling `RP*`/`RD*` that do not get transformed.
2. `validate_ct_geometry` enforces uniform `ImageOrientationPatient`, `PixelSpacing`, and slice spacing within 1% (it does not check that the orientation is axial or head-first). `validate_for_consistency` checks that the RTSTRUCT references the CT's `FrameOfReferenceUID`.
3. `find_point_markers` enumerates all ROIs whose `ContourGeometricType == "POINT"` — these become valid rotation centres alongside `volume` and explicit `x,y,z`. `parse_center_spec` handles all three forms; `interactive_center_prompt` is used when stdin is a TTY and no `--center` is given.
4. CT transform runs through `modifier`'s `build_rigid_transform`, `resample_volume` / `apply_metadata_transform`, and `save_ct_series`. `save_ct_series` returns a `{series_uid, frame_of_reference_uid_used, sop_map}` dict so the RS rewrite can map old→new SOP UIDs.
5. `transform_rtstruct` applies `T` to every `ContourData` triple (rigid → no point-cloud distortion; 3D volume preserved), rewrites every `ReferencedSOPInstanceUID` via `sop_map`, updates `RTReferencedSeriesSequence.SeriesInstanceUID` to point at the new CT, optionally mints a new `FrameOfReferenceUID` (--new-frame-of-reference), and inserts a synthetic POINT-type ROI named `Drehpunkt` at `centre + (tx,ty,tz)` so the planner sees the rotation centre at a glance.
6. Aria-visible metadata: `StructureSetLabel` truncated to DICOM SH (16 chars) with the suffix preserved; `StructureSetDescription` (VR ST; the tool truncates to 64 chars) carries the transform string via `build_transform_description`; `SeriesDescription` mirrors the CT's; `SeriesNumber += 1000` so the transformed series is distinct from the original.
7. `--self-test` runs three metadata-mode checks: an identity round trip (all `ContourData` within 1e-4 mm of the input) and pairwise-distance preservation under a 15° Z- and a 5° X-rotation (within 1e-3 mm). `--verify` after a real run computes per-ROI centroids and compares with `T @ centroid_orig` (centroids are linear under rigid transforms). `--dry-run` validates inputs and prints `T` + planned output paths without writing.

After a real (non-`--dry-run`) write, unless `--no-viz` is given, it calls `visualizer.run_case_visualization` with the in-memory original + transformed RS, `T`, the rotation centre, the `Drehpunkt` position, and the POINT-marker list, writing the three before/after plots into the case output dir. The call is wrapped in a try/except so a visualisation failure never invalidates the already-written transform.

Default FoR behaviour is **keep** (original FoR is preserved on transformed CT+RS), with a runtime warning that legacy plans/doses sharing the FoR will be auto-overlaid by Aria. Use `--new-frame-of-reference` to mint a fresh FoR when this is undesired.

### `visualizer.py` — plots and statistics
Takes either an analysis-results dict (from `analyzer.run_analysis`) or, via its CLI, runs the analyzer first and then renders. Outputs (into `output/`): `volumes.png` (log scale, grouped by category), `shape_metrics.png` (category heatmap of sphericity/solidity/elongation, single-component anatomical structures only), `distances.png` (lollipop of min-distance + HD95 for Target↔serial-OAR pairs, 3/5 mm thresholds), `centroids_3d.png` (category-coloured, numbered targets), plus the SRS plots `proximity_matrix.png`, `nearest_critical_oar.png`, `sphericity_vs_elongation.png`, `gtv_ptv_margin.png`, and `statistics.txt`. All single-RTSTRUCT plots key off the `analyzer.classify_structure` categories (Target / serial-OAR / parallel-OAR / helper / external / marker) via the shared `CATEGORY_COLORS` / `_iter_structures` / `_category_sort_key` helpers, so helper/union/optimisation structures and POINT markers / external body are kept out of the clinical plots. The CLI accepts the same `--targets`/`--oars` filters as the analyzer.

This module also hosts the **case-transform visualisation** used by `case_modifier` (it has no CLI of its own). `run_case_visualization(orig_ds, new_ds, center, drehpunkt_pos, translation, T, output_dir, markers=…, geom=…, volume_hu=…, ct_surface=…)` compares the original RTSTRUCT to the transformed one and emits: `transform_3d.html` (interactive Plotly — original vs transformed contour point clouds, with the rotation centre/`Drehpunkt`, translation vector, POINT markers, axis triad, and an opt-in CT body surface all as independently legend-toggleable traces), `transform_overview.png` (static tri-planar axial/coronal/sagittal before/after), and `displacement.png` (per-ROI centroid displacement `‖T(c)−c‖` vs the `‖t‖` reference line). It works from contour points alone — no CT pixels needed unless `ct_surface=True`, which calls `modifier._extract_surface`. `_structure_pointclouds` skips POINT-type contours so markers/`Drehpunkt` don't pollute the structure clouds.

## Notes for editing

- DICOM patient coordinate system is **LPS** (X=left, Y=posterior, Z=superior); all distances/translations are in mm, rotations in degrees.
- `data/` and `output/` are gitignored — don't commit DICOM files or generated artifacts.
- Assign a new `SOPInstanceUID` only via `modifier.set_sop_instance_uid`: it also updates the file-meta `MediaStorageSOPInstanceUID`, which must match and which pydicom's `save_as` does not sync.
- Some user-facing strings and argparse help text are in German; keep that consistent within each module rather than mixing languages.
- Windows `MAX_PATH` (260): Eclipse-style input names (`RS.<64-char uid>.dcm`) plus deep output folders can exceed it. For example, `analyzer --output` writes `<stem>_analysis.json` and fails with `FileNotFoundError` under the long scratchpad path. Keep test and golden work dirs short.
- The README contains substantial mathematical documentation (volume, sphericity, Hausdorff, affine math, interpolation orders) — consult it before changing the geometric formulas, since the implementations are derived from those exact definitions.
- README math: write display equations as ` ```math ` fenced blocks, not `$$…$$`, and keep inline `$…$` free of backslash-punctuation (`\{`, `\|`, `\,`, `\;`, `\!`, `\\`; use `\lbrace`, `\lVert … \rVert` etc.). GitHub's Markdown parser strips those backslash escapes before MathJax runs, which silently corrupts or breaks the formula.
