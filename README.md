# DICOM File Modifier

A comprehensive toolkit for analyzing, modifying, and visualizing DICOM RT Structure Set files used in radiotherapy planning. This project provides tools to process target volumes (PTV, CTV, GTV) and organs at risk (OAR) from DICOM files, compute geometric metrics, and generate visualizations.

> **Disclaimer — research / educational tool only.** This software is **not a medical device**. It has not been validated or certified for clinical use and must **not** be used to create, modify, or verify data for actual patient treatment. It writes and rewrites clinical DICOM (CT pixel data, RTSTRUCT contours, UID references); all output must be independently verified by a qualified medical physicist before any clinical application. Use at your own risk.

## Contents

- [RTSTRUCT Analyzer](#rtstruct-analyzer-documentation) — geometric analysis of structure sets
- [CT Rigid Body Transformer](#ct-rigid-body-transformer-documentation) — rigid transform of a CT series
- [Case Modifier](#case-modifier-documentation) — lockstep transform of a CT + its RTSTRUCT
- [RTSTRUCT Visualizer](#rtstruct-visualizer-documentation) — plots and statistics
- [Dose Indices](#dose-index-documentation) — Paddick CI, GI, ICRU 83 HI and isodose contours from RTSTRUCT + RTDOSE
- [References](#references)

## Features

- **Analyzer**: Extract and compute geometric properties (volume, centroid, shape metrics, distances) from RTSTRUCT files
- **Modifier**: Rigid body transformation of CT DICOM series (translation + rotation), either resampled onto the original axial grid or HU-exact by rewriting only the geometry tags, with 3D visualisation
- **Case Modifier**: Lockstep rigid body transformation of a CT series **and** its companion RTSTRUCT in a single pass — contour points are transformed alongside the pixel data, UID references are rewritten so the new RS links to the new CT, and an optional `Drehpunkt` POINT marker is inserted at the rotation centre for easy identification in the TPS
- **Visualizer**: Generate plots and statistics from analysis results
- **Dose Indices**: Compute Paddick conformity index, gradient index, ICRU 83 homogeneity index and DVH statistics for each target from an RTSTRUCT + RTDOSE pair (prescription from the RTPLAN), and write a separate RTSTRUCT with the isodose ROIs and the intersection / underdosed / spill helper contours for verification in the TPS; an Eclipse-compatibility mode reproduces the TPS raster conventions

## Project Structure

```
dicom-file-modifier/
├── data/                    # Input DICOM files (not synced)
│   └── <case-id>/           # One folder per case
│       ├── CT/              # CT slices (*.dcm)
│       ├── RS*.dcm          # RTSTRUCT file
│       ├── RD*.dcm          # RTDOSE file (dose_indices)
│       └── RP*.dcm          # RTPLAN file (dose_indices: prescription; optional)
├── output/                  # Analysis results and modified files (not synced)
├── dicom_file_modifier/     # Python package
│   ├── __init__.py
│   ├── analyzer.py          # RTSTRUCT analysis module
│   ├── modifier.py          # CT rigid body transformer
│   ├── case_modifier.py     # CT + RTSTRUCT lockstep transformer
│   ├── visualizer.py        # Visualisation module
│   ├── dose.py              # Dose-grid numerics library (RTDOSE, fine grid, sampling, DVH statistics)
│   ├── dose_indices.py      # Dose index computation (CI/GI/HI) + reports
│   └── rtstruct_writer.py   # Isodose / helper-contour RTSTRUCT export
├── requirements.txt         # Python dependencies
└── README.md                # This file
```

**Case folder convention:** `case_modifier` expects a folder of the form `data/<case-id>/` containing a `CT/` subfolder with the CT DICOM slices plus exactly one `RS*.dcm` file at the case root. Sibling `RP*.dcm` (RTPLAN) and `RD*.dcm` (RTDOSE) files are detected and reported but **not** transformed. `dose_indices` uses the same folder layout: it reads `RS*.dcm` and `RD*.dcm` (plus `RP*.dcm` for the prescription) and needs `CT/` only for the slice references of the exported RTSTRUCT and for the Eclipse-compatibility raster.

## Installation

1. Clone the repository
2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

### Analyzer

Analyse RTSTRUCT files and compute geometric metrics:

```bash
python -m dicom_file_modifier.analyzer data/<case-id>/RS.dcm --output output/
# Synthetic consistency tests of the analyzer geometry (holes/XOR, z-gaps); exit 0 = pass
python -m dicom_file_modifier.analyzer --self-test
```

This will generate:
- `RS_analysis.json` (`<rtstruct-stem>_analysis.json`): Computed metrics for all structures
- Console output with summary statistics

### Visualizer

Generate plots and statistics directly from a RTSTRUCT file:

```bash
python -m dicom_file_modifier.visualizer data/<case-id>/RS.dcm --output output/
# with explicit structure selection:
python -m dicom_file_modifier.visualizer data/<case-id>/RS.dcm \
    --targets GTV,PTV --oars Hirnstamm,Rueckenmark --output output/
```

This creates:
- `volumes.png`: Horizontal bar chart (log scale) of all analysed structure volumes, grouped and colour-coded by category
- `shape_metrics.png`: Category-grouped heatmap of sphericity, solidity, elongation (anatomical structures)
- `distances.png`: Lollipop of min-distance + HD95 for Target↔serial-OAR pairs, with 3/5 mm thresholds
- `proximity_matrix.png`, `nearest_critical_oar.png`, `sphericity_vs_elongation.png`, `gtv_ptv_margin.png`: additional SRS / multi-metastasis plots
- `centroids_3d.png`: 3D scatter of structure centroids in patient space, marker size ∝ √volume
- `statistics.txt`: Full numerical summary of all metrics

### Modifier – CT Rigid Body Transformer

Apply a rigid body transformation (translation + rotation) to a CT DICOM series:

```bash
python -m dicom_file_modifier.modifier data/<case-id>/CT \
    --tx 10 --ty 0 --tz -5 \
    --rx 0  --ry 0 --rz 15 \
    --output output/ct_transformed
```

This produces:
- `CT_0000.dcm … CT_NNNN.dcm`: Transformed CT series on the original axial grid (with `--method metadata` the pixels are unchanged and the slices become oblique instead)
- `visualization_3d.html`: Interactive 3D comparison (original vs. transformed)

**All CLI options:**

| Option | Default | Description |
|---|---|---|
| `--tx/ty/tz` | 0 mm | Translation along X / Y / Z axis |
| `--rx/ry/rz` | 0 ° | Rotation around X / Y / Z axis (intrinsic XYZ Euler) |
| `--method` | `resample` | `resample`: standard axial output; `metadata`: exact HU preservation |
| `--order` | `1` | Interpolation order: 0 = nearest neighbour, 1 = linear, 3 = cubic |
| `--output` | `output/ct_transformed` | Output directory |
| `--save-viz` | *(auto)* | Path for HTML visualisation file |
| `--no-viz` | off | Skip visualisation |

### Case Modifier – CT + RTSTRUCT Lockstep Transformer

Apply the **same** rigid body transformation to a CT series and its companion RTSTRUCT in one go, so the contours stay anchored to the moved anatomy:

```bash
# List POINT-type markers in the RTSTRUCT (potential rotation centres)
python -m dicom_file_modifier.case_modifier data/<case-id> --list-markers

# Identity-transform self-test (validates the round-trip; exit 0 = pass)
python -m dicom_file_modifier.case_modifier data/<case-id> --self-test

# Real transformation: 10 mm shift + 15° rotation around marker "HS1"
python -m dicom_file_modifier.case_modifier data/<case-id> \
    --tx 10 --ty 0 --tz -5 --rx 0 --ry 0 --rz 15 \
    --center marker:HS1 --output output/run1 --verify
```

This produces `output/run1/<case-id>_RB/`:
- `CT/CT_0000.dcm … CT_NNNN.dcm` — transformed CT series (new `SeriesInstanceUID` and per-slice `SOPInstanceUID`s)
- `RS_RB.dcm` — transformed RTSTRUCT (every `ContourData` triple multiplied by `T`, all SOP/Series UID references rewritten to point at the new CT, plus a synthetic `Drehpunkt` POINT-ROI inserted at the rotation centre)
- `transform_3d.html` — interactive 3D before/after comparison of the contour geometry (see [Case Transform Visualisation](#case-transform-visualisation))
- `transform_overview.png` — static tri-planar (axial / coronal / sagittal) before/after projection
- `displacement.png` — per-ROI centroid displacement bar chart

Visualisation is on by default; pass `--no-viz` to skip it (e.g., for batch/CI runs), or `--viz-ct-surface` to additionally extract the CT body surface into the 3D view.

**All CLI options:**

| Option | Default | Description |
|---|---|---|
| `--tx/ty/tz` | 0 mm | Translation along X / Y / Z axis (LPS) |
| `--rx/ry/rz` | 0 ° | Rotation around X / Y / Z axis (intrinsic XYZ Euler) |
| `--method` | `resample` | Same as modifier: `resample` (default) or `metadata` |
| `--order` | `1` | Interpolation order (resample only): 0/1/3 |
| `--center` | *(prompt)* | Rotation centre. `volume`, `marker:NAME`, or `x,y,z` (LPS, mm) |
| `--list-markers` | off | Print all POINT-type ROIs and exit |
| `--non-interactive` | off | Skip the interactive centre prompt; default to volume centre |
| `--label` | `_RB` | Suffix for output folder, RS filename, `StructureSetLabel`, `SeriesDescription` |
| `--rs` | *(auto)* | Explicit path to RTSTRUCT file (when multiple `RS*.dcm` exist) |
| `--new-frame-of-reference` | off | Mint a new `FrameOfReferenceUID` for the transformed pair |
| `--dry-run` | off | Validate inputs and print plan; write nothing |
| `--verify` | off | After write: re-read RS and check centroid linearity per ROI |
| `--self-test` | off | End-to-end identity-transform test; exit 0 = pass |
| `--no-viz` | off | Skip the before/after visualisation plots |
| `--viz-ct-surface` | off | Also extract the CT body surface (marching cubes) into `transform_3d.html` |
| `--output` | `output` | Base output directory |

**Aria/TPS-visible metadata.** The transformed series carries identifying information in the human-readable DICOM tags so a planner can spot it at a glance:
- `SeriesDescription = "<orig>_RB"` (LO, ≤ 64 chars)
- `StructureSetLabel = "<orig>_RB"` (SH, truncated to 16 chars with the suffix preserved)
- `StructureSetName = "<orig>_RB"`
- `StructureSetDescription = "rigid t=(10,0,-5) r=(0,0,15) c=Marker 'HS1' m=res FoR=keep"` (VR ST; the tool truncates it to 64 chars) — translation, rotation, centre label, method and FoR strategy
- `SeriesNumber += 1000` so the transformed series is easy to tell apart from the original (e.g. 2 → 1002)
- `--label` overrides the suffix (e.g., `--label _SHIFT_LR10`)

### Dose Indices – Paddick CI, GI, HI + Isodose RTSTRUCT

Compute plan-quality indices for the target volume(s) of a case from its RTSTRUCT + RTDOSE (prescription from the RTPLAN) and write a separate RTSTRUCT with the isodose and helper contours:

```bash
# Default: 0.25 mm evaluation grid, slab volume model, PIV restricted to the isodose component of the target
python -m dicom_file_modifier.dose_indices data/<case-id> --output output/

# Everything on Eclipse conventions (CT-pixel-aligned grid, half end slabs, global PIV, field isolines)
python -m dicom_file_modifier.dose_indices data/<case-id> --eclipse-compat high --label _ECL

# Explicit files, prescription and isodose levels
python -m dicom_file_modifier.dose_indices --rs RS.dcm --rd RD.dcm --rx 20 --isodose 100,80,50,12Gy --no-rs

# ROI table + prescriptions / analytic phantom self-test / writer self-test (exit 0 = pass)
python -m dicom_file_modifier.dose_indices data/<case-id> --list
python -m dicom_file_modifier.dose_indices --self-test
python -m dicom_file_modifier.rtstruct_writer --self-test
```

This produces `output/<case-id>_IDX/`:
- `RS_<case-id>_IDX.dcm` — separate RTSTRUCT with `ISO_100%_20.0Gy`, `ISO_50%_10.0Gy` and, per target, `<Target>_x_ISO100` (intersection), `<Target>_minus_ISO100` (underdosed) and `ISO100_minus_<Target>` (spill); `--include-target` adds a verbatim copy of the target
- `<case-id>_indices.json` — all components, indices, DVH statistics and settings
- `indices.txt` — the console report; `indices.csv` — one row per target (`--append-csv PATH` appends the rows to a cross-case collection table)

**All CLI options:**

| Option | Default | Description |
|---|---|---|
| `case_dir` / `--rs --rd --rp` | *(auto)* | Case folder with `RS*.dcm`, `RD*.dcm`, optional `RP*.dcm` and `CT/`; or explicit files |
| `--target` | *(auto)* | `NAME[,NAME…]`, exact or unique substring. Default: ROIs classified as targets whose name starts with `PTV`; a RTPLAN `DoseReferenceDescription` prefix narrows the choice to one |
| `--rx` / `--rx-pct-of-max` | *(RTPLAN)* | Prescription in Gy, or as a percentage of the grid maximum (SRS convention) |
| `--isodose` | `100,50` | Isodose levels in % of Rx or absolute (`12Gy`); 100 and 50 are always included (CI, GI) |
| `--grid` | `0.25` | In-plane resolution of the evaluation grid in mm (`1.0`, `0.5`, `0.25`, `0.1`); z stays on the dose planes |
| `--dose-interp` | `linear` | Dose interpolation: trilinear or cubic B-spline |
| `--volume-model` | `slab` | `slab`: every contour slice is a full slab; `eclipse`: the first and last slab of a structure count half (TPS convention) |
| `--piv-scope` | `component` | PIV = only the Rx-isodose component(s) overlapping the target (multi-target plans), or `global` |
| `--iso-contours` | `mask` | Isodose contours from the voxel mask (edge midpoints) or as isolines of the sampled dose field (`field`) |
| `--eclipse-compat` | off | `high` / `default`: evaluation grid on the CT pixel raster (1 or 2 pixels, see below), `eclipse` volume model, global PIV, linear interpolation, field isolines without simplification |
| `--label` | `_IDX` | Suffix for the output folder, RS file name and `StructureSetLabel` |
| `--no-rs` / `--include-target` | off | Skip the RTSTRUCT export / add a copy of the target ROI(s) |
| `--simplify-mm` | `0.1` | Douglas–Peucker tolerance of the exported contours (keep below half the grid resolution) |
| `--transfer-syntax` | `explicit` | Transfer syntax of the exported RTSTRUCT (`explicit` or `implicit` VR little endian) |
| `--max-name-len` | `64` | Maximum ROI name length (DICOM LO limit) |
| `--append-csv` | off | Append the result rows to a collection CSV (header written once) |
| `--list` / `--self-test` | off | ROI table + prescriptions; analytic phantom self-test |
| `--output` | `output` | Base output directory |

See [Dose Index Documentation](#dose-index-documentation) for the definitions, the grid conventions and the Eclipse comparison.

## Dependencies

- pydicom: DICOM file handling
- numpy: Numerical computations
- scipy: Scientific computing (distances, convex hull, interpolation, rotations)
- shapely: 2D geometry operations
- matplotlib: Plotting and visualisation
- scikit-image: Marching cubes surface extraction (analyzer sphericity surface; modifier & case_modifier CT surfaces)
- plotly: Interactive 3D visualisation (modifier, case_modifier)

---

# RTSTRUCT Analyzer Documentation

## Overview

The **RTSTRUCT Analyzer** is a Python tool for geometric analysis of target volumes (*Target Volumes*) and organs at risk (*Organs at Risk, OAR*) from DICOM RT Structure Set files. These files are created during radiotherapy planning, where physicians draw 3D contours on CT images to define the tissue to be irradiated and surrounding structures to be protected.

The script extracts these contour data and computes clinically relevant geometric parameters: volume, centroid, shape metrics, and distances between structures.

## Data Foundation: DICOM RT Structure Set

### Structure of an RTSTRUCT File

An RTSTRUCT file stores contours as ordered sequences of 3D points. Each structure (e.g., *PTV*, *CTV*, *GTV*, or an organ at risk) consists of one or more **contours per CT slice**. Each contour is a closed polygon in the respective axial plane, defined by a sequence of *(x, y, z)* coordinates in millimeters in the DICOM patient coordinate system.

Relevant DICOM fields:

| DICOM Field | Content |
|---|---|
| `StructureSetROISequence` | List of all structures with ROI number and name |
| `RTROIObservationsSequence` | Clinical classification (e.g., PTV, OAR) |
| `ROIContourSequence` | Contour data per structure and slice |
| `ContourData` | Flat list of x₁, y₁, z₁, x₂, y₂, z₂, … coordinates |

### Structure Classification

ICRU Reports 50, 62, and 83 define a hierarchy of target volumes:

- **GTV** (*Gross Tumor Volume*): Macroscopically visible tumour tissue.
- **CTV** (*Clinical Target Volume*): GTV plus the surrounding tissue that may contain subclinical (microscopic) disease.
- **PTV** (*Planning Target Volume*): CTV plus a geometric margin for setup uncertainties and organ motion.

Organs at risk (OAR) are structures whose radiation dose must be limited, e.g., spinal cord, parotid, or bladder. Classification is stored in the `RTROIInterpretedType` field.

## Calculation Methods

### Volume Calculation

The volume of a structure is calculated as a stack of slabs, one per contour slice:

```math
V = \Delta z_{\mathrm{nom}} \cdot \sum_{k=1}^{S} A_k
```

Where $S$ is the number of slices (contour z-positions closer than 0.05 mm are merged into one slice), $A_k$ is the area of slice $k$ (see below) and $\Delta z_{\mathrm{nom}}$ is the nominal slice spacing, i.e. each slice represents a slab one slice thick. The area of a single contour is determined using the **Shoelace formula** (Gauss's area formula, also known as the surveyor's formula), with the vertex index taken cyclically ($x_n = x_0$, $y_n = y_0$):

```math
A = \frac{1}{2} \left| \sum_{j=0}^{n-1} (x_j \, y_{j+1} - x_{j+1} \, y_j) \right|
```

**Several contours on one slice (holes and islands).** DICOM encodes a ring-shaped or hollow structure as an outer contour plus a nested inner contour on the same slice, and a fragmented structure as several disjoint contours. All contours of a slice are combined by **XOR** (exclusive disjunction, the even-odd rule; Shapely `symmetric_difference`). This is the semantics of the DICOM contour type `CLOSED_PLANAR_XOR` and the convention Eclipse uses for its plain `CLOSED_PLANAR` contours: a nested contour is a hole, an island inside a hole counts again, disjoint islands add up, and the order of the contours does not matter. $A_k$ is the area of the resulting region. A hole encoded with the keyhole technique (a single contour with a narrow channel) yields the same area. Two consequences of the XOR rule, both consistent with the TPS: identical duplicate contours cancel to zero area, and partially overlapping contours lose their overlap. Invalid contours (self-intersecting spikes) are repaired with Shapely `make_valid` before the combination.

**Nominal slice spacing.** $\Delta z_{\mathrm{nom}}$ is the median of the differences between consecutive slice positions, taken over those differences that do not exceed 1.5 times the smallest difference. Larger differences are **z-gaps** (a skipped slice, or separate lesions stored in one ROI). They are counted and reported as `n_gaps`, but **not bridged**: no volume is attributed to a z-range without a contour, and the gap does not inflate the spacing of the other slices. (The previous implementation used the mean spacing, which a single 30 mm gap in a 1 mm series inflates by a factor of about 1.7.) A ROI with genuinely variable slice thickness is not supported: its thicker part is counted with $\Delta z_{\mathrm{nom}}$ and flagged as gaps.

**Unit:** Results are converted from mm³ to cm³ (division by 1000).

**Limitation:** This method is a planimetric approximation. It is more accurate with thinner CT slices. A voxel-based volume on the CT's own grid (as a TPS computes it) would require the associated CT image; the analyzer's voxel cross-check uses a local grid of its own (see *Sphericity*). A structure with a single contour slice has no measurable spacing and reports 0 cm³.

### Centroid Calculation

The centroid is calculated as the **area-weighted average** of the slice centroids:

```math
\vec{C} = \frac{\sum_{k=1}^{S} A_k \cdot \vec{c}_k}{\sum_{k=1}^{S} A_k}
```

Where $\vec{c}_k$ is the area centroid (first area moment) of the XOR region of slice $k$, computed with Shapely and placed at the slice's z-position, and $A_k$ is that region's area as weight (see *Volume Calculation*). A hole therefore shifts the centroid away from itself, as it should, instead of attracting it as a positively weighted inner contour would. The area centroid is used rather than the arithmetic mean of the vertices, which would be biased towards densely sampled segments of the contour.

Weighting by area ensures that larger cross-sections contribute more to the overall centroid than small edge slices. The result is three coordinates *(x, y, z)* in millimeters in the DICOM patient coordinate system.

### Shape Metrics

Shape metrics quantify the geometric shape of a structure independently of its absolute size.

#### Sphericity

Sphericity describes how closely a structure resembles a sphere (Wadell's definition [Wadell 1935]). It is defined as the ratio of the surface area of a volume-equivalent sphere to the actual surface area:

```math
\Psi = \frac{\pi^{1/3} \cdot (6V)^{2/3}}{A_{\text{surface}}}
```

A value of 1.0 corresponds to a perfect sphere; smaller values indicate irregular or elongated shapes. By the isoperimetric inequality $\Psi \le 1$ for **any** solid, but this holds only when $V$ and $A$ describe the *same* body. The volume $V$ and surface area $A$ are therefore both taken from a **single consistent voxel mask**: the stacked contours are rasterised onto a local grid (in-plane via `matplotlib.path`, the contours of one plane combined by XOR so that holes stay open, one grid plane on every contour plane at the nominal slice spacing), the volume is the filled voxel volume, and the surface area is obtained from a marching-cubes mesh of that mask (`skimage.measure`). The result is clamped to $(0,1]$. If rasterisation is unavailable the code falls back to a convex-hull-consistent estimate.

Because the surface is measured on a discretised, slice-stacked mask (terraces between slices, marching-cubes facets), the surface area is systematically over-estimated and $\Psi$ is **biased low**: a perfect sphere contoured on 1–3 mm slices scores only about 0.79–0.92, lower for small structures and thick slices. Sphericity values are therefore meaningful for comparing structures of similar size and slice spacing, not as absolute values.

#### Solidity

Solidity (the scikit-image / ImageJ name for this metric) relates the actual volume to the volume of the convex hull:

```math
S = \frac{V_{\text{mask}}}{V_{\text{convex}}}
```

A value of 1.0 means the structure is convex; lower values indicate concave or highly irregular shapes, which may be clinically relevant for tumours that wrap around other structures. For the true body $S \le 1$, because a body is contained in its convex hull. The implementation, however, divides the voxel-mask volume — the same slab stack as the planimetric volume, which extends half a slice beyond the outermost contours — by the convex hull of the contour *vertices*, which ends at the outermost contour planes. The raw ratio is therefore biased upwards (exactly $N/(N-1)$ for a prism spanning $N$ slices; 1.01–1.13 for spheres in synthetic tests) and is clamped to $(0,1]$. As a result convex structures read 1.0, and moderate concavities of small, few-slice structures can be masked. Holes do reduce the mask volume: a ring with radii 20/10 mm over 20 slices reads about 0.79 (= 0.75 × 20/19), where an implementation that fills holes would read 1.0. For multi-component unions / dose shells / optimisation structures (flagged `n_components > 1` / `shape_valid = false`) the convex-hull metrics are not meaningful and are excluded from the shape plots.

#### Elongation

Elongation describes the stretching of a structure along its principal axes. It is determined via **principal component analysis** (PCA) of the 3D point cloud:

```math
E = \sqrt{\frac{\lambda_{\max}}{\lambda_{\min}}}
```

Where $\lambda_{\max}$ and $\lambda_{\min}$ are the largest and smallest eigenvalues of the covariance matrix of the contour vertices, so $E$ is the ratio of the standard deviations along the longest and shortest principal axes. A value of 1.0 corresponds to isotropic extension, larger values show increasing stretching in a preferred direction.

The vertices sample the *surface*, and not uniformly, so $E$ depends on how the structure was contoured. With a constant vertex spacing along each contour (vertex count ∝ perimeter) the small polar slices are under-represented and even a perfect sphere yields $E \approx 1.2$; with a fixed vertex count per contour it yields $E \approx 1.0$. Compare elongation values only between structures contoured the same way.

#### Bounding Box

The axis-aligned bounding box gives the extent of the structure in all three spatial directions:

```math
\Delta x = x_{\max} - x_{\min}, \quad \Delta y = y_{\max} - y_{\min}, \quad \Delta z = z_{\max} - z_{\min}
```

It provides a quick overview of the spatial extent in millimeters.

## Distance Calculations

Distances between structures are central to evaluating the radiotherapy plan: How close is an organ at risk to the target? All distances below are computed between the contour **vertices** of the two structures, i.e. between point samples of their surfaces. They are never negative and do not by themselves detect overlap (see *Minimum Distance*).

### Minimum Distance

The minimum distance between two structures *A* and *B* is defined as:

```math
d_{\min}(A, B) = \min_{a \in A, \, b \in B} \| a - b \|_2
```

Calculation is performed efficiently using a **KD-tree** (`scipy.spatial.cKDTree`): For each point of structure *A*, the nearest neighbour in *B* is found and the global minimum determined.

Clinically, this value is particularly relevant for assessing whether a safety margin between PTV and adjacent organs at risk is maintained. Two properties matter when reading it:

- **It is a surface-to-surface distance and cannot indicate overlap.** Structures whose surfaces intersect give a value close to, but generally not exactly, 0 mm. A structure lying entirely *inside* another gives a clearly positive value although the two overlap completely: an OAR sphere of radius 5 mm centred in a PTV sphere of radius 20 mm reports 15 mm. Overlap has to be checked separately (e.g. via an intersection volume).
- **Only vertices are compared.** The result can over-estimate the distance between the continuous surfaces by up to about the in-plane vertex spacing, or about one slice spacing where the closest approach falls between contour planes. That is the unsafe direction for OAR clearance.

### Hausdorff Distance

The Hausdorff distance [Huttenlocher et al. 1993] is a measure of the maximum deviation between two point sets:

```math
d_H(A, B) = \max\!\Big(\,\sup_{a \in A} \inf_{b \in B} \|a - b\|, \;\sup_{b \in B} \inf_{a \in A} \|a - b\|\,\Big)
```

Intuitively: The Hausdorff distance indicates how far one must travel in the worst case from a point of one structure to the nearest point of the other structure. Like the minimum distance it is set by a single point pair, so it is highly sensitive to outliers. It is the standard worst-case measure for comparing two delineations of the *same* structure; between two different structures (e.g. Target ↔ OAR) it mainly reflects their size and separation rather than a clinically meaningful margin.

### 95th-Percentile Hausdorff Distance (HD95)

Because the raw Hausdorff distance is fixed by a single worst-case point, it is dominated by isolated outliers. The 95th-percentile Hausdorff distance replaces the outer suprema with the 95th percentile of the directed nearest-neighbour distances, then takes the larger of the two directions:

```math
d_{H95}(A, B) = \max\!\Big(\,P_{95}\big\{\inf_{b \in B}\|a - b\| : a \in A\big\}, \; P_{95}\big\{\inf_{a \in A}\|b - a\| : b \in B\big\}\,\Big)
```

This is the robust variant recommended for segmentation comparison [Taha & Hanbury 2015] and is the companion metric reported in `distances.png`. Other implementations take the 95th percentile of the *pooled* distances of both directions instead, which gives slightly different values, so HD95 figures are only comparable between tools that use the same definition.

### Average Symmetric Surface Distance (ASSD)

ASSD is the mean nearest-neighbour distance, averaged symmetrically over both directions:

```math
d_{\text{ASSD}}(A, B) = \frac{\sum_{a \in A} \inf_{b \in B}\|a - b\| + \sum_{b \in B} \inf_{a \in A}\|b - a\|}{|A| + |B|}
```

This is the definition used by [Heimann et al. 2009]. Unlike the (95th-percentile) Hausdorff distance, ASSD reflects the *typical* separation over the whole boundary rather than its extremes [Taha & Hanbury 2015].

### Centroid Distance

The Euclidean distance between the centroids of two structures:

```math
d_C = \| \vec{C}_A - \vec{C}_B \|_2
```

This value provides a rough but robust estimate of spatial separation. It is insensitive to outliers and suitable as a quick comparison value.

## Performance Notes

All distance metrics are computed **exactly on the full point clouds**: a single pair of KD-tree queries (`scipy.spatial.cKDTree`) yields the minimum, Hausdorff, HD95, and ASSD distances per structure pair, deterministically. Sub-sampling is avoided by design: dropping points can only remove the closest pair, never create a closer one, so any thinning biases the minimum distance **upwards** — i.e. it would over-report OAR clearance in the unsafe direction. As a safety net, only clouds exceeding **50,000 points** are thinned to that cap, with a fixed seed so results remain reproducible; typical structures (hundreds to a few tens of thousands of contour points) are processed in full.

## Limitations

- **No inter-structure overlap metrics:** Distance/volume work from contour points (the reported volume is the planimetric slice-stack volume). Sphericity/solidity additionally rasterise each structure onto its **own** local voxel grid, but overlap measures such as **Dice / Jaccard** need both structures on a *common* grid and are not implemented. The distance metrics cannot substitute for them (see *Minimum Distance*).
- **Vertex-based distances:** All distances are computed between contour vertices, a discrete sample of the surfaces, so they approximate the continuous surface distances only to within roughly the vertex and slice spacing.
- **Surface approximation:** Sphericity's surface area is taken from a marching-cubes mesh of the rasterised voxel mask. The mask resolution is bounded for performance, so the surface is a discretised approximation that over-estimates the true area ($\Psi$ biased low, see *Sphericity*); very thin or sub-voxel features may be under-resolved.
- **Planimetric volume:** $\Delta z$ is the nominal slice spacing (median of the smallest cluster of z-differences), so a **gap in z** (a skipped slice, or a union of lesions at different heights stored in one ROI) neither inflates the spacing nor is bridged. Two spheres of radius 10 mm with centres 60 mm apart, stored in one ROI, measure exactly twice one sphere and are reported with `n_components = 2` and `n_gaps = 1`; a sphere with one skipped slice loses that slab (about 7.5 % for a 10 mm radius on 1 mm slices) and is likewise flagged with a gap, whereas a TPS would interpolate the missing slice. A ROI with genuinely variable slice thickness is under-read in its thicker part. Nested inner contours are holes (XOR rule, see *Volume Calculation*); the voxel cross-check in `statistics.txt` uses the same XOR rule and the same nominal spacing.
- **Not a clinical diagnostic tool:** The script serves geometric analysis and does not replace clinical evaluation by a medical physicist or radiation therapist.

## Used Libraries

| Library | Version | Purpose |
|---|---|---|
| `pydicom` | ≥ 2.3 | Reading DICOM files |
| `numpy` | ≥ 1.21 | Numerical calculations and linear algebra |
| `scipy` | ≥ 1.7 | KD-tree, Hausdorff distance, convex hull, interpolation, rotations |
| `shapely` | ≥ 1.8 | 2D polygon operations (validation, area, centroid, XOR of same-plane contours) |
| `matplotlib` | ≥ 3.5 | Static plotting and visualisation |
| `scikit-image` | ≥ 0.19 | Marching cubes surface extraction (analyzer sphericity surface; modifier & case_modifier CT surfaces) |
| `plotly` | ≥ 5.0 | Interactive 3D visualisation (modifier, case_modifier) |

---

# CT Rigid Body Transformer Documentation

## Overview

The **CT Rigid Body Transformer** (`modifier.py`) applies a rigid body transformation — consisting of a translation in three spatial directions and a rotation around three spatial axes — to a CT DICOM series. The result is a new DICOM series for import into a treatment planning system (TPS). The central design goal is to preserve the original Hounsfield Unit (HU) values as accurately as possible while guaranteeing that no geometric distortion of the image is introduced.

## Data Foundation: CT DICOM Geometry

### Patient Coordinate System

DICOM defines a fixed, right-handed **patient coordinate system** (LPS):

| Axis | Direction |
|---|---|
| X | Increases to the patient's **left** |
| Y | Increases **posteriorly** (towards the patient's back) |
| Z | Increases **superiorly** (towards the patient's head) |

All positions and distances in the DICOM standard are specified in millimetres within this coordinate system.

### Relevant DICOM Tags per CT Slice

| Tag | Name | Content |
|---|---|---|
| `(0020,0037)` | `ImageOrientationPatient` (IOP) | Six direction cosines defining row and column orientation |
| `(0020,0032)` | `ImagePositionPatient` (IPP) | 3D position of the centre of the first transmitted pixel (row 0, col 0) in mm |
| `(0028,0030)` | `PixelSpacing` | In-plane pixel size [row spacing, col spacing] in mm |
| `(0028,1053)` | `RescaleSlope` | Linear HU conversion: HU = stored × slope + intercept |
| `(0028,1052)` | `RescaleIntercept` | See above |
| `(0028,0103)` | `PixelRepresentation` | 0 = unsigned, 1 = signed integer |

### Voxel-to-Patient Affine Matrix

The spatial position of any voxel $(k, j, i)$ — where $k$ is the slice index, $j$ the row index, and $i$ the column index — in patient coordinates is given by the following affine transformation:

```math
\begin{pmatrix} x \\ y \\ z \\ 1 \end{pmatrix} = \mathbf{A} \begin{pmatrix} k \\ j \\ i \\ 1 \end{pmatrix}
```

where the $4 \times 4$ affine matrix $\mathbf{A}$ is constructed from the DICOM tags as:

```math
\mathbf{A} = \begin{pmatrix} n_x \cdot \Delta z & F_4 \cdot \Delta r & F_1 \cdot \Delta c & \text{IPP}_x \\ n_y \cdot \Delta z & F_5 \cdot \Delta r & F_2 \cdot \Delta c & \text{IPP}_y \\ n_z \cdot \Delta z & F_6 \cdot \Delta r & F_3 \cdot \Delta c & \text{IPP}_z \\ 0 & 0 & 0 & 1 \end{pmatrix}
```

Here $\mathbf{F} = (F_1, F_2, F_3, F_4, F_5, F_6)$ is the `ImageOrientationPatient` vector, $(F_1, F_2, F_3)$ are the direction cosines of the row direction (increasing column index) and $(F_4, F_5, F_6)$ are the direction cosines of the column direction (increasing row index). The slice normal $\mathbf{n} = (n_x, n_y, n_z) = (F_1, F_2, F_3) \times (F_4, F_5, F_6)$ is computed as the cross product of the two IOP vectors. $\Delta z$ is the slice spacing, $\Delta r$ the row pixel spacing, and $\Delta c$ the column pixel spacing.

The inverse $\mathbf{A}^{-1}$ maps patient coordinates back to voxel indices and is used during resampling.

**Slice ordering assumption.** The loader sorts slices by increasing IPP z and takes the slice axis of $\mathbf{A}$ along $\mathbf{n}$. The two agree only when $\mathbf{n}$ points superiorly, as it does for head-first acquisitions (HFS `[1,0,0,0,1,0]`, HFP `[-1,0,0,0,-1,0]`). For feet-first series (FFS `[-1,0,0,0,1,0]`, FFP `[1,0,0,0,-1,0]`) $\mathbf{n}$ points inferiorly, so the slice axis of $\mathbf{A}$ is reversed and the modelled geometry is wrong; in a synthetic test the computed volume centre lay outside the volume. Such series are currently **not supported**.

### HU Conversion

Stored integer pixel values are converted to Hounsfield Units via a linear mapping defined per slice:

```math
\text{HU} = \text{stored} \times \text{RescaleSlope} + \text{RescaleIntercept}
```

For modern CT scanners the slope is typically 1 and the intercept −1024. By definition of the Hounsfield scale, water is 0 HU and air −1000 HU. With the common 12-bit storage (`BitsStored = 12`, 4096 levels) this gives a range of −1024 to +3071 HU. Cortical bone lies well inside that range, while metal implants usually saturate at the upper end unless the scanner's extended HU scale is used. Values below −1024 HU can occur with signed pixel data, e.g. as padding outside the reconstruction circle.

## Rigid Body Transformation

### Definition

A rigid body transformation in three-dimensional Euclidean space preserves all pairwise distances and angles. It comprises exactly six degrees of freedom: three translational $(t_x, t_y, t_z)$ and three rotational $(r_x, r_y, r_z)$. Formally it is an element of the special Euclidean group $SE(3)$:

```math
\mathbf{T} : \mathbf{p} \mapsto \mathbf{R}\,\mathbf{p} + \mathbf{t}
```

where $\mathbf{R} \in SO(3)$ is a $3 \times 3$ rotation matrix satisfying $\mathbf{R}^T \mathbf{R} = \mathbf{I}$ and $\det(\mathbf{R}) = +1$, and $\mathbf{t} \in \mathbb{R}^3$ is the translation vector. In homogeneous coordinates this becomes the $4 \times 4$ matrix:

```math
\mathbf{T} = \begin{pmatrix} \mathbf{R} & \mathbf{t} \\ \mathbf{0}^T & 1 \end{pmatrix}
```

### Rotation Matrix Construction

The rotation matrix is constructed from three rotation angles using **intrinsic Euler angles in XYZ order** (SciPy's uppercase convention): each successive rotation acts about the axes of the already-rotated, body-fixed frame — first by $r_x$ about the X axis, then by $r_y$ about the new Y axis, then by $r_z$ about the resulting Z axis. The combined rotation matrix is:

```math
\mathbf{R} = \mathbf{R}_x(r_x)\,\mathbf{R}_y(r_y)\,\mathbf{R}_z(r_z)
```

with the elementary rotation matrices:

```math
\mathbf{R}_x(\alpha) = \begin{pmatrix} 1 & 0 & 0 \\ 0 & \cos\alpha & -\sin\alpha \\ 0 & \sin\alpha & \cos\alpha \end{pmatrix}, \quad \mathbf{R}_y(\beta) = \begin{pmatrix} \cos\beta & 0 & \sin\beta \\ 0 & 1 & 0 \\ -\sin\beta & 0 & \cos\beta \end{pmatrix}, \quad \mathbf{R}_z(\gamma) = \begin{pmatrix} \cos\gamma & -\sin\gamma & 0 \\ \sin\gamma & \cos\gamma & 0 \\ 0 & 0 & 1 \end{pmatrix}
```

These are active rotations following the right-hand rule: a positive angle turns counter-clockwise when viewed from the positive end of the axis towards the origin. For example, a positive $r_z$ turns the +X (left) axis towards +Y (posterior).

Equivalently, this is the same rotation matrix as an *extrinsic* ZYX rotation about the fixed patient axes. For a single non-zero angle the intrinsic and extrinsic conventions coincide; they differ only when two or more rotation angles are applied simultaneously.

The implementation uses `scipy.spatial.transform.Rotation.from_euler("XYZ", ...)`, where in SciPy **uppercase letters indicate the intrinsic convention** (lowercase letters would select extrinsic).

### Rotation Centre

The rotation is performed about the **geometric centre** of the CT volume in patient coordinates:

```math
\mathbf{c} = \mathbf{A}\, \begin{pmatrix} (N_z-1)/2 \\ (N_y-1)/2 \\ (N_x-1)/2 \\ 1 \end{pmatrix}
```

Rotating around the volume centre ensures that the patient body remains approximately centred within the output voxel grid and does not drift outside the field of view for small angles. The full forward transformation applied to a point $\mathbf{p}$ is therefore:

```math
\mathbf{p}' = \mathbf{R}\,(\mathbf{p} - \mathbf{c}) + \mathbf{c} + \mathbf{t} = \mathbf{R}\,\mathbf{p} + \underbrace{(-\mathbf{R}\,\mathbf{c} + \mathbf{c} + \mathbf{t})}_{\mathbf{t}_{\text{eff}}}
```

which is stored compactly in the upper-right column of the $4 \times 4$ matrix $\mathbf{T}$.

## Resampling Method (`--method resample`)

### Principle: Inverse Mapping

When a rigid body transformation moves a patient, a new CT scan of the patient in their new position is simulated. Note that the *entire* image content is transformed, including the couch and any immobilisation devices. The result therefore corresponds to moving the whole scene, not the patient relative to a fixed couch, and the couch appears shifted or tilted in the output. The standard approach in medical image processing is **inverse mapping** (also called pull-back or backward mapping): instead of pushing each source voxel into the output grid — which would leave holes — each output voxel asks "where in the source volume did this intensity come from?"

For each output voxel at voxel index $(k, j, i)$, the corresponding patient position is $\mathbf{p}_{\text{out}} = \mathbf{A}[k,j,i,1]^T$. The source position in the original (untransformed) volume is:

```math
\mathbf{p}_{\text{in}} = \mathbf{T}^{-1}\,\mathbf{p}_{\text{out}}
```

Converting back to voxel indices: $\mathbf{v}_{\text{in}} = \mathbf{A}^{-1}\mathbf{p}_{\text{in}}$. Combining these steps, the complete voxel-to-voxel mapping is:

```math
\mathbf{v}_{\text{in}} = \underbrace{\mathbf{A}^{-1} \mathbf{T}^{-1} \mathbf{A}}_{\mathbf{M}}\, \mathbf{v}_{\text{out}}
```

The matrix $\mathbf{M}$ is computed once and then applied to the coordinates of every output voxel (chunk by chunk, see below). The source coordinates $\mathbf{v}_{\text{in}}$ are generally non-integer; the HU value is obtained by **interpolation** in the source volume.

Voxels whose source coordinates fall outside the original volume are assigned −1000 HU (air). This is a fill convention for regions that were never imaged, which may in reality contain anatomy. For example, after a z-shift the vacated slab at one end of the scan range becomes air, although the patient continues beyond it. This matters for dose calculation if beams pass through such regions.

### Why Inverse Rather Than Forward Mapping?

Forward mapping (pushing each source voxel to its new position) suffers from two problems: (1) output voxels that no source voxel maps to remain unfilled (holes), and (2) output voxels that multiple source voxels map to require a compositing strategy. Inverse mapping avoids both problems by construction and is the standard in medical image registration (e.g. SPM, FSL, ITK).

### Memory-Efficient Chunk Processing

For a typical CT volume of $512 \times 512 \times 320$ voxels, a whole-volume homogeneous coordinate array (shape $4 \times N$, float64) would occupy approximately 2.5 GiB of RAM. The inverse mapping needs two such arrays at the same time (the output-voxel coordinates and their mapped source coordinates), plus the integer index grids and SciPy's internal copy of the coordinates. The implementation therefore processes the volume in **chunks of 20 slices at a time**, generating the coordinate arrays only for the current chunk. Each per-chunk $4 \times N$ array occupies approximately 160 MiB, and the measured transient working set per chunk, all coordinate-related arrays included, is about 0.55 GiB. The volume data comes on top of that: the float32 input and output plus a float64 working copy, about 1.25 GiB for $512 \times 512 \times 320$.

### Interpolation

The source-volume lookup at non-integer coordinates requires interpolation. Three orders are available:

#### Nearest Neighbour (order 0)

```math
\text{HU}(\mathbf{v}) = \text{HU}\bigl(\text{round}(\mathbf{v})\bigr)
```

No weighted averaging; the nearest voxel value is taken directly. This guarantees that **only HU values that actually exist in the original volume appear in the output** (plus the −1000 HU fill outside it). The disadvantages are blocky staircase artefacts at structure boundaries and a positional error of up to half a voxel along each axis. It suits cases where the exact discrete values must be kept, e.g. synthetic phantoms or label-like images with only a few distinct values. Dose calculation does not require it: a TPS's HU-to-density calibration curve handles interpolated HU values without difficulty.

#### Trilinear Interpolation (order 1, default)

For a point at fractional position $(k + \delta_k, j + \delta_j, i + \delta_i)$ with $\delta \in [0,1)$, the value is the weighted average of the eight surrounding voxels:

```math
\text{HU}(\mathbf{v}) = \sum_{a \in \{0,1\}} \sum_{b \in \{0,1\}} \sum_{c \in \{0,1\}} w_{abc}\,\text{HU}(k{+}a,\, j{+}b,\, i{+}c)
```

with trilinear weights $w_{abc} = |1-\delta_k-a|\cdot|1-\delta_j-b|\cdot|1-\delta_i-c|$. Trilinear interpolation is continuous and cannot overshoot the original value range: each interpolated value is a convex combination of the eight surrounding voxels, so it lies between their minimum and maximum. Where HU varies smoothly, the interpolation error is typically **< 2 HU** (well below 1 HU in a synthetic soft-tissue test); larger errors are confined to about one voxel at sharp boundaries. Regional mean HU values, which the TPS density lookup and the dose calculation effectively depend on, are preserved. Individual voxel values do change, however. The weighted averaging acts as a mild low-pass filter that reduces uncorrelated image noise (σ = 12 HU → 6.6 HU in a single pass of the test) and slightly blurs edges. A voxel-by-voxel comparison of noisy CT therefore shows differences comparable to the noise level, not < 2 HU. Trilinear interpolation is a common default in image registration software; see [Lehmann et al. 1999] and [Thévenaz et al. 2000] for surveys of interpolation methods in medical imaging.

#### Tricubic B-Spline Interpolation (order 3)

A $C^2$-continuous piecewise-cubic **interpolant** that passes exactly through the original voxel values. Its B-spline coefficients are computed by a separable recursive prefilter applied once to the entire volume before sampling. It blurs less than trilinear interpolation, preserving fine detail and noise texture better (σ = 12 HU → 9.8 HU in the same test), and it has the smallest error in smooth regions away from high-contrast edges. It is, however, not bounded by the neighbouring values: near sharp boundaries it overshoots (Gibbs-like ringing) and produces HU values outside the original range. The overshoot is about 11 % of the step height for an ideal step edge. For the PSF-blurred edges of real CT it is a few percent; a synthetic test with realistic edge blur went about 10–20 HU beyond the input range. Recommended when interpolation quality matters more than strict range preservation and the extra computation and memory (a float64 coefficient volume) are acceptable.

### HU Value Preservation: Verification

`modifier.py` prints the minimum and maximum HU values of the original and transformed volumes as a quick sanity check (it does not fail on a mismatch):

```
HU range original:     [-1024, 3071] HU
HU range transformed:  [-1024, 3071] HU
```

For trilinear (order 1) interpolation the output is guaranteed **not to exceed** the original range (no overshoot); the extreme values may shrink slightly, and regions without source data are filled with −1000 HU. For cubic (order 3) interpolation the range can widen by the overshoot described above, typically by 10–20 HU at the PSF-blurred edges of real CT. Nearest-neighbour (order 0) interpolation introduces no new values.

Note that a voxel-wise difference map between original and transformed volumes is **not a meaningful quality metric** in this context: after a rigid body motion the same tissue appears at different voxel positions in the two volumes, so direct subtraction compares different anatomical structures. A correct interpolation quality assessment requires a round-trip test (apply T then T⁻¹) or comparison within a co-registered reference frame.

### Output DICOM Structure

The resampled volume is written back to DICOM using the metadata of the original slices. The output geometry (IPP, IOP, pixel spacing, slice spacing) is **identical** to the input, meaning the output slices occupy the same spatial positions as the original. The transformed patient body appears shifted/rotated *within* the fixed voxel grid. Every output voxel shows whatever lay at its source position: regions the body moved out of therefore show the surrounding air, or the −1000 HU fill where the source position is outside the original volume.

The HU-to-stored conversion is the exact inverse of the loading step:

```math
\text{stored} = \text{round}\!\left(\frac{\text{HU} - \text{RescaleIntercept}}{\text{RescaleSlope}}\right)
```

Values are clipped to the valid int16 range $[-32768, 32767]$ and stored as signed 16-bit integers (`PixelRepresentation = 1`), which is common for CT. (Some scanners instead store unsigned pixels with `PixelRepresentation = 0` and a matching intercept; the output here is written as signed int16 regardless.) The original `RescaleSlope` and `RescaleIntercept` are preserved unchanged so the output is correctly calibrated in any TPS.

Each output file receives a new `SOPInstanceUID` and `SeriesInstanceUID` (generated with `pydicom.uid.generate_uid()` which produces DICOM-conformant UID strings) so that the transformed series is recognised as an independent series by the TPS and PACS, while the `StudyInstanceUID` and patient demographics remain identical for correct study association. The new `SOPInstanceUID` is also written to the file-meta header (`MediaStorageSOPInstanceUID`), because DICOM Part 10 requires the two to match.

## Metadata-Only Method (`--method metadata`)

### Principle

For the metadata-only approach, **no pixel data is modified**. Instead, the spatial meaning of each slice is updated by transforming the DICOM positional tags:

**ImagePositionPatient (IPP):** The position of the first pixel of each slice is mapped through the forward transformation:

```math
\text{IPP}'_k = \mathbf{T}\,\begin{pmatrix}\text{IPP}_k \\ 1\end{pmatrix}
```

**ImageOrientationPatient (IOP):** The six direction cosines (row direction and column direction) are rotated by $\mathbf{R}$:

```math
\mathbf{F}'_{\text{row}} = \mathbf{R}\,\mathbf{F}_{\text{row}}, \qquad \mathbf{F}'_{\text{col}} = \mathbf{R}\,\mathbf{F}_{\text{col}}
```

Since $\mathbf{R} \in SO(3)$ preserves norms, the transformed direction cosines remain unit vectors and their cross product continues to define the correct slice normal.

### Advantages and Limitations

The metadata-only method guarantees **exact HU preservation** because no interpolation is performed. It is also significantly faster since no resampling loop is needed. However, after an arbitrary rotation the IOP vectors are no longer aligned with the standard axial orientation `[1,0,0,0,1,0]`. Many treatment planning systems require the planning CT to be axial. They reject oblique datasets, or accept them only as secondary images for registration. Such a check rejects *any* non-axial orientation regardless of the rotation magnitude, so a small angle is no guarantee of acceptance. Support differs between systems and versions, so verify it for your TPS; when in doubt, use `--method resample` for guaranteed-axial output.

## No-Distortion Guarantee

A geometric distortion would occur if different regions of the anatomy were scaled, sheared, or mapped non-linearly. The transformation used here is strictly rigid: $\mathbf{R}$ is orthonormal ($\det \mathbf{R} = 1$, $\lVert\mathbf{R}\mathbf{v}\rVert = \lVert\mathbf{v}\rVert$ for all $\mathbf{v}$), so all distances and angles between any two anatomical points are preserved exactly **in patient space** — the space in which the no-distortion guarantee is defined. The voxel-index mapping $\mathbf{M} = \mathbf{A}^{-1}\mathbf{T}^{-1}\mathbf{A}$, by contrast, has an orthonormal linear part only for isotropic voxels: with anisotropic spacing (e.g. 1 × 1 mm in-plane, 2 mm slices) a physical rotation legitimately appears as a combination of rotation, scaling, and shear in index space, because the same millimetre displacement spans a different number of voxels along each axis. This is a property of the grid representation, not an anatomical distortion. Interpolation is needed in any case, because the mapped source positions are generally non-integer; it affects the image *sharpness* (see *Interpolation*), not the geometry.

## 3D Visualisation

### Surface Extraction: Marching Cubes

The 3D body surface is extracted from the CT volume using the **Marching Cubes algorithm** (Lorensen & Cline, 1987). For each $2 \times 2 \times 2$ cube of adjacent voxels, the algorithm determines which of the 256 possible configurations of voxels above/below the iso-threshold is present and places triangular surface patches accordingly. This produces a polygonal mesh approximating the iso-surface at the chosen HU threshold.

Two surfaces are extracted:
- **Body surface**: threshold −300 HU, separating soft tissue from air
- **Bone surface**: threshold +400 HU, isolating cortical bone and other dense material. Cortical bone is typically ≳ +700 HU. Most cancellous (trabecular) bone, at roughly +100 to +400 HU, lies below this level, while strongly contrast-enhanced vessels and metal lie above it.

A **downsampling factor of 2** is applied before running Marching Cubes on the body surface (every second voxel in each direction) to reduce computation time and triangle count; the bone surface uses a factor of 3. This lowers the spatial resolution of the surface mesh but has no effect on the underlying DICOM data.

### Coordinate Conversion

The Marching Cubes algorithm returns vertices in voxel index space $(k, j, i)$. These are converted to patient coordinates in millimetres by the affine matrix $\mathbf{A}$:

```math
\mathbf{p}_{\text{patient}} = \mathbf{A}\, \begin{pmatrix} k \\ j \\ i \\ 1 \end{pmatrix}
```

For the metadata-only method, where pixel data is unchanged, the transformed surface is obtained by applying the forward transformation $\mathbf{T}$ to the original vertices directly — no second surface extraction from a resampled volume is needed.

### Interactive Visualisation

The extracted meshes are rendered using **Plotly's Mesh3d** trace with Gouraud-style lighting. The layers (original body, transformed body, original bone, rotation centre) are independently toggleable via the legend; the coordinate axes at the rotation centre are always shown. The scene uses `aspectmode='data'` to ensure that spatial distances are displayed without distortion, i.e., 1 mm in X, Y, and Z corresponds to the same pixel length on screen.

## Summary of Design Decisions

| Decision | Choice | Rationale |
|---|---|---|
| Mapping direction | Inverse (pull-back) | Avoids holes and compositing problems inherent to forward mapping |
| Default interpolation | Trilinear (order 1) | Common default; no overshoot; regional mean HU preserved, small error in smooth tissue |
| Rotation convention | Intrinsic XYZ (≡ extrinsic ZYX) | SciPy's uppercase `"XYZ"`; coincides with the extrinsic convention for any single-axis rotation, so an $r_z$-only correction is always a pure axial-plane rotation |
| Rotation centre | Geometric centre of the CT grid | Body stays within the output FOV for typical small-angle corrections |
| Output geometry | Same grid as input | Standard axial IOP preserved (for axial input); broadly TPS-compatible |
| HU back-conversion | Per-slice slope/intercept | Preserves original calibration; no re-calibration artefacts |
| UIDs | New Series + SOP UIDs | TPS/PACS recognises output as independent series; no confusion with original |

---

# Case Modifier Documentation

## Overview

The **Case Modifier** (`case_modifier.py`) extends the CT Rigid Body Transformer to operate on a **complete case** — both the CT series and its companion RTSTRUCT — in a single invocation. Where `modifier.py` only moves the pixel data, `case_modifier.py` additionally transforms every contour point in the structure set so that ROIs continue to enclose the same anatomy after rigid motion. UID references between RS and CT are rewritten so the transformed pair links to each other.

The module is built on top of `modifier.py` and `analyzer.py`; it adds no new geometric mathematics, only orchestration, RTSTRUCT-specific bookkeeping, and metadata management.

## Why the RTSTRUCT Must Move with the CT

A DICOM RTSTRUCT stores ROI contours as flat lists of (x, y, z) coordinates in the **patient frame** (LPS, mm). When the CT is rigidly transformed, the patient anatomy at original coordinate $\mathbf{p}$ now appears at $\mathbf{p}' = \mathbf{T}\mathbf{p}$. A contour that was drawn on a structure at $\mathbf{p}$ remains stored at $\mathbf{p}$ in the RTSTRUCT — so without a corresponding RS update, the contour no longer encloses the moved anatomy. After applying the same $\mathbf{T}$ to every `ContourData` triple, the contours follow the anatomy exactly.

Because $\mathbf{T}$ is a rigid transformation, $\lVert\mathbf{R}\mathbf{p}_1 - \mathbf{R}\mathbf{p}_2\rVert = \lVert\mathbf{p}_1 - \mathbf{p}_2\rVert$ for all pairs of points, and consequently:

- Contour shapes are not distorted (no shearing, no anisotropic scaling).
- The **true 3D volume** is preserved exactly, since a rigid transform has $|\det \mathbf{R}| = 1$. Note, however, that the volume *re-measured* afterwards by the planimetric Shoelace + mean-slice-spacing formula is invariant only for translations and pure Z-rotations. For X/Y rotations that tilt the contour planes by an angle $\theta$, it shrinks by roughly $\cos^2\theta$ (−0.8 % at 5°, −3 % at 10°, −7 % at 15°). The Shoelace formula measures the area projected onto the axial plane ($A\cos\theta$), and the z-spacing of the tilted contours shrinks by the same factor. This is an artefact of the per-axial-slice convention (see [The `resample` vs `metadata` Trade-Off](#the-resample-vs-metadata-trade-off)), not a change of the true volume.
- Centroids transform linearly: $\mathbf{T}(\overline{\mathbf{p}}) = \overline{\mathbf{T}(\mathbf{p}_i)}$, which is exploited by `--verify` as a per-ROI sanity check.

The implementation reads `ContourData`, reshapes the flat list into an $N \times 3$ matrix, applies $\mathbf{T}$ to every row in homogeneous coordinates, formats the result back into the flat string list with six-decimal precision, and writes it back. `NumberOfContourPoints` is unchanged.

## UID Bookkeeping: Why This Is the Hard Part

A DICOM RTSTRUCT references its companion CT through three layers of UIDs:

1. `FrameOfReferenceUID` (top level + per-ROI in `StructureSetROISequence`) — the patient coordinate system.
2. `ReferencedFrameOfReferenceSequence[*].RTReferencedStudySequence[*].RTReferencedSeriesSequence[*].SeriesInstanceUID` — the CT series the RS belongs to.
3. Per-contour `ContourImageSequence[*].ReferencedSOPInstanceUID` — the specific CT slice each contour was drawn on.

When `save_ct_series` writes the transformed CT, every output slice receives a fresh `SOPInstanceUID` and the series receives a fresh `SeriesInstanceUID`. Without rewriting, the original RS would now reference SOPs and a series UID that no longer exist — the planning system would report "missing image references" or refuse to link the structure set at all.

`save_ct_series` returns, among other things, a mapping `sop_map = {old_sop: new_sop, ...}` for all transformed slices. `transform_rtstruct` walks the RS and substitutes:

- Per-contour `ContourImageSequence[*].ReferencedSOPInstanceUID` via the SOP map (raises `KeyError` if a referenced SOP is not in the map — this is a hard error, since it means the RS is referencing slices outside the input CT folder).
- Top-level `RTReferencedSeriesSequence[*].SeriesInstanceUID` to the new CT series UID.
- Top-level `RTReferencedSeriesSequence[*].ContourImageSequence[*].ReferencedSOPInstanceUID` via the SOP map.

The RS itself receives a fresh `SOPInstanceUID` (again also in the file-meta header) and `SeriesInstanceUID` so it is recognised as a new structure set.

## FrameOfReferenceUID Strategy

The `FrameOfReferenceUID` (DICOM tag `(0020,0052)`) declares "all series with this UID share the same patient coordinate system." Two strategies are supported:

- **Default (keep):** the original FoR is preserved on both the transformed CT and the transformed RS. They link correctly to each other because both carry the same FoR. **Caveat:** any existing RTPLAN, RTDOSE, or sibling RTSTRUCT that also references this FoR will be auto-overlaid by Aria/Eclipse onto the transformed CT, even though they were planned in the un-transformed coordinate system. The case modifier prints a clear runtime warning whenever this default is used.
- **`--new-frame-of-reference`:** a fresh FoR UID is minted and applied to both the transformed CT and the transformed RS. They still link to each other, but they are now in a coordinate system distinct from the original CT's. Existing plans/doses keep their old FoR and are not auto-overlaid. Any cross-frame comparison must go through an explicit registration object — this is the safer choice for clinical workflows where the transformed dataset is meant to represent a different geometric situation rather than an alternate view of the same one. It is **not** the tool's default, however: the default is **keep**, for backward compatibility with existing plan/dose linkage.

Note that *new* is also the choice consistent with DICOM semantics. Series that share a `FrameOfReferenceUID` are declared spatially related in one patient coordinate system. Because the transformed anatomy no longer coincides with the original, keeping the FoR asserts an identity registration that does not hold. *keep* is a pragmatic compatibility default, not the DICOM-conformant one.

## Variable Rotation Centre

Rotation is performed about a configurable centre $\mathbf{c}$, with the offset folded into the same $4 \times 4$ matrix used for the CT:

```math
\mathbf{p}' = \mathbf{R}\,(\mathbf{p} - \mathbf{c}) + \mathbf{c} + \mathbf{t}
```

Three centre-selection modes are provided:

- **`--center volume`** (default for non-interactive runs) — the geometric centre of the CT volume in patient coordinates. Same as `modifier.py`.
- **`--center marker:NAME`** — the position of a POINT-type ROI from the RTSTRUCT, looked up case-insensitively. POINT contours are typically used to mark fiducials, isocentres, or registration points (`ContourGeometricType == "POINT"`, single (x,y,z) triplet in `ContourData`). Listing them with `--list-markers` gives the user the marker names and their LPS positions.
- **`--center x,y,z`** — three comma-separated floats, interpreted as LPS millimetres directly. Useful for arbitrary offsets that do not correspond to an existing marker.

When neither `--center` nor `--non-interactive` is given and stdin is a TTY, the user is presented with a list of markers and may pick one by index, type `v` for the volume centre, or `m` to enter a manual point.

## Drehpunkt Marker

After every successful run the transformed RS contains an extra POINT-type ROI named **`Drehpunkt`** (German for "pivot point" / "rotation centre") at the position $\mathbf{c} + (t_x, t_y, t_z)$. This works because rotation leaves $\mathbf{c}$ invariant — $\mathbf{T}(\mathbf{c}) = \mathbf{R}(\mathbf{c} - \mathbf{c}) + \mathbf{c} + \mathbf{t} = \mathbf{c} + \mathbf{t}$ — so the rotation centre in the transformed coordinate system simply shifts by the translation. Visualised in the TPS, this marker shows at a glance where the rigid motion was anchored, which is otherwise not derivable from the DICOM tags alone.

Implementation details: a new entry is appended to `StructureSetROISequence` (next free `ROINumber`, `ROIName = "Drehpunkt"`), to `RTROIObservationsSequence` (`RTROIInterpretedType = "MARKER"`), and to `ROIContourSequence` (one POINT contour with `NumberOfContourPoints = 1`, yellow display colour). The marker's `ReferencedFrameOfReferenceUID` matches the chosen FoR strategy.

## Aria/TPS-Visible Metadata

To make the transformed series identifiable in the planning system without opening the DICOM headers, the case modifier writes the transform parameters into the human-readable string tags:

| Tag | VR | Max | Content |
|---|---|---|---|
| `SeriesDescription` (CT + RS) | LO | 64 | `<orig>_RB` (or custom suffix) |
| `StructureSetLabel` | SH | 16 | `<orig>_RB` truncated; suffix preserved if `<orig>` is too long |
| `StructureSetName` | LO | 64 | `<orig>_RB` |
| `StructureSetDescription` | ST | 1024 (tool truncates to 64) | `rigid t=(tx,ty,tz) r=(rx,ry,rz) c=<centre> m=<method> FoR=<keep\|new>` |
| `SeriesNumber` | IS | — | original + 1000 |

The label suffix is configurable via `--label`. The transform description always lands in `StructureSetDescription` regardless of the suffix. It names the rotation centre (`Marker 'HS1'`, `Volumenmitte`, `manuell`, `interaktiv`) but does not give its coordinates; the position is recorded by the `Drehpunkt` marker at $\mathbf{c} + \mathbf{t}$.

## Pre-Flight Validation

Before any data is written, the case modifier validates the input:

1. **Folder layout.** The case directory must exist and contain a `CT/` subfolder. Exactly one `RS*.dcm` file must be present at the case root (or `--rs PATH` must be given). Sibling `RP*.dcm`/`RD*.dcm` files trigger a hint warning.
2. **CT geometry.** All slices must share the same `ImageOrientationPatient` (within 1e-3) and `PixelSpacing` (within 1e-4 mm). Slice spacing must be uniform within 1 % relative deviation. At least 2 slices are required. The orientation is only checked for consistency across slices, not for being axial (see *Limitations and Out-of-Scope*).
3. **FoR consistency.** The CT slices must all carry the same `FrameOfReferenceUID`, and the RTSTRUCT must reference this FoR through its `ReferencedFrameOfReferenceSequence`. A mismatch raises a clear error — almost always indicating that the user picked an RS file that does not belong to this CT.

Any failure terminates the run with exit code 2 before any output is written.

## The `resample` vs `metadata` Trade-Off

Rigid-body rotation of a CT + RTSTRUCT pair has a subtle geometric implication that depends on which method is used. The contour points themselves are always transformed correctly by the same $\mathbf{T}$ as the CT — point-to-point distances are preserved exactly, which is the definitive rigid-body invariant. The trade-off is about how the transformed dataset is *represented*:

**`--method metadata`:** the CT's `ImagePositionPatient`/`ImageOrientationPatient` are rotated by $\mathbf{R}$, so the slice planes themselves become tilted in patient space. The contour points (which were originally on axial slices at constant $z$) are also rotated; their new positions lie exactly on the new tilted slice planes. **The relative geometry between CT slices and RS contours is preserved exactly.** No clipping, no discretisation artefacts.
- ✓ True rigid-body fidelity. Pure 3D mathematics, no resampling.
- ✓ HU values byte-identical to the source.
- ✗ Output IOP is no longer `[1,0,0,0,1,0]`. Many TPS require an axial planning CT, and some PACS viewers mishandle oblique CT/RTSTRUCT; check your system.

**`--method resample` (default):** the CT's voxel grid stays axial, and pixels are inverse-mapped from the source. The contour points are still rotated by $\mathbf{T}$ in 3D, which means **after a non-axial rotation, the contours are tilted polygons that no longer lie on the axial slice planes** of the new CT. DICOM RTSTRUCT is conventionally per-axial-slice, so this is technically off-spec. Two consequences:
- *Clipping at the grid boundaries.* The output grid is fixed, so anatomy moved beyond it is lost. Grid regions whose source lies outside the original scan are filled with −1000 HU. A contour point that stays inside the grid always has the correct anatomy beneath it, because it samples exactly its original source position. A contour point moved *outside* the grid has no image data at all. The case modifier detects such points automatically and prints a per-ROI warning at the end of the run.
- *Volume measurements via the per-axial-slice formula (XOR area × nominal slice spacing) are no longer invariant.* The transformation itself preserves volume in 3D, but a tool that interprets the output RS as per-axial-slice contours (including this repo's `analyzer.py`) will measure a volume reduced by roughly $\cos^2\theta$ (−3 % at 10°, −7 % at 15°). This is a measurement artefact of the formula, not a transformation error.
- *Z-only rotations are immune* to the tilt and volume effects, because $r_z$ leaves contour Z values unchanged and each contour stays in its original axial slice. They can still move anatomy out of the in-plane grid when the rotation centre is far from the image centre or a translation is added.

**Practical guidance:**
- Pure translations: contours and anatomy stay exactly aligned with both methods. Resample keeps standard axial output, but on its fixed grid it discards whatever is shifted out and fills the vacated margin with air (for a z-shift, one end of the scan range).
- Z-rotation only: both methods are exactly equivalent in geometric fidelity; resample is preferable for TPS compatibility.
- X/Y rotation: prefer **metadata** for analytic correctness if your TPS accepts oblique CT/RS; otherwise accept the resample artefacts and verify on critical structures.

## Verification Modes

Three orthogonal modes are provided to verify correctness:

- **`--dry-run`** — runs all input validation, computes the transform matrix, resolves the rotation centre, and prints a plan including the resolved centre, the 4×4 matrix $\mathbf{T}$, and the planned output paths. No files are written.
- **`--verify`** — after writing the transformed RS, re-reads it from disk, computes per-ROI vertex centroids (the mean of all contour points), and compares them with $\mathbf{T} \cdot \overline{\mathbf{p}}_{\text{orig}}$. Vertex centroids transform linearly under rigid (indeed any affine) motion, so a deviation well above the ~1e-6 mm rounding floor (e.g. > 1e-4 mm) indicates a bug in the transform pipeline. The maximum norm and the worst ROI are reported; the check is informational and does not change the exit code.
- **`--self-test`** — runs three independent checks, each writing to a temporary directory with `--method metadata` (so the CT resampling path is not exercised), and reports PASS/FAIL each:
  1. *Identity round-trip*: zero translation, zero rotation; every contour point matches the original within 1e-4 mm. Catches accidental side effects in the pipeline.
  2. *Z-rotation pairwise distance drift* (15°): point-to-point distances within each ROI (up to 200 randomly sampled points per ROI) are preserved within 1e-3 mm. Validates the trivially-volume-preserving direction.
  3. *X-rotation pairwise distance drift* (5°): same test on a non-axial rotation. This is the **definitive rigid-body check** — distances are coordinate-system-invariant and are preserved even when the planar Shoelace formula no longer yields invariant volumes.

  Pairwise distances are the right invariant to check: they are preserved by every rigid transform regardless of orientation, whereas the analyzer's planar-axial volume formula is not (see *resample vs metadata trade-off* above).

The case modifier additionally runs an automatic **clipping check** after every real transformation in resample mode: it counts how many transformed contour points fall outside the output voxel grid, prints a per-ROI table sorted by severity, and recommends `--method metadata` when clipping is observed. In metadata mode the check is geometrically void (the slice grid rotates with the contours) and is skipped.

## Numerical Verification

Run on a representative clinical case (a few hundred CT slices, several dozen ROIs including POINT markers), the implementation achieves the following typical precision. The first three rows are reproducible on any case via `--self-test`, the centroid row via `--verify`:

| Test | Result | Tolerance |
|---|---|---|
| Identity round-trip — max contour-point deviation | 5 × 10⁻⁷ mm | < 1 × 10⁻⁴ mm |
| Z-rotation 15° — max pairwise distance drift | 1.3 × 10⁻⁶ mm | < 1 × 10⁻³ mm |
| X-rotation 5° — max pairwise distance drift | 1.3 × 10⁻⁶ mm | < 1 × 10⁻³ mm |
| Centroid linearity under non-trivial $\mathbf{T}$ | < 1 × 10⁻⁶ mm (worst ROI) | < 1 × 10⁻⁴ mm (guideline) |
| Marker fixpoint (rotation about itself) | 0 mm | < 1 × 10⁻⁶ mm |
| Drehpunkt position | exact | — |
| FoR consistency CT ↔ RS (both modes) | preserved | — |

The numerical floor (~1e-7…1e-6 mm) is set by the float-to-six-decimal-string round trip in `ContourData`, not by the transform mathematics. Note that **volume preservation as measured by the planar-axial Shoelace formula is *not* a tolerance to claim** — it holds exactly only for Z-rotations and translations. For X/Y rotations the formula under-reads by roughly $\cos^2\theta$ (see above), an artefact of the per-slice convention, not a transformation error.

## Limitations and Out-of-Scope

- **RTPLAN and RTDOSE are not transformed.** Sibling `RP*.dcm` and `RD*.dcm` files are detected and reported but their geometry is left untouched. Transforming a plan would require also transforming beam isocentres, gantry angles, and couch positions, which is out of scope for this tool.
- **Head-first axial CT input is assumed but not enforced.** `validate_ct_geometry` rejects non-uniform orientation, pixel spacing and slice spacing, but it does not check that the orientation is axial or that the slice positions advance along the slice normal. Two kinds of series therefore pass validation although their voxel geometry is modelled incorrectly: gantry-tilted series, whose slice positions do not advance along the normal, and feet-first series, whose normal points inferiorly (see *Voxel-to-Patient Affine Matrix*). Resample them to a standard head-first axial grid first.
- **Single RTSTRUCT per case.** When multiple `RS*.dcm` files are present, the user must select one with `--rs`.
- **No dose recomputation.** The transformed CT can be re-imported into a TPS for fresh dose calculation, but no dose is recomputed in this tool.
- **No re-projection of tilted contours onto axial slices in resample mode.** A semantically clean solution for X/Y rotations under resample would extract a 3D mesh from the original axial contours, rotate the mesh, and slice it with the new axial planes. This is not implemented; metadata mode is the recommended path when contour-on-axial-slice fidelity matters.

---

# RTSTRUCT Visualizer Documentation

## Overview

The **RTSTRUCT Visualizer** (`visualizer.py`) generates a set of purpose-built static plots from the analysis results produced by the RTSTRUCT Analyzer. Each plot is designed to answer a specific clinical question and avoids chart types that are misleading at the typical structure counts encountered in radiotherapy planning (typically 5–40 structures per case). The visualizer calls `run_analysis` internally, so only the RTSTRUCT file path is required — no intermediate JSON file is needed.

## Output Files and Clinical Motivation

### 1. `volumes.png` – Structure Volume Bar Chart

**Chart type:** Horizontal bar chart on a logarithmic volume axis, grouped by category. Targets come first (GTV/PTV of the same lesion adjacent), then serial OARs, parallel OARs and, hatched, helper structures; each OAR/helper group is sorted by descending volume. The external contour and POINT markers are not shown.

A bar chart is used because each bar represents one named anatomical structure (a histogram would only be meaningful for many samples drawn from an unknown distribution, which is not the case here). The log axis keeps small SRS targets (≈ 0.03 cm³) visible next to the whole brain (> 1000 cm³). The chart makes immediately apparent which structures are largest, how Targets and OARs compare in size, and whether any structure has an unexpectedly small or large volume (which can indicate a contouring error). Colours follow the shared category palette: Targets blue, serial OARs red, parallel OARs orange, helpers grey. Volume values are annotated on each bar.

**Clinical relevance:** Volume is a primary descriptor for both Target coverage and OAR sparing (e.g., mean dose to a parallel organ correlates with the irradiated volume). Unexpected outliers in volume are a common QA flag for contouring errors.

### 2. `shape_metrics.png` – Shape Metrics Heatmap

**Chart type:** Category-grouped heatmap table (rows = anatomical structures, columns = sphericity / solidity / elongation), with a category colour strip on the left.

A heatmap is used because it stays readable at the typical ~40 structures per case (grouped bar charts become an unreadable tangle of rotated x-axis labels at that count) and lets outliers (e.g. a highly elongated spinal cord) be spotted at a glance. Sphericity and solidity use a 0→1 diverging colormap (green = round / convex), elongation a sequential map (darker = more stretched). **Only single-component anatomical structures** (Targets + OARs with `shape_valid = true`) are shown; multi-component helper/union/shell structures are excluded because their convex-hull metrics are meaningless (they remain listed in `statistics.txt`).

**Sphericity** quantifies how closely the structure resembles a sphere: round PTVs score highest, about 0.8–0.9 rather than 1 because of the discretisation bias (see *Sphericity*); complex shapes score lower and may need more beam arrangements. **Solidity** (volume / convex-hull volume) detects concavity. As a rough heuristic, a value below ~0.85 suggests the structure wraps around other anatomy, e.g. a C-shaped PTV around the brainstem (small structures read high, see *Solidity*). **Elongation** (sqrt of largest-to-smallest PCA eigenvalue ratio) measures directional stretching; high elongation with low sphericity in the spinal cord confirms its cylindrical nature, while unexpectedly high elongation in a GTV may indicate a drawing artefact.

### 3. `distances.png` – Critical Target↔OAR Distance Lollipop

**Chart type:** Horizontal lollipop plot — one row per Target↔**serial-OAR** pair, sorted ascending by minimum distance (most critical on top). The filled dot is the minimum distance, the open marker is HD95, and a connecting line spans the two. Dashed lines mark the 3 mm and 5 mm thresholds (pragmatic planning-margin conventions used here as visual guides, not a formal dosimetric standard); rows below 5 mm are highlighted red.

**Design decisions:**
- **Only Target↔serial-OAR pairs** (brainstem, cord, optic nerves, chiasm, pituitary) are shown — these are the proximity-critical, dose-limiting organs. Category-aware filtering keeps out clinically meaningless *containment* pairs. For a PTV and the union structure that contains it, the distance is ~0 mm by construction because they share boundary points. For a PTV inside the whole brain, the surface-to-surface distance says nothing about clearance. If no serial OAR is found, all Target↔OAR pairs are shown instead. Broader proximity is covered by `proximity_matrix.png`.
- **HD95 instead of raw Hausdorff** as the companion metric: the raw maximum Hausdorff is dominated by single outliers, so the 95th-percentile variant is shown for robustness.
- Min/HD95/Hausdorff/ASSD/centroid distances for all pairs remain available in `statistics.txt`.

All distances are computed **exactly** on the full point clouds (KD-tree, deterministic; only clouds above 50,000 points are thinned to that cap with a fixed seed — see the analyzer's *Performance Notes*).

### 4. `centroids_3d.png` – Spatial Centroid Map

**Chart type:** 3D scatter plot using matplotlib's mpl_toolkits.mplot3d.

Structure centroids are plotted in the DICOM patient coordinate system (X=left, Y=posterior, Z=superior), coloured by category (Target / serial OAR / parallel OAR). Marker size scales with the square root of structure volume. Instead of labelling all ~40 centroids (previously an unreadable label tangle), only the serial OARs are labelled and the Targets are numbered 1…N with a side legend mapping number → lesion (POINT markers and the external body contour are excluded). The 3D box aspect ratio is set proportional to the data extent in each axis, so equal millimetre distances render as equal lengths.

**Clinical relevance:** This plot answers the question "where is everything relative to everything else?" at a glance. It is particularly useful for multi-metastasis cases (several GTV/PTV pairs plus critical OARs): one can verify that all GTV/PTV pairs are spatially co-located and that the OARs (brainstem, optic chiasm, cochleae) are in the expected positions.

**Limitation:** The 3D scatter is static (not interactive). The repo's only interactive 3D views are the Plotly HTML files of the modifier (`visualization_3d.html`, CT surfaces) and the case modifier (`transform_3d.html`, contour point clouds); neither shows centroids.

### 5. `statistics.txt` – Numerical Summary

A structured plain-text file, grouped by category (Targets, serial OARs, parallel OARs, helper structures), with per-structure details (contour and slice count, nominal slice spacing with a z-gap flag, volume + voxel cross-check, equivalent-sphere diameter, max 3D diameter, centroid, bounding box, sphericity, solidity, elongation, and a flag for multi-component / invalid shape metrics). It also contains aggregated Target↔OAR distance statistics (mean/std/min/max of min, HD95, Hausdorff and ASSD over *all* Target↔OAR pairs, including containment pairs such as Target↔whole brain), the most critical pairs, and a GTV→PTV margin check (minimum distance, HD95 and ASSD per lesion). Suitable for copy-paste into reports or further spreadsheet analysis.

### 6. Additional clinical plots (multi-metastasis / SRS)

- **`proximity_matrix.png`** – heatmap of minimum distance (mm) for every PTV (rows) × critical OAR (columns: serial OARs plus eyes/lenses/hippocampi), colour-banded (red < 2, orange 2–5, yellow 5–10, light-green 10–20, green ≥ 20 mm). One glance shows which metastasis threatens which organ.
- **`nearest_critical_oar.png`** – per-PTV triage bar of the single nearest serial OAR (labelled with the organ and distance), sorted, with 3 mm / 5 mm threshold lines.
- **`sphericity_vs_elongation.png`** – scatter of shape character for anatomical structures (x = elongation, y = sphericity, marker size ∝ √volume, colour = category); round targets cluster top-left, elongated serial OARs (cord, optic nerves) sit far right.
- **`gtv_ptv_margin.png`** – per-lesion GTV vs PTV volume (log axis) with the implied isotropic margin (difference of the equivalent-sphere radii), a quick check that every metastasis received a consistent GTV→PTV expansion. The implied margin is exact only for spherical lesions. For a less round lesion it over-estimates the true margin by roughly a factor $1/\Psi$, because the volume gained per millimetre of margin equals the lesion's surface area, which exceeds that of the volume-equivalent sphere by $1/\Psi$.

### Structure classification

Before any plotting, each ROI is classified into **Target / serial-OAR / parallel-OAR / helper / external / marker** (`classify_structure`), so every plot shares one colour vocabulary and helper/union/optimisation structures (e.g. `h_PTV_gesamt`, dose shells, `opt system`) and POINT markers / the external body contour are kept out of the clinical plots and metrics. GTV/PTV pairs are matched by a lesion key derived from the ROI name.

## Case Transform Visualisation

When `case_modifier` applies a rigid transform to a CT + RTSTRUCT pair, it calls `visualizer.run_case_visualization` to produce a *before/after* comparison of the contour geometry. These plots need only the contour points of the original and transformed RTSTRUCT — no CT pixel data — so they are cheap to generate in both `resample` and `metadata` modes. The (expensive) CT body surface is opt-in via `--viz-ct-surface`.

### A. `transform_3d.html` – Interactive Before/After

**Chart type:** Plotly `Scatter3d` point clouds plus marker/axis overlays.

The original contour points (grey-blue) and transformed contour points (red) are drawn in the same DICOM patient coordinate system, so the displacement and any rotation are directly visible. Every layer is an **independently toggleable legend entry**, which is the natural place to switch the reference geometry on and off:

- **Konturen (Original)** / **Konturen (Transformiert)** — the two point clouds.
- **Drehpunkt / Rotationszentrum** + **Translationsvektor** — the chosen rotation centre and a dashed line to its post-transform position `c + t`.
- **POINT-Marker (Original)** / **POINT-Marker (Transformiert)** — every POINT-type ROI, with names, original and `T`-mapped (off by default).
- **DICOM-Achsen** — an L/P/S triad at the rotation centre (off by default).
- **CT-Koerperoberflaeche** — only present with `--viz-ct-surface`; hidden until toggled.

The scene uses `aspectmode='data'` so 1 mm is the same on-screen length in X, Y, and Z. Point clouds are deterministically sub-sampled (at most 600 points per ROI and 9000 per layer) so the file stays responsive even with dozens of ROIs.

### B. `transform_overview.png` – Static Tri-Planar Projection

**Chart type:** Three orthographic 2D scatter panels (axial X-Y, coronal X-Z, sagittal Y-Z).

A headless, report/CI-friendly companion to the HTML. Each panel overlays the original (grey) and transformed (red) contour points with equal aspect, marks the rotation centre (green star), draws the translation vector as an arrow, and shows the POINT markers. A pure translation appears as a uniform shift; a rotation shows up as a visibly tilted/rotated red cloud relative to the grey one.

### C. `displacement.png` – Per-ROI Centroid Displacement

**Chart type:** Horizontal bar chart, sorted descending, limited to the 40 most-displaced ROIs.

For each ROI the displacement $\lVert\mathbf{T}(\bar{\mathbf{p}}) - \bar{\mathbf{p}}\rVert$ of its vertex centroid $\bar{\mathbf{p}}$ is plotted, with a dashed reference line at the pure-translation magnitude $\lVert\mathbf{t}\rVert$. A point $\mathbf{p}$ is displaced by $(\mathbf{R} - \mathbf{I})(\mathbf{p} - \mathbf{c}) + \mathbf{t}$, where $\mathbf{c}$ is the rotation centre. The rotational part has magnitude $2 d_\perp \sin(\theta/2)$, with $\theta$ the net rotation angle and $d_\perp$ the distance of $\mathbf{p}$ from the rotation axis through $\mathbf{c}$. Structures on or near that axis therefore move by ≈ $\lVert\mathbf{t}\rVert$. Structures farther from it get an additional lever-arm displacement that can add to, or partly cancel, the translation, so their bars may lie on either side of the reference line. This makes it a quick QA check: an unexpectedly large displacement flags a structure that swings a long way under the requested rotation.

## Implementation Notes

### Structure Filtering and Auto-Detection

If `--targets` or `--oars` are not specified on the CLI, structures are classified automatically by `analyzer.classify_structure`, which combines the ROI name, the DICOM `RTROIInterpretedType`, and the contour geometry type into the categories Target / serial-OAR / parallel-OAR / helper / external / marker. Name-based rules deliberately take precedence over the DICOM tag for helper structures, because exported structure sets frequently mistag union PTVs as `GTV` and optimisation shells as `ORGAN`. POINT markers, the external body contour, and helper/union/optimisation structures are kept out of the clinical plots; they remain listed in `statistics.txt`.

### Scalability

The heatmap and lollipop layouts are chosen to stay readable at ~40 structures and beyond (grouped bar charts become unreadable above ~30). The distance lollipop is limited to the 25 most critical pairs regardless of total structure count; the full pairwise data remain in the proximity matrix and `statistics.txt`.

### Matplotlib Backend

`matplotlib.use("Agg")` is set at module level so the visualizer can run on headless servers (e.g., CI pipelines, remote compute nodes) without an X display. All output is written to files; no interactive window is opened.

---

# Dose Index Documentation

## Overview

`dose_indices` evaluates a dose distribution against one or more target volumes. It reads the RTSTRUCT (contours), the RTDOSE (3D dose grid) and, if present, the RTPLAN (prescription), rasterises the target and the isodose regions onto one common fine grid, and computes the conformity, gradient and homogeneity indices together with the DVH statistics of the target. The isodose regions and the helper regions used by the indices are written back as contours into a separate RTSTRUCT, so a planner can overlay them on the CT in the treatment planning system and see *exactly* which voxels the numbers were computed from. Everything is computed on the same grid with the same sampling, so the exported contours and the reported volumes are consistent by construction.

## Data Foundation: RTDOSE Grid

An RTDOSE object stores the dose as a multi-frame image: `pixel_array` has the shape *(frames, rows, columns)*, its integer values are converted to Gy with `DoseGridScaling`, and the frames are placed along the slice normal by the `GridFrameOffsetVector` (GFOV). The GFOV is either *relative* (first value 0, offsets from `ImagePositionPatient`) or *absolute* (z coordinates); the tool accepts both and requires the offsets to be equidistant. `DoseUnits` must be `GY`; a `DoseSummationType` other than `PLAN` (e.g., a single beam or fraction) is accepted with a warning, because the indices then describe a partial dose. The frame of reference must match the RTSTRUCT.

The voxel-to-patient mapping uses the same affine convention as the CT modules (see *Voxel-to-Patient Affine Matrix*), with the frame step from the GFOV:

```math
\mathbf{p} = \mathbf{A} \begin{pmatrix} k \\ j \\ i \\ 1 \end{pmatrix}, \qquad
\mathbf{A} = \begin{pmatrix} \Delta_k \hat{\mathbf{n}} & \Delta_r \hat{\mathbf{c}} & \Delta_c \hat{\mathbf{r}} & \mathbf{p}_0 \\ 0 & 0 & 0 & 1 \end{pmatrix}
```

where $k$ is the frame index, $\hat{\mathbf{n}}$ the slice normal, $\Delta_k$ the GFOV step and $\mathbf{p}_0$ the position of the first frame. The dose is never resampled onto a coarser grid; all evaluation points are obtained by interpolating the native grid.

**Prescription.** The reference dose $D_{\mathrm{Rx}}$ is taken from the RTPLAN `DoseReferenceSequence` (`TargetPrescriptionDose` of the `TARGET` reference, the total plan dose), unless `--rx` (absolute Gy) or `--rx-pct-of-max` (percentage of the grid maximum, the usual SRS convention) is given. **Target selection.** Without `--target`, all ROIs that the analyzer classifies as targets and whose name starts with `PTV` are evaluated; if the plan's `DoseReferenceDescription` (a 16-character DICOM SH string, i.e. usually a truncated ROI name) is a prefix of exactly one of them, only that one is used.

## Evaluation Grid

All masks live on **one** fine grid so that intersections and Dice coefficients are exact set operations:

- **In-plane resolution** `--grid` (default 0.25 mm). The axes are snapped to multiples of 1 mm so that the 1.0 / 0.5 / 0.25 / 0.1 mm grids nest; voxel centres lie at $x_0 + (i + \tfrac{1}{2}) \Delta$.
- **z axis** = the native dose planes, restricted to planes that also exist in the CT series (the exported contours must reference a CT slice). Contour planes must coincide with grid planes; if they do not, the z axis is refined by an integer factor (up to 8), otherwise the run stops with an error.
- **Bounding box** = union of the target contours and the native voxels reaching the lowest requested isodose level, plus a margin. If an isodose touches the grid boundary a warning is printed.
- **Dose interpolation** `--dose-interp`: trilinear (default) or cubic B-spline (the spline coefficients are computed once on a cropped copy of the native array, as in the modifier). Cubic interpolation follows a smooth dose field more faithfully than trilinear interpolation (see *Resolution and Convergence*) but overshoots near steep gradients; it tends to give a slightly larger prescription isodose volume and a slightly lower Paddick CI.

A 60 mm cube at 0.25 mm in-plane and 1 mm slices is 3.5 million voxels (14 MB of float32 dose per copy). Grids above 40 million voxels are refused with a hint to use a coarser `--grid`.

## Contour Rasterisation and Volume Models

Each target is rasterised with the analyzer's `rasterize_contours`: every contour ring is tested against the voxel centres of its plane and combined by XOR, so nested contours are holes and islands add up (the same even-odd rule as the planimetric volume). The volume is the weighted voxel sum

```math
V = \Delta z \sum_{k} w_k \, N_k \, \Delta^2
```

where $N_k$ is the number of voxels inside the structure on plane $k$ and $w_k$ the *slab weight* of that plane. Two conventions are available:

- **`slab`** (default, consistent with the analyzer): every contour plane represents a full slab, $w_k = 1$.
- **`eclipse`**: the first and the last plane of every contiguous run of contour planes count half, $w_k = 1/2$. This reproduces the volumes Eclipse reports in its DVH statistics; the difference to the planimetric volume is half the area of the first and the last contour times the slice spacing, which is noticeable for small targets on few slices.

The weights are applied to everything derived from the *target* (its volume, the intersection with the isodose, the underdosed part and every DVH statistic), never to the isodose regions themselves, so that $V_{\mathrm{TV} \cap \mathrm{PIV}} \le \min(\mathrm{TV}, \mathrm{PIV})$ always holds. Both models are always evaluated; the report shows the chosen one and a comparison line for the other. Isodose regions are simply the fine-grid voxels whose interpolated dose is $\ge$ the level.

**DVH statistics.** The target's dose values are the interpolated doses at its voxel centres, weighted with $w_k$. $D_x$ is the dose received by at least $x$ % of the weighted target volume (a weighted percentile with linear interpolation; $D_{98}$ is the 2nd percentile, $D_2$ the 98th, $D_{50}$ the median). $V_{95}$ and $V_{100}$ are the weighted fractions of the target receiving at least $0.95 D_{\mathrm{Rx}}$ and $D_{\mathrm{Rx}}$. $D_{\min}$ and $D_{\max}$ are the extreme interpolated samples and therefore depend on the grid resolution.

## Index Definitions

With $\mathrm{TV}$ the target volume, $\mathrm{PIV}$ the prescription isodose volume ($D \ge D_{\mathrm{Rx}}$), $\mathrm{TV}_{\mathrm{PIV}}$ their intersection and $\mathrm{PIV}_{50}$ the volume receiving at least half the prescription dose (all in cm³):

```math
\mathrm{CI}_{\mathrm{Paddick}} = \frac{\mathrm{TV}_{\mathrm{PIV}}^{\,2}}{\mathrm{TV} \cdot \mathrm{PIV}}
= \underbrace{\frac{\mathrm{TV}_{\mathrm{PIV}}}{\mathrm{TV}}}_{\text{coverage}} \cdot
  \underbrace{\frac{\mathrm{TV}_{\mathrm{PIV}}}{\mathrm{PIV}}}_{\text{selectivity}}
```

The Paddick index [Paddick 2000] is identical to the conformation number of van't Riet [van't Riet 1997]; it is 1 only for a perfectly conformal plan and penalises both under-coverage and spill. The RTOG index [Shaw 1993] $\mathrm{CI}_{\mathrm{RTOG}} = \mathrm{PIV}/\mathrm{TV}$ and the Dice coefficient $2 \mathrm{TV}_{\mathrm{PIV}} / (\mathrm{TV} + \mathrm{PIV})$ are reported alongside.

```math
\mathrm{GI} = \frac{\mathrm{PIV}_{50}}{\mathrm{PIV}}, \qquad
\mathrm{GM} = r_{\mathrm{eq}}(\mathrm{PIV}_{50}) - r_{\mathrm{eq}}(\mathrm{PIV}), \qquad
r_{\mathrm{eq}}(V) = \left( \frac{3V}{4\pi} \right)^{1/3}
```

The gradient index [Paddick & Lippitt 2006] measures the dose fall-off outside the target; the gradient measure GM (in cm) is the difference of the equivalent-sphere radii of the two isodose volumes, which is what Eclipse displays.

```math
\mathrm{HI}_{\mathrm{ICRU\,83}} = \frac{D_{2\%} - D_{98\%}}{D_{50\%}}
```

The ICRU 83 homogeneity index is 0 for a perfectly uniform target dose; SRS plans that prescribe to a low isodose line (e.g. 20 Gy at the 80 % line, so $D_{\max}$ = 25 Gy) have values around 0.2 by design. All quantities are `NaN` (written as `null` in the JSON) when a denominator is zero, e.g. when the prescription dose exceeds the grid maximum.

## PIV Scoping for Multi-Target Plans

In a plan with several targets in one dose grid the prescription isodose consists of several disconnected components. With `--piv-scope component` (default) the PIV of a target is the union of the components (6-connectivity) of the Rx isodose that overlap the target; the same rule is applied to the 50 % isodose for the gradient index. If no component overlaps the target the nearest one is used and flagged. The global PIV is always reported too, and with more than one target an additional block evaluates the union of all targets against the global isodoses. `--piv-scope global` uses the whole isodose for every target, which is the textbook definition and what a TPS reports.

## Eclipse Compatibility Mode

Inspecting contours exported by Eclipse reveals how the TPS stores its structures: for structures at *High* resolution every contour vertex has one coordinate exactly on the CT pixel-centre lattice $x_0 + k \cdot \Delta_{\mathrm{px}}$ ($\Delta_{\mathrm{px}}$ = CT pixel spacing) while the other coordinate is freely interpolated, and the median vertex spacing is one CT pixel; structures at *Default* resolution (Eclipse's own `Dose …[Gy]` isodose structures among them) use a lattice of two CT pixels offset by half a pixel. This is the signature of an isoline traced on a raster whose cells are centred on the CT pixels: the contour crosses the grid lines of that raster.

`--eclipse-compat high` (or `default`) therefore sets every parameter that influences the comparison with the TPS at once:

| Parameter | `high` | `default` |
|---|---|---|
| evaluation grid | CT pixel spacing, voxel centres on the CT pixel centres | 2 × pixel spacing, centres offset by half a pixel |
| volume model | `eclipse` (end slabs half) | `eclipse` |
| PIV scope | `global` | `global` |
| dose interpolation | linear | linear |
| isodose contours | isolines of the sampled dose field, no simplification | same |

On a clinical SRS case this mode reproduced Eclipse's DVH statistics of the target (coverage, $D_{98}$ / $D_{50}$ / $D_2$, HI, target volume) to within about 1 %. The Paddick index, however, came out several hundredths above the TPS value; the difference lay entirely in the prescription isodose volume, i.e. in how the TPS samples the isodose, not in the target statistics. The report prints every component so the comparison can be completed once the TPS values are known.

## RTSTRUCT Export

The exported file is built from scratch (not a modified copy): patient, study and frame-of-reference attributes are copied from the original RTSTRUCT so the file imports next to the original CT series, everything else — UIDs, series number (+1000), `StructureSetLabel` (original label plus suffix, cut to 16 characters), dates, a summary of the indices in `StructureSetDescription` — is new. Each ROI is written as `CLOSED_PLANAR` contours, one per ring, on the CT slice planes with a `ContourImageSequence` reference to that slice; holes are separate rings with negative signed area, which is the convention Eclipse itself uses. The isodose ROIs carry the `RTROIInterpretedType` `CONTROL` like Eclipse's own dose structures. The file is written as explicit VR little endian with a complete file-meta header under pydicom's strict validation mode and immediately re-read and checked (`verify_rtstruct`: UID consistency, VR lengths, unique names, planar contours, slice references). Contours are extracted from the masks by marching squares at the 0.5 level, so their vertices lie on the voxel boundaries; the Douglas–Peucker simplification (`--simplify-mm`, default 0.1 mm) stays below half a voxel so that re-rasterising the contours reproduces the mask exactly (the self-test checks Dice = 1.0). With `--iso-contours field` the isodose contours are instead the isolines of the sampled dose field, which gives sub-voxel accuracy and the Eclipse vertex convention.

## Resolution and Convergence

The analytic phantom of the self-test (see *Validation*), evaluated through the same code path as a real case: a spherical target of radius $R$ contoured on 1 mm planes inside a smooth radial dose field whose prescription isodose (radius 10.5 mm resp. 5.8 mm) is shifted by 2 mm against the target, so TV, PIV, their intersection and the DVH have closed forms. The field is sampled from a 1 mm native dose grid like a TPS export; slab volume model, prescription 20 Gy. The rows marked *analytic* are the closed-form values; the second phantom is the same construction with the radius of a small SRS target:

| Target | Grid (in-plane) | Interpolation | TV cm³ | PIV cm³ | TV∩PIV cm³ | CI Paddick | GI | HI |
|---|---|---|---|---|---|---|---|---|
| R = 10 mm | 1.0 mm | linear | 4.224 | 4.776 | 3.820 | 0.723 | 2.46 | 0.298 |
| R = 10 mm | 0.5 mm | linear | 4.210 | 4.806 | 3.822 | 0.722 | 2.47 | 0.299 |
| R = 10 mm | 0.25 mm | linear | 4.190 | 4.801 | 3.808 | 0.721 | 2.48 | 0.296 |
| R = 10 mm | 0.25 mm | cubic | 4.190 | 4.843 | 3.823 | 0.720 | 2.45 | 0.294 |
| R = 10 mm | analytic | — | 4.194 | 4.838 | 3.824 | 0.721 | 2.46 | 0.294 |
| R = 5.5 mm | 1.0 mm | linear | 0.672 | 0.840 | 0.548 | 0.532 | 2.39 | 0.457 |
| R = 5.5 mm | 0.5 mm | linear | 0.680 | 0.792 | 0.538 | 0.537 | 2.54 | 0.462 |
| R = 5.5 mm | 0.25 mm | linear | 0.690 | 0.802 | 0.544 | 0.535 | 2.50 | 0.471 |
| R = 5.5 mm | 0.25 mm | cubic | 0.690 | 0.816 | 0.549 | 0.536 | 2.46 | 0.467 |
| R = 5.5 mm | analytic | — | 0.691 | 0.819 | 0.551 | 0.536 | 2.45 | 0.466 |

At 1 mm in-plane the voxel discretisation alone moves the individual volumes by about 1 % for the 10 mm sphere and by about 3 % for the small target (TV −2.7 %, PIV +2.6 %), and the gradient index of the small target by 0.06. At 0.25 mm both target volumes are within 0.1 % of the closed form; the remaining PIV gap of the linear rows (−0.8 % resp. −2 %) is the interpolation error of the 1 mm dose grid, which the cubic model removes for this smooth field. The Paddick index is more forgiving — the errors of TV, PIV and their intersection partly cancel, so it stays within 0.005 of the closed form on every row — but the comparison with a TPS is made on the components, not on the ratio, and only the fine grid makes them trustworthy individually, which is why 0.25 mm is the default. The remaining spread between interpolation schemes and volume models is not a numerical error but the genuine ambiguity of a sub-millimetre index; a report that hides it is not more accurate, only less honest.

## Validation

`python -m dicom_file_modifier.dose_indices --self-test` builds an analytic phantom in memory — a spherical target of radius 10 mm contoured on 1 mm planes and a smooth radial dose field $D(r) = D_{\max} / (1 + (r/r_0)^6)$ around a centre offset by 2 mm — for which TV, PIV, their intersection (a stack of circle–circle lenses), the DVH and hence every index have closed-form values. The 23 checks cover the RTDOSE loader (relative and absolute GFOV), the sampler (a linear field must be reproduced to 1e-6 Gy by both interpolation orders; points outside the grid give NaN), all volumes and indices at 0.25 mm (within 0.5 % / ±0.005) and 1.0 mm, both volume models, the sampling from a native 0.5 mm grid, a ring-shaped target with holes, the component scoping with two separate hot spots, and the mask → contour → mask round trip. `python -m dicom_file_modifier.rtstruct_writer --self-test` writes synthetic masks (sphere, annulus, islands, border-touching block, single voxel) to a temporary RTSTRUCT and verifies orientation, volumes, name uniqueness and the file-meta consistency.

## Limitations

- **The dose grid is the limit.** A 1 mm dose grid cannot locate an isodose surface better than its interpolation model allows; the fine evaluation grid removes the discretisation error of the volumes, not the uncertainty of the dose itself. Cubic interpolation is offered as a second model, not as the truth.
- **TPS agreement is empirical.** The volume model and the Eclipse-compatibility mode were tuned to reproduce the target statistics of one TPS (Eclipse) on clinical data; the prescription isodose volume of the TPS could not be reproduced from the exported data (see above). Index values should be compared with the TPS on each installation before they are used to judge plans.
- **Contour planes must coincide with dose planes** (up to an integer refinement); non-uniform GFOVs, relative dose units and tilted dose grids are rejected.
- **$D_{\min}$ / $D_{\max}$ are sample extremes** on the evaluation grid and therefore resolution dependent; the percentile values are robust.
- **Component scoping is a heuristic** for multi-target plans; targets that share one isodose component get the same PIV and are flagged.
- **Not a clinical tool.** See the disclaimer at the top: nothing here replaces the plan evaluation in the TPS by a medical physicist.

---

# References

Standards and reports:

- **ICRU Report 50** (1993). *Prescribing, Recording, and Reporting Photon Beam Therapy.*
- **ICRU Report 62** (1999). *Prescribing, Recording and Reporting Photon Beam Therapy (Supplement to ICRU Report 50).*
- **ICRU Report 83** (2010). *Prescribing, Recording, and Reporting Photon-Beam Intensity-Modulated Radiation Therapy (IMRT).*
- **DICOM Standard**, Part 3, Section C.8.8.5 – *Structure Set Module* (and C.8.8.6 – *ROI Contour Module*). [dicomstandard.org](https://www.dicomstandard.org/)
- **DICOM Standard**, Part 3, Sections C.7.6.2, C.7.6.3 – *Image Plane Module, Image Pixel Module.* [dicomstandard.org](https://www.dicomstandard.org/)

Publications:

- Feuvret, L., Noël, G., Mazeron, J.-J., & Bey, P. (2006). *Conformity index: a review.* International Journal of Radiation Oncology, Biology, Physics, 64(2), 333–342.
- Heimann, T., van Ginneken, B., Styner, M. A., et al. (2009). *Comparison and evaluation of methods for liver segmentation from CT datasets.* IEEE Transactions on Medical Imaging, 28(8), 1251–1265.
- Huttenlocher, D. P., Klanderman, G. A., & Rucklidge, W. J. (1993). *Comparing images using the Hausdorff distance.* IEEE Transactions on Pattern Analysis and Machine Intelligence, 15(9), 850–863.
- Lehmann, T. M., Gönner, C., & Spitzer, K. (1999). *Survey: Interpolation methods in medical image processing.* IEEE Transactions on Medical Imaging, 18(11), 1049–1075.
- Lorensen, W. E., & Cline, H. E. (1987). *Marching cubes: A high resolution 3D surface construction algorithm.* ACM SIGGRAPH Computer Graphics, 21(4), 163–169.
- Paddick, I. (2000). *A simple scoring ratio to index the conformity of radiosurgical treatment plans.* Journal of Neurosurgery, 93(Suppl 3), 219–222.
- Paddick, I., & Lippitt, B. (2006). *A simple dose gradient measurement tool to complement the conformity index.* Journal of Neurosurgery, 105(Suppl), 194–201.
- Shaw, E., Kline, R., Gillin, M., et al. (1993). *Radiation Therapy Oncology Group: radiosurgery quality assurance guidelines.* International Journal of Radiation Oncology, Biology, Physics, 27(5), 1231–1239.
- Taha, A. A., & Hanbury, A. (2015). *Metrics for evaluating 3D medical image segmentation: analysis, selection, and tool.* BMC Medical Imaging, 15(1), 29.
- Thévenaz, P., Blu, T., & Unser, M. (2000). *Interpolation revisited.* IEEE Transactions on Medical Imaging, 19(7), 739–758.
- van't Riet, A., Mak, A. C. A., Moerland, M. A., Elders, L. H., & van der Zee, W. (1997). *A conformation number to quantify the degree of conformality in brachytherapy and external beam irradiation: application to the prostate.* International Journal of Radiation Oncology, Biology, Physics, 37(3), 731–736.
- Wadell, H. (1935). *Volume, shape, and roundness of quartz particles.* The Journal of Geology, 43(3), 250–280.




