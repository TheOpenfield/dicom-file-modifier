# DICOM File Modifier

A comprehensive toolkit for analyzing, modifying, and visualizing DICOM RT Structure Set files used in radiotherapy planning. This project provides tools to process target volumes (PTV, CTV, GTV) and organs at risk (OAR) from DICOM files, compute geometric metrics, and generate visualizations.

## Features

- **Analyzer**: Extract and compute geometric properties (volume, centroid, shape metrics, distances) from RTSTRUCT files
- **Modifier**: Rigid body transformation of CT DICOM series (translation + rotation) with HU preservation and 3D visualization
- **Visualizer**: Generate plots and statistics from analysis results

## Project Structure

```
dicom-file-modifier/
├── data/                    # Input DICOM files (not synced)
├── output/                  # Analysis results and modified files (not synced)
├── dicom_file_modifier/     # Python package
│   ├── __init__.py
│   ├── analyzer.py          # RTSTRUCT analysis module
│   ├── modifier.py          # CT rigid body transformer
│   └── visualizer.py        # Visualization module
├── requirements.txt         # Python dependencies
└── README.md               # This file
```

## Installation

1. Clone the repository
2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

### Analyzer

Analyze RTSTRUCT files and compute geometric metrics:

```bash
python -m dicom_file_modifier.analyzer data/rtstruct.dcm --output output/
```

This will generate:
- `analysis_results.json`: Computed metrics for all structures
- Console output with summary statistics

### Visualizer

Generate plots and statistics directly from a RTSTRUCT file:

```bash
python -m dicom_file_modifier.visualizer data/0000000171/RS.dcm --output output/
# with explicit structure selection:
python -m dicom_file_modifier.visualizer data/0000000171/RS.dcm \
    --targets GTV,PTV --oars Hirnstamm,Rueckenmark --output output/
```

This creates:
- `volumes.png`: Horizontal bar chart of all structure volumes, colour-coded by type
- `shape_metrics.png`: Grouped bar chart of sphericity, compactness, elongation per structure
- `distances.png`: Grouped bar chart of min/Hausdorff/centroid distances for the 25 most critical Target–OAR pairs
- `centroids_3d.png`: 3D scatter of structure centroids in patient space, marker size ∝ volume
- `statistics.txt`: Full numerical summary of all metrics

### Modifier – CT Rigid Body Transformer

Apply a rigid body transformation (translation + rotation) to a CT DICOM series:

```bash
python -m dicom_file_modifier.modifier data/0000000171/CT \
    --tx 10 --ty 0 --tz -5 \
    --rx 0  --ry 0 --rz 15 \
    --output output/ct_transformed
```

This produces:
- `CT_0000.dcm … CT_NNNN.dcm`: Transformed CT series, importable into any TPS
- `visualization_3d.html`: Interactive 3D comparison (original vs. transformed)

**All CLI options:**

| Option | Default | Description |
|---|---|---|
| `--tx/ty/tz` | 0 mm | Translation along X / Y / Z axis |
| `--rx/ry/rz` | 0 ° | Rotation around X / Y / Z axis (extrinsic) |
| `--method` | `resample` | `resample`: standard axial output; `metadata`: exact HU preservation |
| `--order` | `1` | Interpolation order: 0 = nearest neighbour, 1 = linear, 3 = cubic |
| `--output` | `output/ct_transformed` | Output directory |
| `--save-viz` | *(auto)* | Path for HTML visualization file |
| `--no-viz` | off | Skip visualization |

## Dependencies

- pydicom: DICOM file handling
- numpy: Numerical computations
- scipy: Scientific computing (distances, convex hull, interpolation, rotations)
- shapely: 2D geometry operations
- matplotlib: Plotting and visualization
- scikit-image: Marching cubes surface extraction (modifier)
- plotly: Interactive 3D visualization (modifier)

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

- **GTV** (*Gross Tumor Volume*): Macroscopically visible tumor tissue.
- **CTV** (*Clinical Target Volume*): GTV plus safety margin for microscopic spread.
- **PTV** (*Planning Target Volume*): CTV plus safety margin for setup uncertainties and organ motion.

Organs at risk (OAR) are structures whose radiation dose must be limited, e.g., spinal cord, parotid, or bladder. Classification is stored in the `RTROIInterpretedType` field.

## Calculation Methods

### Volume Calculation

The volume of a structure is calculated slice-by-slice by summing contour areas:

$$V = \sum_{i=1}^{N} A_i \cdot \Delta z$$

Where *Aᵢ* is the area of the contour in slice *i* and *Δz* is the mean slice spacing. The area of each contour is determined using the **Shoelace formula** (Gauss's trapezoidal rule):

$$A = \frac{1}{2} \left| \sum_{j=0}^{n-1} (x_j \, y_{j+1} - x_{j+1} \, y_j) \right|$$

If a slice contains multiple contours (e.g., for ring-shaped or fragmented structures), their areas are added.

**Unit:** Results are converted from mm³ to cm³ (division by 1000).

**Limitation:** This method is a planimetric approximation. It is more accurate with thinner CT slices. For exact voxel-based volume calculation, the associated CT image would be required as a reference grid.

### Centroid Calculation

The centroid is calculated as the **area-weighted average** of contour centroids:

$$\vec{C} = \frac{\sum_{i=1}^{N} A_i \cdot \vec{c}_i}{\sum_{i=1}^{N} A_i}$$

Where *c⃗ᵢ* is the geometric center of the contour in slice *i* (arithmetic mean of all contour points) and *Aᵢ* is the corresponding contour area as weight.

Weighting by area ensures that larger cross-sections contribute more to the overall centroid than small edge slices. The result is three coordinates *(x, y, z)* in millimeters in the DICOM patient coordinate system.

### Shape Metrics

Shape metrics quantify the geometric shape of a structure independently of its absolute size.

#### Sphericity

Sphericity describes how closely a structure resembles a sphere. It is defined as the ratio of the surface area of a volume-equivalent sphere to the actual surface area:

$$\Psi = \frac{\pi^{1/3} \cdot (6V)^{2/3}}{A_{\text{surface}}}$$

A value of 1.0 corresponds to a perfect sphere; smaller values indicate irregular or elongated shapes. Since the true surface area is not trivial to determine from contour data, the **surface area of the convex hull** of the point cloud is used as an approximation (calculated via `scipy.spatial.ConvexHull`).

#### Compactness

Compactness relates the actual volume to the volume of the convex hull:

$$K = \frac{V}{V_{\text{convex}}}$$

A value close to 1.0 means the structure has few indentations or cavities. Low values indicate concave or highly irregular shapes, which may be clinically relevant for tumors that wrap around other structures.

#### Elongation

Elongation describes the stretching of a structure along its principal axes. It is determined via **principal component analysis** (PCA) of the 3D point cloud:

$$E = \sqrt{\frac{\lambda_{\max}}{\lambda_{\min}}}$$

Where *λ_max* and *λ_min* are the largest and smallest eigenvalues of the covariance matrix. A value of 1.0 corresponds to isotropic (spherical) extension, larger values show increasing stretching in a preferred direction.

#### Bounding Box

The axis-aligned bounding box gives the extent of the structure in all three spatial directions:

$$\Delta x = x_{\max} - x_{\min}, \quad \Delta y = y_{\max} - y_{\min}, \quad \Delta z = z_{\max} - z_{\min}$$

It provides a quick overview of the spatial extent in millimeters.

## Distance Calculations

Distances between structures are central to evaluating the radiotherapy plan: How close is an organ at risk to the target? Do structures overlap?

### Minimum Distance

The minimum distance between two structures *A* and *B* is defined as:

$$d_{\min}(A, B) = \min_{a \in A, \, b \in B} \| a - b \|_2$$

Calculation is performed efficiently using a **KD-tree** (`scipy.spatial.cKDTree`): For each point of structure *A*, the nearest neighbor in *B* is found and the global minimum determined.

A value of 0 mm means the contours touch or overlap. Clinically, this value is particularly relevant for assessing whether a safety margin between PTV and adjacent organs at risk is maintained.

### Hausdorff Distance

The Hausdorff distance is a measure of the maximum deviation between two point sets:

$$d_H(A, B) = \max\!\Big(\,\sup_{a \in A} \inf_{b \in B} \|a - b\|, \;\sup_{b \in B} \inf_{a \in A} \|a - b\|\,\Big)$$

Intuitively: The Hausdorff distance indicates how far one must travel in the worst case from a point of one structure to the nearest point of the other structure. It is thus more sensitive to local outliers than the minimum distance and is well-suited for assessing shape agreement between two structures.

### Centroid Distance

The Euclidean distance between the centroids of two structures:

$$d_C = \| \vec{C}_A - \vec{C}_B \|_2$$

This value provides a rough but robust estimate of spatial separation. It is insensitive to outliers and suitable as a quick comparison value.

## Performance Notes

Large structures can comprise tens of thousands of contour points. To perform distance calculations in acceptable time, point clouds are reduced to a configurable sample size if necessary (default: 5,000 points for minimum distance, 3,000 for Hausdorff distance). This results in a small inaccuracy that is negligible for clinical practice.

## Limitations

- **No voxel-based analysis:** The script works exclusively with contour points from the RTSTRUCT file. For overlap metrics like the **Dice coefficient** or **Conformity Number**, the associated CT would be required as a reference grid to rasterize the contours into a 3D voxel grid.
- **Surface approximation:** Sphericity uses the convex hull as an approximation of the actual surface. For structures with strong indentations, the surface is underestimated and sphericity correspondingly overestimated.
- **Planimetric volume:** Volume calculation assumes equidistant slice spacing. With non-equidistant slices, the mean spacing is used, which can lead to small inaccuracies.
- **Not a clinical diagnostic tool:** The script serves geometric analysis and does not replace clinical evaluation by a medical physicist or radiation therapist.

## Used Libraries

| Library | Version | Purpose |
|---|---|---|
| `pydicom` | ≥ 2.3 | Reading DICOM files |
| `numpy` | ≥ 1.21 | Numerical calculations and linear algebra |
| `scipy` | ≥ 1.7 | KD-tree, Hausdorff distance, convex hull, interpolation, rotations |
| `shapely` | ≥ 1.8 | 2D polygon operations (area calculation) |
| `matplotlib` | ≥ 3.5 | Static plotting and visualization |
| `scikit-image` | ≥ 0.19 | Marching cubes surface extraction (modifier) |
| `plotly` | ≥ 5.0 | Interactive 3D visualization (modifier) |

## Literature

1. **ICRU Report 50** (1993). *Prescribing, Recording, and Reporting Photon Beam Therapy.*
2. **ICRU Report 62** (1999). *Prescribing, Recording and Reporting Photon Beam Therapy (Supplement to ICRU Report 50).*
3. **ICRU Report 83** (2010). *Prescribing, Recording, and Reporting Photon-Beam Intensity-Modulated Radiation Therapy (IMRT).*
4. **DICOM Standard**, Part 3, Section C.8.8.6 – *RT Structure Set Module.* [dicom.nema.org](https://www.dicomstandard.org/)
5. Huttenlocher, D. P., Klanderman, G. A., & Rucklidge, W. J. (1993). *Comparing images using the Hausdorff distance.* IEEE Transactions on Pattern Analysis and Machine Intelligence, 15(9), 850–863.
6. **DICOM Standard**, Part 3, Sections C.7.6.2, C.7.6.3 – *Image Plane Module, Image Pixel Module.*
7. Lehmann, T. M., Gönner, C., & Spitzer, K. (1999). *Survey: Interpolation methods in medical image processing.* IEEE Transactions on Medical Imaging, 18(11), 1049–1075.
8. Lorensen, W. E., & Cline, H. E. (1987). *Marching cubes: A high resolution 3D surface construction algorithm.* ACM SIGGRAPH Computer Graphics, 21(4), 163–169.

---

# CT Rigid Body Transformer Documentation

## Overview

The **CT Rigid Body Transformer** (`modifier.py`) applies a rigid body transformation — consisting of a translation in three spatial directions and a rotation around three spatial axes — to a CT DICOM series. The result is a new DICOM series that can be imported into any treatment planning system (TPS). The central design goal is to preserve the original Hounsfield Unit (HU) values as accurately as possible while guaranteeing that no geometric distortion of the image is introduced.

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
| `(0020,0032)` | `ImagePositionPatient` (IPP) | 3D position of the first pixel (row 0, col 0) in mm |
| `(0028,0030)` | `PixelSpacing` | In-plane pixel size [row spacing, col spacing] in mm |
| `(0028,1053)` | `RescaleSlope` | Linear HU conversion: HU = stored × slope + intercept |
| `(0028,1052)` | `RescaleIntercept` | See above |
| `(0028,0103)` | `PixelRepresentation` | 0 = unsigned, 1 = signed integer |

### Voxel-to-Patient Affine Matrix

The spatial position of any voxel $(k, j, i)$ — where $k$ is the slice index, $j$ the row index, and $i$ the column index — in patient coordinates is given by the following affine transformation:

$$\begin{pmatrix} x \\ y \\ z \\ 1 \end{pmatrix} = \mathbf{A} \begin{pmatrix} k \\ j \\ i \\ 1 \end{pmatrix}$$

where the $4 \times 4$ affine matrix $\mathbf{A}$ is constructed from the DICOM tags as:

$$\mathbf{A} = \begin{pmatrix} n_x \cdot \Delta z & F_4 \cdot \Delta r & F_1 \cdot \Delta c & \text{IPP}_x \\ n_y \cdot \Delta z & F_5 \cdot \Delta r & F_2 \cdot \Delta c & \text{IPP}_y \\ n_z \cdot \Delta z & F_6 \cdot \Delta r & F_3 \cdot \Delta c & \text{IPP}_z \\ 0 & 0 & 0 & 1 \end{pmatrix}$$

Here $\mathbf{F} = (F_1, F_2, F_3, F_4, F_5, F_6)$ is the `ImageOrientationPatient` vector, $(F_1, F_2, F_3)$ are the direction cosines of the row direction (increasing column index) and $(F_4, F_5, F_6)$ are the direction cosines of the column direction (increasing row index). The slice normal $\mathbf{n} = (n_x, n_y, n_z) = (F_1, F_2, F_3) \times (F_4, F_5, F_6)$ is computed as the cross product of the two IOP vectors. $\Delta z$ is the slice spacing, $\Delta r$ the row pixel spacing, and $\Delta c$ the column pixel spacing.

The inverse $\mathbf{A}^{-1}$ maps patient coordinates back to voxel indices and is used during resampling.

### HU Conversion

Stored integer pixel values are converted to Hounsfield Units via a linear mapping defined per slice:

$$\text{HU} = \text{stored} \times \text{RescaleSlope} + \text{RescaleIntercept}$$

For modern CT scanners the slope is typically 1 and the intercept −1024, placing air at −1024 HU and water at 0 HU. The full diagnostic CT range spans approximately −1024 HU (air) to +3071 HU (dense bone / metal).

## Rigid Body Transformation

### Definition

A rigid body transformation in three-dimensional Euclidean space preserves all pairwise distances and angles. It comprises exactly six degrees of freedom: three translational $(t_x, t_y, t_z)$ and three rotational $(r_x, r_y, r_z)$. Formally it is an element of the special Euclidean group $SE(3)$:

$$\mathbf{T} : \mathbf{p} \mapsto \mathbf{R}\,\mathbf{p} + \mathbf{t}$$

where $\mathbf{R} \in SO(3)$ is a $3 \times 3$ rotation matrix satisfying $\mathbf{R}^T \mathbf{R} = \mathbf{I}$ and $\det(\mathbf{R}) = +1$, and $\mathbf{t} \in \mathbb{R}^3$ is the translation vector. In homogeneous coordinates this becomes the $4 \times 4$ matrix:

$$\mathbf{T} = \begin{pmatrix} \mathbf{R} & \mathbf{t} \\ \mathbf{0}^T & 1 \end{pmatrix}$$

### Rotation Matrix Construction

The rotation matrix is constructed from three rotation angles using **extrinsic Euler angles in XYZ order**: the patient is first rotated by $r_x$ around the fixed X axis, then by $r_y$ around the fixed Y axis, and finally by $r_z$ around the fixed Z axis. The combined rotation matrix is:

$$\mathbf{R} = \mathbf{R}_z(r_z)\,\mathbf{R}_y(r_y)\,\mathbf{R}_x(r_x)$$

with the elementary rotation matrices:

$$\mathbf{R}_x(\alpha) = \begin{pmatrix} 1 & 0 & 0 \\ 0 & \cos\alpha & -\sin\alpha \\ 0 & \sin\alpha & \cos\alpha \end{pmatrix}, \quad \mathbf{R}_y(\beta) = \begin{pmatrix} \cos\beta & 0 & \sin\beta \\ 0 & 1 & 0 \\ -\sin\beta & 0 & \cos\beta \end{pmatrix}, \quad \mathbf{R}_z(\gamma) = \begin{pmatrix} \cos\gamma & -\sin\gamma & 0 \\ \sin\gamma & \cos\gamma & 0 \\ 0 & 0 & 1 \end{pmatrix}$$

Extrinsic rotation around fixed world axes is used in preference to intrinsic (body-fixed) rotation because it is more intuitive in the clinical setting: $r_z$ always corresponds to a rotation in the axial plane regardless of the other angles applied.

The implementation uses `scipy.spatial.transform.Rotation.from_euler("XYZ", ...)` where uppercase letters indicate extrinsic convention.

### Rotation Centre

The rotation is performed about the **geometric centre** of the CT volume in patient coordinates:

$$\mathbf{c} = \mathbf{A}\, \begin{pmatrix} (N_z-1)/2 \\ (N_y-1)/2 \\ (N_x-1)/2 \\ 1 \end{pmatrix}$$

Rotating around the volume centre ensures that the patient body remains approximately centred within the output voxel grid and does not drift outside the field of view for small angles. The full forward transformation applied to a point $\mathbf{p}$ is therefore:

$$\mathbf{p}' = \mathbf{R}\,(\mathbf{p} - \mathbf{c}) + \mathbf{c} + \mathbf{t} = \mathbf{R}\,\mathbf{p} + \underbrace{(-\mathbf{R}\,\mathbf{c} + \mathbf{c} + \mathbf{t})}_{\mathbf{t}_{\text{eff}}}$$

which is stored compactly in the upper-right column of the $4 \times 4$ matrix $\mathbf{T}$.

## Resampling Method (`--method resample`)

### Principle: Inverse Mapping

When a rigid body transformation moves a patient, a new CT scan of the patient in their new position is simulated. The standard approach in medical image processing is **inverse mapping** (also called pull-back or backward mapping): instead of pushing each source voxel into the output grid — which would leave holes — each output voxel asks "where in the source volume did this intensity come from?"

For each output voxel at voxel index $(k, j, i)$, the corresponding patient position is $\mathbf{p}_{\text{out}} = \mathbf{A}\,[k,j,i,1]^T$. The source position in the original (untransformed) volume is:

$$\mathbf{p}_{\text{in}} = \mathbf{T}^{-1}\,\mathbf{p}_{\text{out}}$$

Converting back to voxel indices: $\mathbf{v}_{\text{in}} = \mathbf{A}^{-1}\,\mathbf{p}_{\text{in}}$. Combining these steps, the complete voxel-to-voxel mapping is:

$$\mathbf{v}_{\text{in}} = \underbrace{\mathbf{A}^{-1} \mathbf{T}^{-1} \mathbf{A}}_{\mathbf{M}}\, \mathbf{v}_{\text{out}}$$

The matrix $\mathbf{M}$ is computed once before the loop and applied to all voxel coordinates simultaneously. The source coordinates $\mathbf{v}_{\text{in}}$ are generally non-integer; the HU value is obtained by **interpolation** in the source volume.

Voxels whose source coordinates fall outside the original volume boundaries are assigned −1000 HU (air), which corresponds to the physical situation where the patient body has moved out of the detector field.

### Why Inverse Rather Than Forward Mapping?

Forward mapping (pushing each source voxel to its new position) suffers from two problems: (1) output voxels that no source voxel maps to remain unfilled (holes), and (2) output voxels that multiple source voxels map to require a compositing strategy. Inverse mapping avoids both problems by construction and is the standard in medical image registration (e.g. SPM, FSL, ITK).

### Memory-Efficient Chunk Processing

For a typical CT volume of $512 \times 512 \times 320$ voxels, computing the full coordinate array at once would require approximately 2.5 GB of RAM. The implementation therefore processes the volume in **chunks of 20 slices at a time**, generating the coordinate array only for the current chunk. Peak memory usage per chunk is approximately 80 MB for the coordinate arrays plus the volume data itself.

### Interpolation

The source-volume lookup at non-integer coordinates requires interpolation. Three orders are available:

#### Nearest Neighbour (order 0)

$$\text{HU}(\mathbf{v}) = \text{HU}\bigl(\text{round}(\mathbf{v})\bigr)$$

No weighted averaging; the nearest voxel value is taken directly. This guarantees that **only HU values that actually exist in the original volume appear in the output**. The disadvantage is staircase artefacts at structure boundaries (Gibbs phenomenon). Suitable when exact discrete HU preservation is strictly required (e.g., for lookup-table-based TPS dose calculations).

#### Trilinear Interpolation (order 1, default)

For a point at fractional position $(k + \delta_k,\, j + \delta_j,\, i + \delta_i)$ with $\delta \in [0,1)$, the value is the weighted average of the eight surrounding voxels:

$$\text{HU}(\mathbf{v}) = \sum_{a \in \{0,1\}} \sum_{b \in \{0,1\}} \sum_{c \in \{0,1\}} w_{abc}\,\text{HU}(k{+}a,\, j{+}b,\, i{+}c)$$

with trilinear weights $w_{abc} = |1-\delta_k-a|\cdot|1-\delta_j-b|\cdot|1-\delta_i-c|$. Trilinear interpolation is continuous, preserves the value range of the original data exactly (no overshoot), and introduces a maximum HU error of approximately half the local gradient magnitude. For soft tissue with smooth HU gradients, deviations are typically **< 2 HU** — well within clinical relevance thresholds and scanner reproducibility (5–10 HU). Trilinear interpolation is the clinical standard in image registration (e.g. in Eclipse, RayStation).

#### Tricubic B-Spline Interpolation (order 3)

A continuous piecewise cubic approximation, computed via a separable recursive spline filter applied once to the entire volume before sampling. Produces the smoothest output and lowest systematic HU error (typically **< 0.5 HU** for smooth tissue regions), but introduces slight overshoot near sharp boundaries (e.g. bone–air interfaces), which can generate HU values marginally outside the original range. Recommended when maximum interpolation quality is required and the slight computational overhead is acceptable.

### HU Value Preservation: Verification

The implementation verifies HU preservation by comparing the minimum and maximum HU values of the original and transformed volumes:

```
HU range original:     [-1024, 3071] HU
HU range transformed:  [-1024, 3071] HU
```

For trilinear (order 1) interpolation the range is guaranteed to be preserved exactly (no overshoot). For cubic (order 3) interpolation marginal overshoot at high-contrast boundaries is theoretically possible but bounded and clinically irrelevant. Nearest-neighbour (order 0) interpolation is lossless by definition.

Note that a voxel-wise difference map between original and transformed volumes is **not a meaningful quality metric** in this context: after a rigid body motion the same tissue appears at different voxel positions in the two volumes, so direct subtraction compares different anatomical structures. A correct interpolation quality assessment requires a round-trip test (apply T then T⁻¹) or comparison within a co-registered reference frame.

### Output DICOM Structure

The resampled volume is written back to DICOM using the metadata of the original slices. The output geometry (IPP, IOP, pixel spacing, slice spacing) is **identical** to the input, meaning the output slices occupy the same spatial positions as the original. The transformed patient body appears shifted/rotated *within* the fixed voxel grid. Regions that the body moved into contain the resampled HU values; regions it moved out of contain −1000 HU (air).

The HU-to-stored conversion is the exact inverse of the loading step:

$$\text{stored} = \text{round}\!\left(\frac{\text{HU} - \text{RescaleIntercept}}{\text{RescaleSlope}}\right)$$

Values are clipped to the valid int16 range $[-32768, 32767]$ and stored as signed 16-bit integers (`PixelRepresentation = 1`), which is standard for CT. The original `RescaleSlope` and `RescaleIntercept` are preserved unchanged so the output is correctly calibrated in any TPS.

Each output file receives a new `SOPInstanceUID` and `SeriesInstanceUID` (generated with `pydicom.uid.generate_uid()` which produces DICOM-conformant UID strings) so that the transformed series is recognised as an independent series by the TPS and PACS, while the `StudyInstanceUID` and patient demographics remain identical for correct study association.

## Metadata-Only Method (`--method metadata`)

### Principle

For the metadata-only approach, **no pixel data is modified**. Instead, the spatial meaning of each slice is updated by transforming the DICOM positional tags:

**ImagePositionPatient (IPP):** The position of the first pixel of each slice is mapped through the forward transformation:

$$\text{IPP}'_k = \mathbf{T}\,\begin{pmatrix}\text{IPP}_k \\ 1\end{pmatrix}$$

**ImageOrientationPatient (IOP):** The six direction cosines (row direction and column direction) are rotated by $\mathbf{R}$:

$$\mathbf{F}'_{\text{row}} = \mathbf{R}\,\mathbf{F}_{\text{row}}, \qquad \mathbf{F}'_{\text{col}} = \mathbf{R}\,\mathbf{F}_{\text{col}}$$

Since $\mathbf{R} \in SO(3)$ preserves norms, the transformed direction cosines remain unit vectors and their cross product continues to define the correct slice normal.

### Advantages and Limitations

The metadata-only method guarantees **exact HU preservation** because no interpolation is performed. It is also significantly faster since no resampling loop is needed. However, after an arbitrary rotation the IOP vectors are no longer aligned with the standard axial orientation `[1,0,0,0,1,0]`. Some treatment planning systems (particularly older versions) require strictly axial CT imports and will reject datasets with oblique IOP. Modern systems such as RayStation and Eclipse generally handle oblique orientations correctly. For small rotations (< 5°) the deviation from axial is minimal and acceptance is unlikely to be an issue.

## No-Distortion Guarantee

A geometric distortion would occur if different regions of the image were scaled, sheared, or mapped non-linearly. The transformation used here is strictly rigid: $\mathbf{R}$ is orthonormal ($\det \mathbf{R} = 1$, $\|\mathbf{R}\,\mathbf{v}\| = \|\mathbf{v}\|$ for all $\mathbf{v}$), so all distances and angles between any two anatomical points are preserved exactly in patient space. The inverse mapping matrix $\mathbf{M} = \mathbf{A}^{-1}\mathbf{T}^{-1}\mathbf{A}$ is also an affine map with orthonormal linear part, so no shearing or anisotropic scaling is introduced at the voxel level either.

## 3D Visualisation

### Surface Extraction: Marching Cubes

The 3D body surface is extracted from the CT volume using the **Marching Cubes algorithm** (Lorensen & Cline, 1987). For each $2 \times 2 \times 2$ cube of adjacent voxels, the algorithm determines which of the 256 possible configurations of voxels above/below the iso-threshold is present and places triangular surface patches accordingly. This produces a polygonal mesh approximating the iso-surface at the chosen HU threshold.

Two surfaces are extracted:
- **Body surface**: threshold −300 HU, separating soft tissue from air
- **Bone surface**: threshold +400 HU, isolating cortical bone

A **downsampling factor of 2** is applied before running Marching Cubes (every second voxel in each direction) to reduce computation time and triangle count. This halves the spatial resolution of the surface mesh but has no effect on the underlying DICOM data.

### Coordinate Conversion

The Marching Cubes algorithm returns vertices in voxel index space $(k, j, i)$. These are converted to patient coordinates in millimetres by the affine matrix $\mathbf{A}$:

$$\mathbf{p}_{\text{patient}} = \mathbf{A}\, \begin{pmatrix} k \\ j \\ i \\ 1 \end{pmatrix}$$

For the metadata-only method, where pixel data is unchanged, the transformed surface is obtained by applying the forward transformation $\mathbf{T}$ to the original vertices directly — no second surface extraction from a resampled volume is needed.

### Interactive Visualisation

The extracted meshes are rendered using **Plotly's Mesh3d** trace with Gouraud-style lighting. All layers (original body, transformed body, original bone, rotation axes) are independently toggleable via the legend. The scene uses `aspectmode='data'` to ensure that spatial distances are displayed without distortion, i.e., 1 mm in X, Y, and Z corresponds to the same pixel length on screen.

## Summary of Design Decisions

| Decision | Choice | Rationale |
|---|---|---|
| Mapping direction | Inverse (pull-back) | Avoids holes and compositing problems inherent to forward mapping |
| Default interpolation | Trilinear (order 1) | Clinically standard; no overshoot; < 2 HU error in soft tissue |
| Rotation convention | Extrinsic XYZ | Each angle acts on a fixed patient axis, independent of the others |
| Rotation centre | Volume centroid | Body stays within the output FOV for typical small-angle corrections |
| Output geometry | Same grid as input | Standard axial IOP preserved; compatible with all TPS without reconfiguration |
| HU back-conversion | Per-slice slope/intercept | Preserves original calibration; no re-calibration artefacts |
| UIDs | New Series + SOP UIDs | TPS/PACS recognises output as independent series; no confusion with original |

## Literature

6. **DICOM Standard**, Part 3, Sections C.7.6.2, C.7.6.3 – *Image Plane Module, Image Pixel Module.* [dicom.nema.org](https://www.dicomstandard.org/)
7. Lehmann, T. M., Gönner, C., & Spitzer, K. (1999). *Survey: Interpolation methods in medical image processing.* IEEE Transactions on Medical Imaging, 18(11), 1049–1075.
8. Lorensen, W. E., & Cline, H. E. (1987). *Marching cubes: A high resolution 3D surface construction algorithm.* ACM SIGGRAPH Computer Graphics, 21(4), 163–169.
9. Thévenaz, P., Blu, T., & Unser, M. (2000). *Interpolation revisited.* IEEE Transactions on Medical Imaging, 19(7), 739–758.
10. Sled, J. G., Zijdenbos, A. P., & Evans, A. C. (1998). *A nonparametric method for automatic correction of intensity nonuniformity in MRI data.* IEEE Transactions on Medical Imaging, 17(1), 87–97.

---

# RTSTRUCT Visualizer Documentation

## Overview

The **RTSTRUCT Visualizer** (`visualizer.py`) generates a set of purpose-built static plots from the analysis results produced by the RTSTRUCT Analyzer. Each plot is designed to answer a specific clinical question and avoids chart types that are misleading at the typical structure counts encountered in radiotherapy planning (usually 5–30 structures per patient). The visualizer calls `run_analysis` internally, so only the RTSTRUCT file path is required — no intermediate JSON file is needed.

## Output Files and Clinical Motivation

### 1. `volumes.png` – Structure Volume Bar Chart

**Chart type:** Horizontal bar chart, sorted descending by volume.

A histogram (the previous implementation) is only meaningful when many samples are drawn from an unknown distribution, which is never the case here: each bar represents one named anatomical structure. A sorted bar chart makes it immediately apparent which structures are largest, how Targets and OARs compare in size, and whether any structure has an unexpectedly small or large volume (which can indicate a contouring error). Targets are shown in blue, OARs in red. Volume values are annotated on each bar.

**Clinical relevance:** Volume is the primary metric for Target coverage (PTV volume drives monitor unit calculation) and OAR sparing (e.g., mean brain dose correlates with brain volume irradiated). Unexpected outliers in volume are a common QA flag.

### 2. `shape_metrics.png` – Shape Metrics Comparison

**Chart type:** Three grouped bar charts side by side (one panel per metric).

Plotting all three shape metrics (sphericity, compactness, elongation) in a single figure with a consistent structure ordering allows direct visual comparison. A reference line at 1.0 is drawn for sphericity and compactness, since 1.0 is the theoretical maximum for a convex, sphere-like structure.

**Sphericity** quantifies how closely the structure resembles a sphere: values near 1 indicate round, regular shapes; values well below 1 indicate irregular or elongated structures. Clinically, low sphericity in a PTV may indicate a complex shape requiring more beam arrangements.

**Compactness** (volume / convex hull volume) detects concavity: a value below ~0.85 suggests the structure wraps around other anatomy (e.g., a C-shaped PTV around the brainstem). This is important when choosing between conformal arc and IMRT/VMAT techniques.

**Elongation** (sqrt of largest-to-smallest PCA eigenvalue ratio) measures directional stretching. High elongation combined with low sphericity in an OAR such as the spinal cord confirms its cylindrical nature, which is expected. Unexpectedly high elongation in a GTV may indicate a drawing artefact.

### 3. `distances.png` – Target–OAR Distance Bar Chart

**Chart type:** Grouped bar chart with three distance types per pair, limited to the 25 Target–OAR pairs with the smallest minimum distance (most critical first).

**Design decisions:**
- **Only Target–OAR pairs** are shown. Target–Target and OAR–OAR distances are less clinically meaningful for plan optimisation; they can always be retrieved from `statistics.txt`.
- **Sorted ascending by minimum distance**: the most critical proximity relationships appear on the left.
- **Cap at 25 pairs**: prevents figure overflow. With 65 structures, all pairwise combinations produce over 2000 entries — a bar chart of that size is unreadable (198,660 × 600 px) and clinically useless.
- **5 mm threshold line**: AAPM TG-218 and most institutional protocols flag structures within 5 mm of a Target boundary as requiring explicit dose–volume constraint review. This line immediately highlights which OARs fall into the critical zone.

**Three distance types** are shown simultaneously per pair:
- *Minimum distance*: the closest point between two contour surfaces; 0 mm means overlap.
- *Hausdorff distance*: the worst-case separation; indicates maximum excursion of one structure towards the other.
- *Centroid distance*: robust gross separation; useful as a sanity check (should always exceed minimum distance).

### 4. `centroids_3d.png` – Spatial Centroid Map

**Chart type:** 3D scatter plot using matplotlib's mpl_toolkits.mplot3d.

All structure centroids are plotted in the DICOM patient coordinate system (X=left, Y=posterior, Z=superior). Marker size scales with the square root of structure volume, so large structures are visually prominent without dominating. Targets use filled circles, OARs use triangles.

**Clinical relevance:** This plot answers the question "where is everything relative to everything else?" at a glance. It is particularly useful for multi-metastasis cases (e.g., the example dataset has 7 brain metastases): one can verify that all GTV/PTV pairs are spatially co-located and that the OARs (brainstem, optic chiasm, cochleae) are in the expected positions.

**Limitation:** The 3D scatter is static (not interactive). For interactive exploration, the CT Transformer's `visualization_3d.html` (Plotly) is the better tool.

### 5. `statistics.txt` – Numerical Summary

A structured plain-text file with per-structure details (volume, centroid, bounding box, all three shape metrics) and aggregated distance statistics (mean, std, min, max for each distance type, plus a full pairwise table). Suitable for copy-paste into clinical reports or further spreadsheet analysis.

## Implementation Notes

### Structure Filtering and Auto-Detection

If `--targets` or `--oars` are not specified on the CLI, the visualizer relies on the DICOM `RTROIInterpretedType` tag to classify structures. Recognised target types: `PTV`, `CTV`, `GTV`, `TV`. Recognised OAR types: `OAR`, `ORGAN`, `AVOIDANCE`. Structures of type `MARKER`, `EXTERNAL`, or `SUPPORT` are silently ignored in all plots (they appear in the structure list but are not analysed unless explicitly named).

### Scalability

The plots are designed to remain readable up to approximately 30 structures. Above that, the shape metrics panel becomes crowded and the 3D centroid map overlapping labels may need manual adjustment. The distance plot is always limited to 25 pairs regardless of total structure count.

### Matplotlib Backend

`matplotlib.use("Agg")` is set at module level so the visualizer can run on headless servers (e.g., CI pipelines, remote compute nodes) without an X display. All output is written to files; no interactive window is opened.

## Literature

11. **AAPM Task Group 218** (2021). *Tolerance limits and methodologies for IMRT measurement-based verification QA.* Medical Physics, 48(10).
12. Taha, A. A., & Hanbury, A. (2015). *Metrics for evaluating 3D medical image segmentation: analysis, selection, and tool.* BMC Medical Imaging, 15(1), 29.
