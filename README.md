# DICOM File Modifier

A comprehensive toolkit for analyzing, modifying, and visualizing DICOM RT Structure Set files used in radiotherapy planning. This project provides tools to process target volumes (PTV, CTV, GTV) and organs at risk (OAR) from DICOM files, compute geometric metrics, and generate visualizations.

## Features

- **Analyzer**: Extract and compute geometric properties (volume, centroid, shape metrics, distances) from RTSTRUCT files
- **Modifier**: Modify DICOM structure sets (planned for future release)
- **Visualizer**: Generate plots and statistics from analysis results

## Project Structure

```
dicom-file-modifier/
├── data/                    # Input DICOM files (not synced)
├── output/                  # Analysis results and modified files (not synced)
├── dicom_file_modifier/     # Python package
│   ├── __init__.py
│   ├── analyzer.py          # RTSTRUCT analysis module
│   ├── modifier.py          # DICOM modification module (future)
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

Generate plots and statistics from analysis results:

```bash
python -m dicom_file_modifier.visualizer output/analysis_results.json --output output/
```

This creates:
- `volume_histogram.png`: Volume distribution of structures
- `distance_scatter.png`: Scatter plot of inter-structure distances
- `sphericity_vs_volume.png`: Shape analysis plot
- `statistics.txt`: Summary statistics

### Modifier

*Coming soon* - Tools for modifying DICOM structure sets.

## Dependencies

- pydicom: DICOM file handling
- numpy: Numerical computations
- scipy: Scientific computing (distances, convex hull)
- shapely: 2D geometry operations
- matplotlib: Plotting and visualization

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
| `scipy` | ≥ 1.7 | KD-tree, Hausdorff distance, convex hull |
| `shapely` | ≥ 1.8 | 2D polygon operations (area calculation) |
| `matplotlib` | ≥ 3.5 | Optional visualization |

## Literature

1. **ICRU Report 50** (1993). *Prescribing, Recording, and Reporting Photon Beam Therapy.*
2. **ICRU Report 62** (1999). *Prescribing, Recording and Reporting Photon Beam Therapy (Supplement to ICRU Report 50).*
3. **ICRU Report 83** (2010). *Prescribing, Recording, and Reporting Photon-Beam Intensity-Modulated Radiation Therapy (IMRT).*
4. **DICOM Standard**, Part 3, Section C.8.8.6 – *RT Structure Set Module.* [dicom.nema.org](https://www.dicomstandard.org/)
5. Huttenlocher, D. P., Klanderman, G. A., & Rucklidge, W. J. (1993). *Comparing images using the Hausdorff distance.* IEEE Transactions on Pattern Analysis and Machine Intelligence, 15(9), 850–863.
