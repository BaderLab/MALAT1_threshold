# MALATin

**MALATin** is an R package for identifying and filtering low-MALAT1 cells in single-cell RNA-sequencing (scRNA-seq) data. Low MALAT1 expression is associated with cells lacking an intact nucleus, including empty droplets, cell fragments, and mature erythrocytes. The package automatically estimates a data-driven MALAT1 expression threshold using a two-component gamma mixture model and returns the cells that pass this quality-control filter.

This package can be used to extract cells from a raw matrix of cells + empty droplets (e.g. as an alternative to EmptyDrops or CellRanger), or on a filtered matrix of cells but before clustering and differential expression analysis.

For a detailed description of the method, please see our preprint:

Clarke et al. (2024). MALAT1 expression identifies low-quality cells in single-cell RNA sequencing datasets. BioRxiv. https://doi.org/10.1101/2024.07.14.603469

If you encounter unexpected behaviour or have suggestions for improvements, please open a GitHub issue.

---

## Installation

Install the development version from GitHub using **remotes** or **devtools**:

# install.packages("remotes")
remotes::install_github("BaderLab/malatin")

---


## Quick start

### Filter an existing vector of _MALAT1_ counts

```
library(malatin)

data(example_counts)
data(example_barcodes)

# From an object already filtered with e.g. CellRanger or EmptyDrops
result <- malat1_thres_filtered(counts = example_counts, barcodes = example_barcodes)
# From a raw sample that has not yet been processed by e.g. CellRanger or EmptyDrops
# result <- malat1_thres_raw(counts = example_counts, barcodes = example_barcodes)

result$threshold head(result$cells)
```

The returned object is a list containing
- `threshold`: estimated _MALAT1_ threshold
- `cells`: barcodes of cells passing the filter
- `plot`: histogram with fitted mixture model (optional)
- `histogram`: histogram of MALAT1 expression (optional)

---

## Read CellRanger output directly

For raw Cell Ranger output:

```
result <- malat1_read_10X_raw(raw_cellranger_folder = "filtered_feature_bc_matrix/")
```

or from an HDF5 file:

```
result <- malat1_read_10X_h5_raw(raw_cellranger_h5 = "filtered_feature_bc_matrix.h5")
```

---

## Output

The package estimates a minimum MALAT1 expression threshold by fitting a two-component gamma mixture model to the MALAT1 distribution.

Depending on the function arguments, the returned object may include

- the estimated threshold,
- the barcodes of cells passing the filter,
- a histogram of MALAT1 expression,
- the fitted mixture model overlaid on the histogram.

These plots can be used to visually assess whether the inferred threshold is appropriate for your dataset.

---

## Included example data

Included example data

The package contains a small example dataset:

- `example_counts`
- `example_barcodes`

These datasets are used throughout the package examples and can also be used for testing installations.

---

## Citation

If you use MALATin in your research, please cite:

> Clarke et al. (2024). MALAT1 expression identifies low-quality cells in single-cell RNA sequencing datasets. BioRxiv. https://doi.org/10.1101/2024.07.14.603469

---

## Issues

If you encounter bugs or have suggestions for improvements, please open a GitHub Issue.

We are happy to help troubleshoot unusual datasets and improve the robustness of the package.

