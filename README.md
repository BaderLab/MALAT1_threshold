# How to apply a _MALAT1_ threshold for your scRNA-seq object in R

For a detailed explanation of our findings or citation of this work, please see our preprint on BioRxiv: https://doi.org/10.1101/2024.07.14.603469

Low _MALAT1_ expression is associated with a lack of a nucleus in single-cell RNA-sequencing data. Cells without nuclei are likely either empty droplets filled with ambient RNA, cell fragments, or mature erythrocytes. Our function `define_malat1_threshold` takes a vector of normalized _MALAT1_ expression, and outputs a minimum threshold value that can be used to filter your scRNA-seq object.

We hope to develop a package to allow a user to easily access this function. In the meantime, you can use this the function by either pasting the code directly into your R script, or cloning the GitHub repo, moving the `malat1_function.R` script into your analysis directory, and adding `source("malat1_function.R")` to the top of your script to access the function.

To use this function, isolate the normalized _MALAT1_ expression values from your scRNA-seq object. In a Seurat object, this may look like:

```
norm_counts <- sobj@assays$RNA@data["MALAT1",]
```

This can be fed into the _MALAT1_ threshold function which will return the minimum _MALAT1_ value that each cell should contain:

```
threshold <- define_malat1_threshold(norm_counts)
```

This threshold value can be used to flag or filter cells from your single-cell object. The code below flags cells that don't pass the threshold by using `TRUE` values to represent good cells, and `FALSE` to represent cells that don't pass the filter:

```
malat1_threshold <- norm_counts > threshold
sobj$malat1_threshold <- malat1_threshold
sobj$malat1_threshold <- factor(sobj$malat1_threshold, levels = c("TRUE","FALSE"))
DimPlot(sobj, group.by = "malat1_threshold")
```

From this, you can use the result to remove cells from your object:

```
good_cells <- WhichCells(sobj, expression = malat1_threshold == TRUE)
good_sobj <- subset(sobj, cells = good_cells)
```

## Example analysis

We can demonstrate using this function with the Tabula Muris Senis Pancreas dataset which can be downloaded as a Seurat object from [cellxgene](https://cellxgene.cziscience.com/collections/0b9d8a04-bb9d-44da-aa27-705bb65b54eb). The data are also described in this paper:

The Tabula Muris Consortium. A single-cell transcriptomic atlas characterizes ageing tissues in the mouse. Nature 583, 590–595 (2020). https://doi.org/10.1038/s41586-020-2496-1

Here is a brief look at the cells in the dataset:

<img width="400" alt="tabula_muris_senis_pancreas_celltypes" src="https://github.com/user-attachments/assets/bc4a6e67-640c-44d3-b699-0d63860d83bc">

Here is _MALAT1_ projected onto the UMAP, and a histogram of the normalized _MALAT1_ values for that dataset. You can see that certain cells in the pancreatic acinar cell cluster have especially low _MALAT1_ values and may be suffering from some quality issues. Cells that may actually be empty droplets would be those in the lower _MALAT1_ expression peak in the histogram, in addition to those with a peak at zero:

<img width="400" alt="tabula_muris_senis_pancreas_malat1_umap" src="https://github.com/user-attachments/assets/5ed0ba68-efc4-45cd-bdf4-d4a9e403a6f2">
<img width="397" alt="tabula_muris_senis_pancreas_malat1_hist_noLine" src="https://github.com/user-attachments/assets/5ba81742-0295-4605-a8f0-7384f108dd32">

This function fits a density function to the histogram, and models a quadratic to the highest _MALAT1_ expression peak above the normalized expression value of two. It finds this peak by analysing local minima and maxima that appear on the density function. The lower x-intercept of this quadratic is used to define the minimum _MALAT1_ threshold.

The function outputs the following plots: (1) The density plot with local minima annotated. (2) The density plot with local maxima annotated. (3) The points of the density function in black, with points highlighted in blue covering the range of the data that the quadratic is fit to, with the quadratic fit overtop of the points in red (below). (4) The histogram of normalized _MALAT1_ counts with the red line indicating the minimum threshold value (below).

<img width="400" alt="tabula_muris_senis_pancreas_malat1_quad" src="https://github.com/user-attachments/assets/96e7ad95-1b39-461f-b561-7f4c3a4efc3f">
<img width="372" alt="tabula_muris_senis_pancreas_malat1_hist" src="https://github.com/user-attachments/assets/ace9f285-5b80-4387-b388-ff1b24f6e92c">

Using the code above, we can see which cells passed the filter. Most of the cells that failed the filter are, in fact, pancreatic acinar cells (highlighted below as "FALSE" for having not passed the filter):

<img width="400" alt="tabula_muris_senis_pancreas_malat1_dimplot" src="https://github.com/user-attachments/assets/64049f75-6b54-44b3-95a1-620e61f6b346">

## Troubleshooting

This analysis relies on the assumption that there is a _MALAT1_ peak above the normalized value of two. If such a peak (i.e. local maximum above two) does not exist, the function may call an error. This is probably a good sign to take a closer look at your data anyway, but you can also lower this value by adjusting the parameter `chosen_min`.

Other parameters that can be modified are `bw`, `lwd`, and `breaks`. Increasing or decreasing `bw` to say 0.5 or 0.01 respectively will change the plotting of the density function, with higher values creating a function with fewer inflection points (i.e. a "less curvy" function). Modifying `lwd` changes the thickness of the line on the final plotted histogram, and `breaks` is the number of buckets used in the histogram.

Worst case scenario, if the function doesn't work for some weird, confusing reason, you can always eyeball your _MALAT1_ values to try and figure out if there is something fishy going on with your data. You can manually choose your own threshold by looking at the histogram, or just pick out clusters of concerning cells by projecting _MALAT1_ onto your UMAP.



