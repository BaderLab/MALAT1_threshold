# How to apply a _MALAT1_ threshold for your scRNA-seq object in R

Low _MALAT1_ expression is associated with a lack of a nucleus in single-cell RNA-sequencing data. Cells without nuclei are likely either empty droplets filled with ambient RNA, cell fragments, or mature erythrocytes. Our function `define_malat1_threshold` takes a vector of normalized _MALAT1_ expression, and outputs a minimum threshold value that can be used to filter your scRNA-seq object.

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
