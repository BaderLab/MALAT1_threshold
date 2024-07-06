# How to apply a _MALAT1_ threshold for your scRNA-seq object in R

Low _MALAT1_ expression is associated with a lack of a nucleus in single-cell RNA-sequencing data. Cells without nuclei are likely either empty droplets filled with ambient RNA, cell fragments, or mature erythrocytes. Our function `define_malat1_threshold` takes a vector of normalized _MALAT1_ expression, and outputs a minimum threshold value that can be used to filter your scRNA-seq object.
