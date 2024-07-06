# MALAT1_threshold

Low _MALAT1_ expression is associated with a lack of a nucleus in single-cell RNA-sequencing data. Cells without nuclei are likely either empty droplets filled with ambient RNA, cell fragments, or mature erythrocytes. This function takes a vector of normalized _MALAT1_ expression, and outputs a minimum threshold value that can be used to filter your scRNA-seq object.
