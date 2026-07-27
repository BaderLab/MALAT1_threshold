#' Filter cells based on MALAT1 expression from a raw CellRanger folder.
#'
#' @param raw_cellranger_folder Folder containing barcodes.tsv, genes.tsv, and matrix.mtx
#' @param unique.features Logical to make feature names unique (default TRUE)
#' @param malat1_name Gene symbol used for MALAT1 (default "MALAT1"); use "Malat1" for mouse or Ensembl gene names if required
#' @param min_malat1 The minimum MALAT1 count required for a value to be included in the mixture model fitting.
#' @param print_plots Logical indicating whether or not to print MALAT1 histogram with model fit.
#' @param return_plots Logical indicating whether or not to return MALAT1 histogram with model fit.
#' @export

malat1_read_10X_raw <- function(raw_cellranger_folder,
                                unique.features = TRUE,
                                malat1_name = "MALAT1",
                                min_malat1 = 4,
                                print_plots = TRUE,
                                return_plots = FALSE) {

  sobj.data <- Seurat::Read10X(raw_cellranger_folder,
                               unique.features = unique.features)

  if (!malat1_name %in% rownames(sobj.data)) {
    stop(paste(malat1_name,"not found in dataset."))
    }

  counts <- sobj.data[malat1_name,]

  result <- malat1_thres_raw(counts, names(counts),
                             min_malat1 = min_malat1,
                             print_plots = print_plots,
                             return_plots = return_plots)

  return(result)

}

#' Filter cells based on MALAT1 expression from a raw CellRanger h5 file.
#'
#' @param raw_cellranger_h5 CellRanger output: raw_feature_bc_matrix_h5.h5.
#' @param use.names Row names with feature names rather than ID numbers (default TRUE).
#' @param unique.features Logical to make feature names unique (default TRUE)
#' @param malat1_name Gene symbol used for MALAT1 (default "MALAT1"); use "Malat1" for mouse or Ensembl gene names if required
#' @param min_malat1 The minimum MALAT1 count required for a value to be included in the mixture model fitting.
#' @param print_plots Logical indicating whether or not to print MALAT1 histogram with model fit.
#' @param return_plots Logical indicating whether or not to return MALAT1 histogram with model fit.
#' @export

malat1_read_10X_h5_raw <- function(raw_cellranger_h5,
                                   use.names = TRUE,
                                   unique.features = TRUE,
                                   malat1_name = "MALAT1",
                                   min_malat1 = 4,
                                   print_plots = TRUE,
                                   return_plots = FALSE) {

  sobj.data <- Seurat::Read10X_h5(raw_cellranger_h5,
                                  use.names = use.names,
                                  unique.features = unique.features)

  if (!malat1_name %in% rownames(sobj.data)) {
    stop(paste(malat1_name,"not found in dataset."))
  }

  counts <- sobj.data[malat1_name,]

  result <- malat1_thres_raw(counts, names(counts),
                             min_malat1 = min_malat1,
                             print_plots = print_plots,
                             return_plots = return_plots)

  return(result)

}
