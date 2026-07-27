#' Filter cells based on MALAT1 expression from a raw count matrix containing cells and droplets
#'
#' @param raw_counts A vector of MALAT1 counts.
#' @param barcodes A vector of cell barcodes.
#' @param min_malat1 The minimum MALAT1 count required for a value to be included in the mixture model fitting.
#' @param print_plots Logical indicating whether or not to print MALAT1 histogram with model fit.
#' @param return_plots Logical indicating whether or not to return MALAT1 histogram with model fit.
#' @export
#' @examples
#'
#' data(example_counts)
#' data(example_barcodes)
#'
#' result <- malat1_thres_raw(
#'   example_counts,
#'   example_barcodes,
#'   print_plots = TRUE
#' )
#'
#' result$threshold
#' result$cells[1:5]

malat1_thres_raw <- function(raw_counts,
                             barcodes,
                             min_malat1 = 4,
                             print_plots = TRUE,
                             return_plots = FALSE) {

  if(length(raw_counts) != length(barcodes)) {
    stop("Counts vector and barcodes vector must be the same length.")
  }

  names(raw_counts) <- barcodes

  malat1_vec <- log1p(raw_counts[raw_counts >= min_malat1])

  result <- run_gamlss(malat1_vec,
                       original_counts = raw_counts,
                       min_malat1 = min_malat1,
                       log_transform = TRUE,
                       print_plots = print_plots,
                       return_plots = return_plots)

  result[["cells"]] <- names(malat1_vec)

  return(result)

}

#' Filter cells based on MALAT1 expression from a filtered dataset of cells that may include empty droplets
#'
#' @param counts A vector of MALAT1 counts.
#' @param barcodes A vector of cell barcodes.
#' @param log_transform Whether to apply a log1p transformation to MALAT1 counts before fitting the mixture model.
#' @param min_malat1 The minimum MALAT1 count required for a value to be included in the mixture model fitting.
#' @param breaks Number of bins used to generate the MALAT1 histogram.
#' @param print_plots Logical indicating whether or not to print MALAT1 histogram with model fit.
#' @param return_plots Logical indicating whether or not to return MALAT1 histogram with model fit.
#' @export
#' @examples
#'
#' data(example_counts)
#' data(example_barcodes)
#'
#' result <- malat1_thres_filtered(
#'   example_counts,
#'   example_barcodes,
#'   print_plots = TRUE
#' )
#'
#' result$threshold
#' result$cells[1:5]

malat1_thres_filtered <- function(counts,
                                  barcodes,
                                  log_transform = TRUE,
                                  min_malat1 = 1,
                                  breaks = 100,
                                  print_plots = TRUE,
                                  return_plots = FALSE) {

  if(length(counts) != length(barcodes)) {
    stop("Counts vector and barcodes vector must be the same length.")
  }

  names(counts) <- barcodes

  if(log_transform) {
    malat1_vec <- log1p(counts)
  } else {
    malat1_vec <- counts
  }

  # Before truncating vector, plot histogram of transformed MALAT1 expression
  malat1_df <- data.frame(y = malat1_vec)
  p_hist <- ggplot2::ggplot(malat1_df, ggplot2::aes(x = y)) +
    ggplot2::geom_histogram(bins = breaks,
                            fill = "grey90",
                            color = "white") +
    ggplot2::labs(x = "Transformed MALAT1 expression",
                  y = "Number of cells") +
    ggplot2::theme_classic()

  if (isTRUE(x = print_plots)) {
    print(p_hist)
  }

  if(log_transform) {
    malat1_vec <- log1p(counts[counts >= min_malat1])
  } else {
    malat1_vec <- malat1_vec[malat1_vec >= min_malat1]
  }

  result <- run_gamlss(malat1_vec,
                       min_malat1 = min_malat1,
                       original_counts = counts,
                       log_transform = log_transform,
                       print_plots = print_plots,
                       return_plots = return_plots)

  if (isTRUE(x = return_plots)) {
    result[["histogram"]] <- p_hist
  }

  result[["cells"]] <- names(malat1_vec)

  return(result)

}
