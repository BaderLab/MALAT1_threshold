#' @import gamlss
#' @importFrom stats optimize
NULL

run_gamlss <- function(malat1_vec,
                       counts,
                       min_malat1 = 1,
                       normalized = FALSE,
                       print_plots = TRUE,
                       return_plots = FALSE) {

  success <- FALSE
  attempt <- 1
  max_attempts <- 5

  if (length(malat1_vec) < 10) {
    stop("Too few MALAT1 values remain after filtering")
  }

  malat1_df <- data.frame(y=malat1_vec)

  ctr <- gamlss.mx::MX.control(plot = FALSE)

  while(!success && attempt <= max_attempts) {

    fit <- tryCatch(
      {
        gamlss.mx::gamlssMX(y ~ 1,
                            data = malat1_df,
                            family = get("GA", envir = asNamespace("gamlss")),
                            K = 2,
                            control = ctr)
      },
      error = function(e) {
        message("Trying again with greater absolute minimum:",
                e$message)
        NULL
      })

    if(!is.null(fit)) {
      success <- TRUE

    } else {

      attempt <- attempt + 1
      min_malat1 <- min_malat1 + 1

      if(normalized) {
        malat1_vec <- counts[counts >= min_malat1]
      } else {
        malat1_vec <- log1p(counts[counts >= min_malat1])
      }

      if (length(malat1_vec) < 10) {
        stop("No MALAT1 values remain after applying the minimum expression threshold")
      }

      malat1_df <- data.frame(y=malat1_vec)
    }
  }

  if(!success){
    stop("Model failed after maximum attempts")
  }

  mean1 <- fit$models[[1]]$mu.coefficients
  mean2 <- fit$models[[2]]$mu.coefficients

  if(fit$models[[1]]$mu.link == "log") mean1 <- exp(mean1)
  if(fit$models[[2]]$mu.link == "log") mean2 <- exp(mean2)

  lower <- min(mean1, mean2)
  upper <- max(mean1, mean2)

  fn <- gamlss.mx::getpdfMX(fit,
                            observation=1)

  threshold <- optimize(fn,
                        interval=c(lower, upper))$minimum

  p_gamlss <- ggplot2::ggplot(malat1_df,
                              ggplot2::aes(x = malat1_df$y)) +

    # The Histogram
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                            bins = 50,
                            fill = "grey90",
                            color = "white") +
    # The gamma curve
    ggplot2::stat_function(fun = fn,
                           color = "darkgreen",
                           linewidth = 1.5) +
    # Threshold
    ggplot2::geom_vline(xintercept = threshold,
                        color = "darkorange",
                        linewidth = 1,
                        linetype = "dashed") +

    ggplot2::labs(title = paste("MALAT1 Mixture Fit"),
         x = "Log Normalized Counts",
         y = "Density") +

    ggplot2::theme_classic()

  if (isTRUE(x = print_plots)) {
    print(p_gamlss)
  }

  if (isTRUE(x = return_plots)) {
    result <- list("threshold" = threshold,
                "plot" = p_gamlss)
  } else {
    result <- list("threshold" = threshold)
  }

  return(result)

}

malat1_thres_raw <- function(raw_counts,
                             barcodes,
                             min_malat1 = 4,
                             print_plots = TRUE,
                             return_plots = FALSE) {

  if(length(raw_counts) != length(barcodes)) {
    stop("Counts vector and barcodes vector must be the same length.")
  }

  malat1_vec <- log1p(raw_counts[raw_counts >= min_malat1])

  result <- run_gamlss(malat1_vec,
                       counts = raw_counts,
                       min_malat1 = min_malat1,
                       normalized = FALSE,
                       print_plots = print_plots,
                       return_plots = return_plots)

  cells <- barcodes[malat1_vec > result[["threshold"]]]

  result[["cells"]] <- cells

  return(result)

}

malat1_thres_filtered <- function(counts,
                                  barcodes,
                                  normalized = FALSE,
                                  min_malat1 = 1,
                                  breaks = 100,
                                  print_plots = TRUE,
                                  return_plots = FALSE) {

  if(length(counts) != length(barcodes)) {
    stop("Counts vector and barcodes vector must be the same length.")
  }

  if(normalized) {
    malat1_vec <- counts
  } else {
    malat1_vec <- log1p(counts)
  }

  # Before truncating vector, plot histogram of normalized MALAT1 expression
  malat1_df <- data.frame(y = malat1_vec)
  p_hist <- ggplot2::ggplot(malat1_df, ggplot2::aes(x = malat1_df$y)) +
    ggplot2::geom_histogram(bins = breaks,
                            fill = "grey90",
                            color = "white") +
    ggplot2::labs(x = "Normalized MALAT1 expression",
                  y = "Number of cells") +
    ggplot2::theme_classic()

  if (isTRUE(x = print_plots)) {
    print(p_hist)
  }

  if(normalized) {
    malat1_vec <- malat1_vec[malat1_vec >= min_malat1]
  } else {
    malat1_vec <- log1p(counts[counts >= min_malat1])
  }

  result <- run_gamlss(malat1_vec,
                       min_malat1 = min_malat1,
                       counts = counts,
                       normalized = normalized,
                       print_plots = print_plots,
                       return_plots = return_plots)

  if (isTRUE(x = return_plots)) {
    result[["histogram"]] <- p_hist
  }

  cells <- barcodes[malat1_vec > result[["threshold"]]]

  result[["cells"]] <- cells

  return(result)

}

malat1_read_10X_raw <- function(raw_cellranger_folder,
                                unique.features = TRUE,
                                malat1_name = "MALAT1",
                                min_malat1 = 4,
                                print_plots = TRUE,
                                return_plots = FALSE) {

  sobj.data <- Seurat::Read10X(raw_cellranger_folder,
                               unique.features = unique.features)

  counts <- sobj.data[malat1_name,]

  result <- malat1_thres_raw(counts, names(counts),
                             min_malat1 = 4,
                             print_plots = TRUE,
                             return_plots = FALSE)

  return(result)

}

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

  counts <- sobj.data[malat1_name,]

  result <- malat1_thres_raw(counts, names(counts),
                             min_malat1 = 4,
                             print_plots = TRUE,
                             return_plots = FALSE)

  return(result)

}
