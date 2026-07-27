#' @import gamlss
#' @import gamlss.dist
#' @importFrom stats optimize
#' @keywords internal
NULL

run_gamlss <- function(malat1_vec,
                       original_counts,
                       min_malat1 = 1,
                       log_transform = TRUE,
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
                            family = GA,
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

      if(log_transform) {
        malat1_vec <- log1p(original_counts[original_counts >= min_malat1])
      } else {
        malat1_vec <- original_counts[original_counts >= min_malat1]
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
                              ggplot2::aes(x = y)) +

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
                  x = "Log-Transformed Counts",
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
