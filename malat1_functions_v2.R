run_gamlss <- function(malat1_vec, min_malat1 = min_malat1, 
                       counts = counts, normalized = normalized) {

  success <- FALSE
  attempt <- 1
  max_attempts <- 5

  malat1_df <- data.frame(y=malat1_vec)

  ctr <- MX.control(plot = FALSE)

  while(!success && attempt <= max_attempts) {
    
    fit <- tryCatch(
      {
        gamlssMX(y ~ 1, data = malat1_df, family = GA, K = 2, control = ctr)
      },
      error = function(e) {
        print("Trying again with greater absolute minimum")
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

  fn <- getpdfMX(fit, observation=1)

  threshold <- optimize(fn, interval=c(lower, upper))$minimum

  p <- ggplot(malat1_df, aes(x = y)) +
  
    # The Histogram
    geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "grey90", color = "white") +
    # The gamma curve     
    stat_function(fun = fn, color = "darkgreen", linewidth = 1.5) +
    # Threshold
    geom_vline(xintercept = threshold, color = "darkorange", linewidth = 1, linetype = "dashed") +
  
    labs(title = paste("MALAT1 Mixture Fit"),
         x = "Log Normalized Counts", y = "Density") +
  
    theme_classic()
  
  print(p)

  return(threshold)
  
}

malat1_thres_raw <- function(raw_counts, malat1_name = "MALAT1", min_malat1 = 4) {

  run_gamlss(raw_counts, min_malat1 = min_malat1)
  
}

malat1_thres_filtered <- function(counts, normalized = FALSE, min_malat1 = 1, breaks = 100) {
  
  if(normalized) {
    malat1_vec <- counts
  } else {
    malat1_vec <- log1p(counts)
  }

  # Before truncating vector, plot histogram of normalized MALAT1 expression
  print(hist(malat1_vec, breaks = breaks,
     xlab = "Normalized MALAT1 expression",
     ylab = "Number of cells"))

  if(normalized) {
    malat1_vec <- malat1_vec[malat1_vec >= min_malat1]
  } else {
    malat1_vec <- log1p(counts[counts >= min_malat1])
  }

  run_gamlss(malat1_vec, min_malat1 = min_malat1)
  
}
