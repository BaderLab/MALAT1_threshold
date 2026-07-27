# Function to calculate MALAT1 threshold

# Takes the following parameters as input:
# - counts: A vector of normalized MALAT1 counts; this is ideally from a single unintegrated sample consisting of multiple cell types
# - bw: The "bandwidth" value when plotting the density function to the MALAT1 distribution; default is bw = 0.1,
# but this parameter should be lowered (e.g. to 0.01) if you run the function and the line that's produced doesn't look like it's
# tracing the shape of the histogram accurately (this will make the line less "stiff" and more fitted to the data)
# - lwd: The "line width" fed to the abline function which adds the vertical red line to the output plots; default is 2, and it
# can be increased or decreased depending on the user's plotting preferences
# - breaks: The number of bins used for plotting the histogram of normalized MALAT1 values; default is 100
# - chosen_min: The minimum MALAT1 value cutoff above which a MALAT1 peak in the density function should be found. This value is
# necessary to determine which peak in the density function fitted to the MALAT1 distribution is likely representative of what we
# would expect to find in real cells. This is because some samples may have large numbers of cells or empty droplets with lower than
# expected normalized MALAT1 values, and therefore have a peak close to or at zero. Ideally, "chosen_min" would be manually chosen
# after looking at a histogram of MALAT1 values, and be the normalized MALAT1 value that cuts out all of the cells that look like they
# stray from the expected distribution (a unimodal distribution above zero). The default value is 1 as this works well in many test
# cases, but different types of normalization may make the user want to change this parameter (e.g. Seurat's original normalization
# function generates different results to their SCT function) which may change the MALAT1 distribution).
# Increase or decrease chosen_min depending on where your MALAT1 peak is
# - smooth: The "smoothing parameter" fed into the "smooth.spline" function that adjusts the trade-off between the smoothness of the
# line fitting the histogram, and how closely it fits the histogram; the default is 1, and can be lowered if it looks like the line
# is underfitting the data, and raised in the case of overfitting. The ideal scenario is for the line to trace the histogram in a way
# where the only inflection point(s) are between major peaks, e.g. separating the group of poor-quality cells or empty droplets with
# lower normalized MALAT1 expression from higher-quality cells with higher normalized MALAT1 expression.
# - abs_min: The absolute lowest value allowed as the MALAT1 threshold. This parameter increases the robustness of the function if
# working with an outlier data distribution (e.g. an entire sample is poor quality so there is a unimodal MALAT1 distribution that is
# very low but above zero, but also many values close to zero) and prevents a resulting MALAT1 threshold of zero. In the case where
# a calculated MALAT1 value is zero, the function will return 0.3 by default.
# - rough_max: A rough value for the location of a MALAT1 peak if a peak is not found. This is possible if there are so few cells with
# higher MALAT1 values, that a distribution fitted to the data finds no local maxima. For example, if a sample only has poor-quality
# cells such that all have near-zero MALAT1 expression, the fitted function may look similar to a positive quadratic function which
# has no local maxima. In this case, the function searches for the closest MALAT1 value to the default value, 2, to use in place of
# a real local maximum.

define_malat1_threshold_ggplot2 <- function(counts, bw = 0.1, lwd = 2, breaks = 100,
                                    chosen_min = 1, smooth = 1, abs_min = 0.3,
                                    rough_max = 2, print_plots = TRUE, return_plots = FALSE) {
  tryCatch({
    # Calculate the density values
    density_data <- density(counts, bw = bw,
                           from = min(counts),
                            to = max(counts))
    # Fit a smooth line to the data
    fit <- smooth.spline(density_data$x, density_data$y, spar = smooth)
    # Predict the first derivative
    first_derivative <- predict(fit, density_data$x, deriv = 1)

    # Find the local maxima
    local_maxima <- density_data$x[which(diff(sign(first_derivative$y)) == -2)]
    if(length(local_maxima) == 0) {
      # Find x-val closest to rough_max
      local_maxima <- density_data$x[which(abs(density_data$x - rough_max) == min(abs(density_data$x - rough_max)))]
    }
    # Plot the density with the local maxima
    plot(density_data, main = "Density Plot with Local Maxima")
    lines(fit, col = "blue")
    points(local_maxima, predict(fit, local_maxima)$y, col = "red", pch = 19)

    # Create data frames for plotting
    density_df <- data.frame(x = density_data$x, y = density_data$y)
    fit_df <- data.frame(x = density_data$x, y = predict(fit, density_data$x)$y)
    local_maxima_df <- data.frame(x = local_maxima, y = predict(fit, local_maxima)$y)

    # Plot the density with the local maxima using ggplot2
    p1 <- ggplot() +
      geom_line(data = density_df, aes(x = x, y = y), color = "black") +
      geom_line(data = fit_df, aes(x = x, y = y), color = "blue") +
      geom_point(data = local_maxima_df, aes(x = x, y = y), color = "red", size = 3) +
      ggtitle("Density Plot with Local Maxima") +
      cowplot::theme_cowplot() +
      ylim(-0.05, NA)

    # Find the local minima
    local_minima <- density_data$x[which(diff(sign(first_derivative$y)) == 2)]
    if(length(local_minima) == 0) {
      local_minima <- abs_min
    }
    # Plot the density with the local minima
    plot(density_data, main = "Density Plot with Local Minima")
    lines(fit, col = "blue")
    points(local_minima, predict(fit, local_minima)$y, col = "red", pch = 19)

    # Create data frame for local minima
    local_minima_df <- data.frame(x = local_minima, y = predict(fit, local_minima)$y)

    # Plot the density with the local minima using ggplot2
    p2 <- ggplot() +
      geom_line(data = density_df, aes(x = x, y = y), color = "black") +
      geom_line(data = fit_df, aes(x = x, y = y), color = "blue") +
      geom_point(data = local_minima_df, aes(x = x, y = y), color = "red", size = 3) +
      ggtitle("Density Plot with Local Minima") +
      cowplot::theme_cowplot() +
      ylim(-0.05, NA)

    # Find biggest local maximum greater than desired value (default 1)
    biggest_y <- max(density_data$y[density_data$x %in% local_maxima[local_maxima > chosen_min]])
    maxi <- density_data$x[density_data$y == biggest_y]

    # Find local minimum closest to the left of that
    local_minima <- local_minima[local_minima < maxi]
    if(length(local_minima) < abs_min) {
      local_minima <- abs_min
    }
    mini <- local_minima[(maxi - local_minima) == min(maxi - local_minima)]

    # Calculate delta to get range of values to isolate for quadratic calculation
    delta <- maxi - mini

    # Subset dataframe
    df <- as.data.frame(cbind(density_data$x, density_data$y))
    colnames(df) <- c("x","y")
    subset_df <- df[df[,1] >= (maxi - delta) & df[,1] <= (maxi + delta), ]

    # Fit a quadratic model (y ~ x + I(x^2))
    quad_model <- lm(y ~ poly(x, 2, raw = TRUE), data = subset_df)

    # Plot quadratic
    plot(df$x, df$y,
         xlab = "Normalized MALAT1 expression",
         ylab = "Density value")
    # Add the extracted subset data
    points(subset_df$x, subset_df$y, pch = 16, col = "blue")
    # Add the fitted quadratic curve
    curve(predict(quad_model, newdata = data.frame(x = x)), add = TRUE, col = "red", lwd = 2)

    p3 <- ggplot(df, aes(x = x, y = y)) +
      geom_line(color = "black") +
      geom_point(data = subset_df, aes(x = x, y = y), color = "blue", size = 2) +
      stat_function(fun = function(x) predict(quad_model, newdata = data.frame(x = x)), color = "red", size = 1) +
      labs(
        title = "Density Plot with Quadratic Fit",
        x = "Normalized MALAT1 expression",
        y = "Density value"
      ) +
      cowplot::theme_cowplot() +
      ylim(-0.05, NA)

    # Grab intercept
    # Extract coefficients from the summary
    coefficients <- summary(quad_model)$coefficients
    # Extract coefficients
    c <- coefficients[1, 1]
    b <- coefficients[2, 1]
    a <- coefficients[3, 1]
    # Calculate discriminant
    discriminant <- b^2 - 4 * a * c
    # Grab first intercept
    x_intercept1 <- (-b + sqrt(discriminant)) / (2 * a)
    if(x_intercept1 < 0) {
      x_intercept1 <- abs_min
    }

    # Plot histogram with threshold
    hist(counts, breaks = breaks,
         xlab = "Normalized MALAT1 expression",
         ylab = "Number of cells")
    abline(v = x_intercept1, col = "red", lwd = 2)

    p4 <- ggplot(data = data.frame(counts), aes(x = .data[["counts"]])) +
      geom_histogram(color = "black", bins = breaks, fill = "dodgerblue") +
      cowplot::theme_cowplot() +
      geom_vline(xintercept = x_intercept1, color = "red", linewidth = lwd) +
      ggtitle("MALAT1")

    plots <- wrap_plots(p1, p2, p3, p4, ncol = 2)

    if (isTRUE(x = print_plots)) {
      print(plots)
    }

    if (isTRUE(x = return_plots)) {
      res <- list("threshold" = x_intercept1,
                  "plots" = plots)
      return(res)
    } else {
      return(x_intercept1)
    }

  }, error = function(e) {
    # Code to execute if an error occurs
    message(" An error occurred: Please make sure you have use a vector of normalized counts as input. This may also indicate that you have no high MALAT1 peaks, meaning this particular sample may be poor quality. Please check your histogram of normalized MALAT1 counts to investigate: ", e$message)
    hist(counts, breaks = breaks,
         xlab = "Normalized MALAT1 expression",
         ylab = "Number of cells")

    ggplot(data = data.frame(counts), aes(x = .data[["counts"]])) +
      geom_histogram(color = "black", bins = breaks, fill = "dodgerblue") +
      cowplot::theme_cowplot() +
      ggtitle("MALAT1")

    return(2)
  })
}
