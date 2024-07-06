# Function to calculate MALAT1 threshold

# Counts is vector of normalized MALAT1 counts
# Lower bw if minimum doesn't look totally accurate
# Change lwd to change thickness of red line
# Change breaks to change bins of histogram
# Change chosen_min if you have a very low malat1 peak, lower than 2, that you want to keep
define_malat1_threshold <- function(counts, bw = 0.1, lwd = 2, breaks = 100, chosen_min = 2) {
  tryCatch({
    # Calculate the density values
    density_data <- density(counts, bw = bw)
    # Fit a smooth line to the data
    fit <- smooth.spline(density_data$x, density_data$y, spar = 0.5)
    # Predict the first derivative
    first_derivative <- predict(fit, density_data$x, deriv = 1)
    
    # Find the local minima
    local_minima <- density_data$x[which(diff(sign(first_derivative$y)) == 2)]
    # Plot the density with the local minima
    plot(density_data, main = "Density Plot with Local Minima")
    lines(fit, col = "blue")
    points(local_minima, predict(fit, local_minima)$y, col = "red", pch = 19)
    
    # Find the local maxima
    local_maxima <- density_data$x[which(diff(sign(first_derivative$y)) == -2)]
    # Plot the density with the local maxima
    plot(density_data, main = "Density Plot with Local Maxima")
    lines(fit, col = "blue")
    points(local_maxima, predict(fit, local_maxima)$y, col = "red", pch = 19)
    
    # Find biggest local maximum greater than desired value (default 2)
    biggest_y <- max(density_data$y[density_data$x %in% local_maxima[local_maxima > chosen_min]])
    maxi <- density_data$x[density_data$y == biggest_y]
    
    # Find local minimum closest to the left of that
    local_minima <- local_minima[local_minima < maxi]
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
    
    # Plot histogram with threshold
    hist(counts, breaks = breaks,
         xlab = "Normalized MALAT1 expression",
         ylab = "Number of cells")
    abline(v = x_intercept1, col = "red", lwd = 2)
    
    return(x_intercept1)
  }, error = function(e) {
    # Code to execute if an error occurs
    message(" An error occurred: Please make sure you have use a vector of normalized counts as input. This may also indicate that you have no high MALAT1 peaks, meaning this particular sample may be poor quality. Please check your histogram of normalized MALAT1 counts to investigate: ", e$message)
    hist(counts, breaks = breaks,
         xlab = "Normalized MALAT1 expression",
         ylab = "Number of cells")
    return(2)
  })
}
