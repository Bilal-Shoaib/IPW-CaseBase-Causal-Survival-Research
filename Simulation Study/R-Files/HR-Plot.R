# ==============================================================================
# Plot Hazard Ratio Distributions
# ==============================================================================

## Plot the sampling distributions of the estimated hazard ratios from the
## MS-Cox and Case-Base estimators and compare them with the true marginal HR.
plot.hr <- function(
    true.hr,
    cox.hr,
    cb.hr,
    save.plot = FALSE
) {
  
  # ----------------------------------------------------------------------------
  # Set Up Optional PDF Output
  # ----------------------------------------------------------------------------
  
  ## Open a PDF device when the plot is requested to be saved
  if (save.plot) {
    pdf("hr.pdf", width = 8, height = 5)
  }
  
  # ----------------------------------------------------------------------------
  # Plot Estimated HR Distributions
  # ----------------------------------------------------------------------------
  
  ## Plot the sampling distribution of the MS-Cox HR estimates
  hist(
    cox.hr,
    breaks = 50,
    col = rgb(1, 0, 0, 0.5),
    border = "red",
    xlab = "Hazard Ratio",
    main = "Hazard Ratio Comparison"
  )
  
  ## Overlay the sampling distribution of the Case-Base HR estimates
  hist(
    cb.hr,
    breaks = 50,
    col = rgb(0, 0, 1, 0.5),
    border = "blue",
    add = TRUE
  )
  
  # ----------------------------------------------------------------------------
  # Add True HR Reference
  # ----------------------------------------------------------------------------
  
  ## Mark the true marginal HR with a vertical reference line
  abline(
    v = true.hr,
    col = "black",
    lwd = 3,
    lty = 2
  )
  
  # ----------------------------------------------------------------------------
  # Add Legend
  # ----------------------------------------------------------------------------
  
  ## Identify the two estimator distributions and the true marginal HR
  legend(
    "topright",
    legend = c(
      "MS-Cox HR",
      "Case Base HR",
      "True HR"
    ),
    fill = c(
      rgb(1, 0, 0, 0.5),
      rgb(0, 0, 1, 0.5),
      NA
    ),
    border = c(
      "red",
      "blue",
      NA
    ),
    lty = c(
      NA,
      NA,
      2
    ),
    lwd = c(
      NA,
      NA,
      3
    ),
    col = c(
      NA,
      NA,
      "black"
    )
  )
  
  # ----------------------------------------------------------------------------
  # Close Optional PDF Output
  # ----------------------------------------------------------------------------
  
  ## Close the PDF device after the plot has been written
  if (save.plot) {
    dev.off()
  }
}

# ==============================================================================
# Plot Hazard Ratio Estimation Errors
# ==============================================================================

## Plot the sampling distributions of the estimation errors for the MS-Cox and
## Case-Base HR estimators. An error of zero indicates an estimate equal to the
## true marginal HR.
plot.hr.errors <- function(
    true.hr,
    cox.hr,
    cb.hr,
    save.plot = FALSE
) {
  
  # ----------------------------------------------------------------------------
  # Calculate HR Estimation Errors
  # ----------------------------------------------------------------------------
  
  ## Calculate the deviation of each estimated HR from the true marginal HR
  cox.hr.error <- cox.hr - true.hr
  cb.hr.error <- cb.hr - true.hr
  
  # ----------------------------------------------------------------------------
  # Set Up Optional PDF Output
  # ----------------------------------------------------------------------------
  
  ## Open a PDF device when the plot is requested to be saved
  if (save.plot) {
    pdf("hr-errors.pdf", width = 8, height = 5)
  }
  
  # ----------------------------------------------------------------------------
  # Plot HR Error Distributions
  # ----------------------------------------------------------------------------
  
  ## Plot the sampling distribution of the MS-Cox HR estimation errors
  hist(
    cox.hr.error,
    breaks = 50,
    col = rgb(1, 0, 0, 0.5),
    border = "red",
    xlab = "HR Errors",
    main = "HR Errors Comparison"
  )
  
  ## Overlay the sampling distribution of the Case-Base HR estimation errors
  hist(
    cb.hr.error,
    breaks = 50,
    col = rgb(0, 0, 1, 0.5),
    border = "blue",
    add = TRUE
  )
  
  # ----------------------------------------------------------------------------
  # Add Zero-Error Reference
  # ----------------------------------------------------------------------------
  
  ## Mark zero error, corresponding to an estimate exactly equal to the true
  ## marginal HR
  abline(
    v = 0,
    col = "black",
    lwd = 3,
    lty = 2
  )
  
  # ----------------------------------------------------------------------------
  # Add Legend
  # ----------------------------------------------------------------------------
  
  ## Identify the estimation-error distributions for the two estimators
  legend(
    "topleft",
    legend = c(
      "MS-Cox HR Error",
      "Case Base HR Error"
    ),
    fill = c(
      rgb(1, 0, 0, 0.5),
      rgb(0, 0, 1, 0.5)
    ),
    border = c(
      "red",
      "blue"
    )
  )
  
  # ----------------------------------------------------------------------------
  # Close Optional PDF Output
  # ----------------------------------------------------------------------------
  
  ## Close the PDF device after the plot has been written
  if (save.plot) {
    dev.off()
  }
}