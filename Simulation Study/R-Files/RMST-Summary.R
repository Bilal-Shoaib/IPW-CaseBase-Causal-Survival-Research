# ==============================================================================
# RMST Simulation Performance Summary
# ==============================================================================

## Summarize the finite-sample performance of an RMST estimator across
## simulation replicates. The function calculates the Monte Carlo mean,
## empirical standard deviation, mean estimated standard error, bias, RMSE,
## and confidence-interval coverage at the specified restriction times.
##
## The marginal RMST values are treated as the target estimand against which
## the simulation estimates are evaluated.
get.rmst.summary <- function(
    tau.points,
    rmst.matrix,
    rmst.se.matrix,
    marginal.rmst.means,
    alpha = 0.05
) {
  
  ## Number of simulation replicates and standard-normal critical value for
  ## the nominal (1 - alpha) confidence intervals
  n <- ncol(rmst.matrix)
  z <- qnorm(1 - alpha / 2)
  
  # ----------------------------------------------------------------------------
  # Point Estimation and Variability
  # ----------------------------------------------------------------------------
  
  ## Monte Carlo mean of the estimated RMST contrast across simulation
  ## replicates at each selected restriction time
  rmst.means <- rowMeans(
    rmst.matrix[tau.points, ]
  )
  
  ## Empirical standard deviation of the RMST estimates across simulation
  ## replicates, representing the observed Monte Carlo variability
  rmst.sd <- apply(
    rmst.matrix[tau.points, ],
    1,
    sd
  )
  
  ## Average model-based standard error reported by the estimator across
  ## simulation replicates
  rmst.se.mean <- rowMeans(
    rmst.se.matrix[tau.points, ]
  )
  
  # ----------------------------------------------------------------------------
  # Bias and RMSE
  # ----------------------------------------------------------------------------
  
  ## Monte Carlo bias relative to the marginal RMST target
  rmst.biases <- rmst.means - marginal.rmst.means
  
  ## Monte Carlo standard error of the estimated bias
  rmst.bias.mcse <- rmst.sd / sqrt(n)
  
  ## Root mean squared error relative to the marginal RMST target, combining
  ## both systematic bias and sampling variability
  rmst.rmse <- sqrt(
    rowMeans(
      (
        rmst.matrix[tau.points, ] -
          marginal.rmst.means
      )^2
    )
  )
  
  # ----------------------------------------------------------------------------
  # Confidence-Interval Coverage
  # ----------------------------------------------------------------------------
  
  ## Construct nominal (1 - alpha) Wald confidence intervals for the RMST
  ## contrast using the estimated standard errors from each simulation
  lower <- rmst.matrix[tau.points, ] -
    z * rmst.se.matrix[tau.points, ]
  
  upper <- rmst.matrix[tau.points, ] +
    z * rmst.se.matrix[tau.points, ]
  
  ## Estimate empirical coverage as the proportion of simulation replicates
  ## whose confidence interval contains the marginal RMST target
  coverage <- rowMeans(
    lower <= marginal.rmst.means &
      upper >= marginal.rmst.means
  )
  
  ## Monte Carlo standard error of the estimated coverage probability, treating
  ## coverage as a binomial proportion across simulation replicates
  coverage.mcse <- sqrt(
    coverage * (1 - coverage) / n
  )
  
  # ----------------------------------------------------------------------------
  # Assemble Summary Table
  # ----------------------------------------------------------------------------
  
  ## Initialize the summary matrix with one row per restriction time and one
  ## column for each simulation performance metric
  k <- length(tau.points)
  
  rmst.summary <- matrix(
    NA_real_,
    nrow = k,
    ncol = 8
  )
  
  colnames(rmst.summary) <- c(
    "mean",
    "sd",
    "se.mean",
    "bias",
    "bias.mcse",
    "rmse",
    "coverage",
    "coverage.mcse"
  )
  
  ## Populate the summary matrix with the calculated performance measures
  rmst.summary[, "mean"] <- rmst.means
  rmst.summary[, "sd"] <- rmst.sd
  rmst.summary[, "se.mean"] <- rmst.se.mean
  rmst.summary[, "bias"] <- rmst.biases
  rmst.summary[, "bias.mcse"] <- rmst.bias.mcse
  rmst.summary[, "rmse"] <- rmst.rmse
  rmst.summary[, "coverage"] <- coverage
  rmst.summary[, "coverage.mcse"] <- coverage.mcse
  
  return(rmst.summary)
}

# ==============================================================================
# Print RMST Simulation Summary
# ==============================================================================

## Print the simulation performance summary for an RMST estimator in a
## human-readable format. Results are displayed separately for each selected
## restriction time, including the Monte Carlo mean, bias, RMSE, estimated
## standard error, empirical standard deviation, and confidence-interval
## coverage.
print.rmst.summary <- function(
    rmst.summary,
    tau.values,
    model
) {
  
  # ----------------------------------------------------------------------------
  # Print Model Header
  # ----------------------------------------------------------------------------
  
  ## Display the name of the estimator being summarized
  cat("\n========================================\n")
  cat(model, "\n")
  cat("========================================\n")
  
  cat("\nRestricted Mean Survival Time Summary\n")
  cat("----------------------------------------\n")
  
  # ----------------------------------------------------------------------------
  # Print Results at Each Restriction Time
  # ----------------------------------------------------------------------------
  
  ## Display the simulation performance metrics separately at each selected
  ## restriction time
  for (i in seq_along(tau.values)) {
    
    cat("\nTau =", tau.values[i], "\n")
    cat("----------------------------------------\n")
    
    ## Monte Carlo mean of the estimated RMST contrast
    cat(
      "Mean RMST:",
      rmst.summary[i, "mean"],
      "\n"
    )
    
    ## Monte Carlo bias and its Monte Carlo standard error
    cat(
      "Bias:",
      rmst.summary[i, "bias"],
      " | MC-SE:",
      rmst.summary[i, "bias.mcse"],
      "\n"
    )
    
    ## Root mean squared error of the RMST estimator
    cat(
      "RMSE:",
      rmst.summary[i, "rmse"],
      "\n"
    )
    
    ## Compare the average model-based standard error with the empirical
    ## standard deviation across simulation replicates
    cat(
      "Mean Estimated SE:",
      rmst.summary[i, "se.mean"],
      " | Empirical SD:",
      rmst.summary[i, "sd"],
      "\n"
    )
    
    ## Empirical confidence-interval coverage and its Monte Carlo standard error
    cat(
      "Coverage:",
      rmst.summary[i, "coverage"],
      " | MC-SE:",
      rmst.summary[i, "coverage.mcse"],
      "\n"
    )
  }
  
  cat("\n")
}

# ==============================================================================
# Print True Marginal RMST
# ==============================================================================

## Print the true marginal RMST contrast at each selected restriction time.
## The true value is obtained by averaging the simulated marginal RMST values
## across the corresponding time point.
print.true.rmst <- function(
    true.rmst,
    tau.points,
    tau.values
) {
  
  # ----------------------------------------------------------------------------
  # Print Header
  # ----------------------------------------------------------------------------
  
  ## Identify the values being reported as the true marginal RMST
  cat("\n========================================\n")
  cat("\nTrue Mean Restricted Mean Survival Time\n")
  cat("========================================\n")
  
  
  # ----------------------------------------------------------------------------
  # Print True RMST at Each Restriction Time
  # ----------------------------------------------------------------------------
  
  ## Display the true marginal RMST corresponding to each selected restriction
  ## time
  for (i in seq_along(tau.values)) {
    cat(
      tau.values[i],
      ":",
      mean(true.rmst[tau.points[i], ]),
      "\n"
    )
  }
  
  cat("\n")
}