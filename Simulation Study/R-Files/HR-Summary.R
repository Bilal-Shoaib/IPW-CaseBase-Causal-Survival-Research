# ==============================================================================
# Summarize Hazard Ratio Simulation Performance
# ==============================================================================

## Calculate simulation performance measures for a hazard ratio estimator,
## including bias, RMSE, average estimated standard error, confidence-interval
## coverage, and either Type I error or power.
get.hr.summary <- function(
    true.hr,
    psi.vector,
    se.vector,
    alpha = 0.05,
    power = TRUE
) {
  
  # ----------------------------------------------------------------------------
  # Set Up Simulation Quantities
  # ----------------------------------------------------------------------------
  
  ## Number of simulation replicates and normal critical value for the
  ## confidence intervals
  n <- length(psi.vector)
  z <- qnorm(1 - alpha / 2)
  
  ## Transform log-hazard-ratio estimates to the hazard-ratio scale
  hr.vector <- exp(psi.vector)
  
  # ----------------------------------------------------------------------------
  # Point-Estimation Performance
  # ----------------------------------------------------------------------------
  
  ## Monte Carlo mean, empirical standard deviation, and bias of the estimated
  ## hazard ratios
  hr.mean <- exp(mean(psi.vector))
  hr.emp.sd <- sd(hr.vector)
  hr.bias <- hr.mean - true.hr
  
  ## Monte Carlo standard error of the estimated bias
  hr.bias.mcse <- hr.emp.sd / sqrt(n)
  
  ## Mean squared error and root mean squared error relative to the true HR
  hr.mse <- mean((hr.vector - true.hr)^2)
  hr.rmse <- sqrt(hr.mse)
  
  # ----------------------------------------------------------------------------
  # Standard Error and Confidence-Interval Performance
  # ----------------------------------------------------------------------------
  
  ## Average model-based standard error across simulation replicates
  hr.se.mean <- mean(se.vector)
  
  ## Construct Wald confidence intervals on the HR scale
  hr.ci.lower <- exp(psi.vector - z * se.vector)
  hr.ci.upper <- exp(psi.vector + z * se.vector)
  
  # ----------------------------------------------------------------------------
  # Type I Error / Power and Coverage
  # ----------------------------------------------------------------------------
  
  ## Calculate the proportion of confidence intervals excluding the null HR of
  ## one. This corresponds to Type I error when the true HR equals one and to
  ## empirical power when the true HR differs from one.
  t1.error.or.power <- mean(
    hr.ci.lower > 1 | hr.ci.upper < 1
  )
  
  ## Calculate empirical coverage of the confidence intervals for the true HR
  hr.coverage <- mean(
    hr.ci.lower <= true.hr &
      true.hr <= hr.ci.upper
  )
  
  ## Monte Carlo standard error of the empirical coverage probability
  hr.coverage.mcse <- sqrt(
    hr.coverage * (1 - hr.coverage) / n
  )
  
  # ----------------------------------------------------------------------------
  # Store Summary Statistics and Confidence Intervals
  # ----------------------------------------------------------------------------
  
  ## Store the simulation performance measures and individual confidence
  ## intervals in a list for downstream analysis and reporting
  results <- list(
    summary = c(
      mean               = hr.mean,
      emp.sd             = hr.emp.sd,
      bias               = hr.bias,
      bias.mcse          = hr.bias.mcse,
      mse                = hr.mse,
      rmse               = hr.rmse,
      se.mean            = hr.se.mean,
      t1.error.or.power  = t1.error.or.power,
      coverage           = hr.coverage,
      coverage.mcse      = hr.coverage.mcse
    ),
    ci = cbind(
      lower = hr.ci.lower,
      upper = hr.ci.upper
    )
  )
  
  return(results)
}

# ==============================================================================
# Print Hazard Ratio Simulation Summary
# ==============================================================================

## Print the main simulation performance measures for a hazard ratio estimator
## in a human-readable format.
print.hr.summary <- function(
    hr.stats.summary,
    model
) {
  
  # ----------------------------------------------------------------------------
  # Print Model Header
  # ----------------------------------------------------------------------------
  
  ## Identify the estimator corresponding to the reported results
  cat("\n========================================\n")
  cat(model, "\n")
  cat("========================================\n")
  
  cat("\nHazard Ratio Summary\n")
  cat("----------------------------------------\n")
  
  # ----------------------------------------------------------------------------
  # Print Performance Measures
  # ----------------------------------------------------------------------------
  
  ## Monte Carlo mean and bias of the estimated hazard ratio
  cat(
    "Mean HR:",
    hr.stats.summary["mean"],
    "\n"
  )
  
  ## Bias and its Monte Carlo standard error
  cat(
    "Bias:",
    hr.stats.summary["bias"],
    " | MC-SE:",
    hr.stats.summary["bias.mcse"],
    "\n"
  )
  
  ## Mean squared error and root mean squared error
  cat(
    "MSE:",
    hr.stats.summary["mse"],
    "\n"
  )
  
  cat(
    "RMSE:",
    hr.stats.summary["rmse"],
    "\n"
  )
  
  ## Compare the average estimated standard error with the empirical standard
  ## deviation across simulation replicates
  cat(
    "Mean Estimated SE:",
    hr.stats.summary["se.mean"],
    " | Empirical SD:",
    hr.stats.summary["emp.sd"],
    "\n"
  )
  
  ## Report the proportion of confidence intervals excluding the null HR of one
  cat(
    "Type I Error/Power:",
    hr.stats.summary["t1.error.or.power"],
    "\n"
  )
  
  ## Report empirical confidence-interval coverage and its Monte Carlo standard
  ## error
  cat(
    "Coverage:",
    hr.stats.summary["coverage"],
    " | MC-SE:",
    hr.stats.summary["coverage.mcse"],
    "\n"
  )
  
  cat("\n")
}

# ==============================================================================
# Print True Hazard Ratio
# ==============================================================================

## Print the true marginal hazard ratio used as the reference value in the
## simulation study.
print.true.hr <- function(true.hr) {
  
  cat("\n========================================\n")
  cat("True HR:", true.hr, "\n")
  cat("========================================\n")
}

# ==============================================================================
# Calculate and Print HR Simulation Summaries
# ==============================================================================

## Calculate and print simulation performance summaries for both the MS-Cox and
## Case-Base estimators using their estimated log-HR coefficients and standard
## errors.
summarize.print.hr <- function(
    true.hr,
    cox.psi, cox.psi.se,
    cb.psi, cb.psi.se
) {
  
  # ----------------------------------------------------------------------------
  # Calculate Estimator-Specific Summaries
  # ----------------------------------------------------------------------------
  
  ## Summarize the simulation performance of the MS-Cox hazard ratio estimator
  cox.psi.stats <- get.hr.summary(
    true.hr,
    cox.psi,
    cox.psi.se
  )
  
  ## Summarize the simulation performance of the Case-Base hazard ratio
  ## estimator
  cb.psi.stats <- get.hr.summary(
    true.hr,
    cb.psi,
    cb.psi.se
  )
  
  # ----------------------------------------------------------------------------
  # Print Estimator Results
  # ----------------------------------------------------------------------------
  
  ## Print the MS-Cox simulation summary
  print.hr.summary(
    cox.psi.stats$summary,
    "MS-Cox"
  )
  
  ## Print the Case-Base simulation summary
  print.hr.summary(
    cb.psi.stats$summary,
    "Case Base"
  )
  
  # ----------------------------------------------------------------------------
  # Return Summary Objects
  # ----------------------------------------------------------------------------
  
  ## Return both estimator-specific results for downstream analysis
  return(list(
    cox.psi.stats,
    cb.psi.stats
  ))
}