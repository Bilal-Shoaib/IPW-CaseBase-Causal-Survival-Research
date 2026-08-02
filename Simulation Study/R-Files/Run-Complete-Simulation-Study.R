# ==============================================================================
# Run Complete Simulation Study
# ==============================================================================

## Run the complete simulation study for a specified sample size, treatment
## effect, and censoring level. The function evaluates both RMST and hazard
## ratio performance for the MS-Cox and Case-Base estimators, with optional
## generation and saving of diagnostic plots.
run.complete.ss <- function(
    n,
    beta.A,
    censoring.limit,
    marginal.estimates,
    generate.plots = FALSE,
    save.plots = FALSE
) {
  
  # ----------------------------------------------------------------------------
  # Determine Whether Plots Are Required
  # ----------------------------------------------------------------------------
  
  ## Generate plots when either display or file output has been requested
  need.plots <- generate.plots || save.plots
  
  # ----------------------------------------------------------------------------
  # Run Simulation Replicates
  # ----------------------------------------------------------------------------
  
  ## Generate the simulation results for the specified data-generating scenario
  ss.out <- simulate(
    n, r, l,
    treatment.predictor.coefs,
    intercept,
    outcome.predictor.coefs,
    beta.A,
    censoring.limit,
    constant.HR
  )
  
  # ==============================================================================
  # Restricted Mean Survival Time (RMST)
  # ==============================================================================
  
  # ----------------------------------------------------------------------------
  # Plot RMST Results
  # ----------------------------------------------------------------------------
  
  ## Generate RMST plots comparing the two estimators with the true marginal
  ## RMST curve when plotting has been requested
  if (need.plots) {
    plot.all.rmst(
      marginal.estimates$rmst,
      ss.out$cox.rmst,
      ss.out$cb.rmst,
      save.plots = save.plots
    )
  }
  
  # ----------------------------------------------------------------------------
  # Select Restriction Times for Summary
  # ----------------------------------------------------------------------------
  
  ## Identify integer restriction times greater than zero for tabular summaries
  tau.points <- which(
    (t.grid != 0) &
      (t.grid == floor(t.grid))
  )
  
  # ----------------------------------------------------------------------------
  # Calculate True Marginal RMST Values
  # ----------------------------------------------------------------------------
  
  ## Calculate the Monte Carlo mean and standard deviation of the marginal RMST
  ## values used as the true reference values
  marginal.rmst.means <- rowMeans(
    marginal.estimates$rmst[tau.points, ]
  )
  
  marginal.rmst.sd <- apply(
    marginal.estimates$rmst[tau.points, ],
    1,
    sd
  )
  
  # ----------------------------------------------------------------------------
  # Summarize RMST Simulation Performance
  # ----------------------------------------------------------------------------
  
  ## Calculate simulation performance measures for the MS-Cox RMST estimator
  cox.rmst.summary <- get.rmst.summary(
    tau.points,
    ss.out$cox.rmst,
    ss.out$cox.rmst.se,
    marginal.rmst.means
  )
  
  ## Calculate simulation performance measures for the Case-Base RMST estimator
  cb.rmst.summary <- get.rmst.summary(
    tau.points,
    ss.out$cb.rmst,
    ss.out$cb.rmst.se,
    marginal.rmst.means
  )
  
  ## Define the restriction times reported in the simulation summary
  tau.values <- seq(1, 4, by = 1)
  
  # ----------------------------------------------------------------------------
  # Print RMST Results
  # ----------------------------------------------------------------------------
  
  ## Print the true marginal RMST values at the selected restriction times
  print.true.rmst(
    marginal.estimates$rmst,
    tau.points,
    tau.values
  )
  
  ## Print the MS-Cox RMST simulation performance summary
  print.rmst.summary(
    cox.rmst.summary,
    tau.values,
    "MS-Cox RMST"
  )
  
  ## Print the Case-Base RMST simulation performance summary
  print.rmst.summary(
    cb.rmst.summary,
    tau.values,
    "Case Base RMST"
  )
  
  # ==============================================================================
  # Hazard Ratio (HR)
  # ==============================================================================
  
  # ----------------------------------------------------------------------------
  # Determine the True Marginal HR
  # ----------------------------------------------------------------------------
  
  ## Calculate the true marginal hazard ratio by averaging the simulated
  ## marginal log-HR values and transforming to the HR scale
  true.hr <- exp(mean(marginal.estimates$psi))
  
  # ----------------------------------------------------------------------------
  # Plot HR Results
  # ----------------------------------------------------------------------------
  
  if (need.plots) {
    
    ## Plot the distribution of the marginal HR values used to determine the
    ## simulation reference value
    if (save.plots) {
      pdf("true-hr.pdf", width = 8, height = 5)
    }
    
    hist(
      exp(marginal.estimates$psi),
      breaks = 30,
      xlab = "HR",
      main = "Marginal HR"
    )
    
    if (save.plots) {
      dev.off()
    }
    
    
    ## Compare the sampling distributions of the MS-Cox and Case-Base HR
    ## estimates against the true marginal HR
    plot.hr(
      true.hr,
      exp(ss.out$cox.psi),
      exp(ss.out$cb.psi),
      save.plot = save.plots
    )
    
    
    ## Compare the estimation errors of the two HR estimators relative to the
    ## true marginal HR
    plot.hr.errors(
      true.hr,
      exp(ss.out$cox.psi),
      exp(ss.out$cb.psi),
      save.plot = save.plots
    )
  }
  
  # ----------------------------------------------------------------------------
  # Print True HR
  # ----------------------------------------------------------------------------
  
  ## Report the true marginal HR used as the reference value for evaluating
  ## estimator performance
  print.true.hr(true.hr)
  
  # ----------------------------------------------------------------------------
  # Summarize HR Simulation Performance
  # ----------------------------------------------------------------------------
  
  ## Calculate and print simulation performance measures for both HR estimators
  all.hr.stats <- summarize.print.hr(
    true.hr,
    ss.out$cox.psi,
    ss.out$cox.psi.se,
    ss.out$cb.psi,
    ss.out$cb.psi.se
  )
  
  # ----------------------------------------------------------------------------
  # Plot HR Confidence Intervals
  # ----------------------------------------------------------------------------
  
  ## Visualize confidence intervals from both estimators across simulation
  ## replicates, including the true HR and null value
  if (need.plots) {
    plot.all.hr.ci(
      true.hr,
      exp(ss.out$cox.psi),
      all.hr.stats[[1]],
      exp(ss.out$cb.psi),
      all.hr.stats[[2]],
      seed = seed,
      save.plot = save.plots
    )
  }
}