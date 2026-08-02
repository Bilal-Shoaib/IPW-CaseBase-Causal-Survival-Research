# ==============================================================================
# Run Simulation Study
# ==============================================================================

## Run the complete simulation study for a specified data-generating mechanism.
## Each simulation replicate independently generates data, estimates stabilized
## IPWs, fits the MS-Cox and IPW case-base models, and calculates the treatment
## effect and RMST estimates.
simulate <- function(
    n, r, l,
    treatment.predictor.coefs,
    intercept,
    outcome.predictor.coefs,
    beta.A,
    censoring.limit,
    constant.HR
) {
  
  # ----------------------------------------------------------------------------
  # Parallel Processing Setup
  # ----------------------------------------------------------------------------
  
  ## Allocate approximately 60% of the available CPU cores to the parallel
  ## simulation to balance computational speed with system resource usage
  n.cores <- round(parallel::detectCores() * 0.6)
  
  cl <- parallel::makeCluster(n.cores)
  doParallel::registerDoParallel(cl)
  
  ## Ensure that the parallel cluster is stopped when the function exits,
  ## including if an error occurs during the simulation
  on.exit(
    parallel::stopCluster(cl),
    add = TRUE
  )
  
  # ----------------------------------------------------------------------------
  # Run Simulation Replicates
  # ----------------------------------------------------------------------------
  
  ## Run 1,000 independent simulation replicates in parallel. The RNG option
  ## ensures reproducible random-number generation across parallel workers.
  result <- foreach(
    i = 1:1000,
    .packages = c("survival", "casebase"),
    .export = c(
      "run.one.sim",
      "generate.data",
      "vcovHC", "vcov",
      "FSH", "expand.dot.formula",
      "t.grid", "rmst.fn",
      "cox.rmst.list",
      "cb.rmst.list",
      "shape"
    ),
    .options.RNG = seed
  ) %dorng% {
    
    ## Execute a single simulation replicate using the specified
    ## data-generating mechanism and estimation settings
    run.one.sim(
      n,
      r,
      l,
      treatment.predictor.coefs,
      intercept,
      outcome.predictor.coefs,
      beta.A,
      censoring.limit,
      constant.HR,
      i
    )
  }
  
  # ----------------------------------------------------------------------------
  # Organize Simulation Results
  # ----------------------------------------------------------------------------
  
  ## Convert the list of replicate-specific results into a compact structure
  ## suitable for summarizing estimator performance across simulations
  return(
    refactor(result)
  )
}