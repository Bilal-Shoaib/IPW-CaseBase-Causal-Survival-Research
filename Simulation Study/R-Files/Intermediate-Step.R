# ==============================================================================
# Intermediate Simulation Step
# ==============================================================================

## Repeatedly estimate the large-sample marginal treatment effect and RMST
## contrast across independent simulated datasets. The computation is
## parallelized across multiple worker processes to reduce runtime.
intermediate.step <- function(
    r, l, 
    intercept,
    outcome.predictor.coefs,
    beta.A,
    censoring.limit,
    constant.HR
) {
  
  # ----------------------------------------------------------------------------
  # Parallel Computing Setup
  # ----------------------------------------------------------------------------
  
  ## Create a two-worker parallel cluster and register it for use by foreach
  cl <- parallel::makeCluster(2)
  doParallel::registerDoParallel(cl)
  
  ## Ensure that the cluster is stopped when the function exits, including if
  ## an error occurs during the simulation
  on.exit(
    parallel::stopCluster(cl),
    add = TRUE
  )
  
  # ----------------------------------------------------------------------------
  # Repeated Marginal Effect Estimation
  # ----------------------------------------------------------------------------
  
  ## Repeat the marginal-effect calculation over 100 independently generated
  ## datasets. The random-number-generation option ensures reproducible
  ## results across parallel workers.
  result <- foreach(
    i = 1:100,
    .packages = "survival",
    .export = c(
      "generate.data",
      "get.marginal.estimates",
      "t.grid",
      "rmst.fn"
    ),
    .options.RNG = seed
  ) %dorng% {
    
    ## Generate the large simulated dataset and obtain the corresponding
    ## marginal log hazard ratio and RMST contrast
    get.marginal.estimates(
      r, l,
      intercept,
      outcome.predictor.coefs,
      beta.A,
      censoring.limit,
      constant.HR
    )
  }
  
  # ----------------------------------------------------------------------------
  # Result Restructuring
  # ----------------------------------------------------------------------------
  
  ## Convert the list of simulation results into the required structured
  ## format for subsequent simulation analyses
  return(
    refactor(result)
  )
}