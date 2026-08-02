# ==============================================================================
# Large-Sample Marginal Effect Estimates
# ==============================================================================

## Generate a large simulated dataset and estimate the marginal treatment
## effect using an unadjusted Cox model. The resulting estimates provide
## large-sample benchmarks for evaluating the performance of the estimators
## used in the simulation study.
get.marginal.estimates <- function(
    r, l,
    intercept,
    outcome.predictor.coefs,
    beta.A,
    censoring.limit,
    constant.HR
) {
  
  # ----------------------------------------------------------------------------
  # Generate Large Simulated Dataset
  # ----------------------------------------------------------------------------
  
  ## Generate a large dataset to approximate the marginal treatment effect.
  ## Setting all treatment-model coefficients to zero removes confounding from
  ## treatment assignment while retaining the covariate effects in the outcome
  ## model.
  dat <- generate.data(
    n                         = 10000000,
    r                         = r,
    l                         = l,
    treatment.predictor.coefs = rep(0, 10),
    intercept                 = intercept,
    outcome.predictor.coefs   = outcome.predictor.coefs,
    beta.A                    = beta.A,
    censoring.limit           = censoring.limit,
    constant.HR               = constant.HR
  )
  
  # ----------------------------------------------------------------------------
  # Marginal Cox Model
  # ----------------------------------------------------------------------------
  
  ## Fit an unadjusted Cox model to estimate the marginal treatment effect.
  ## Under the constant-HR specification, the coefficient of A provides the
  ## large-sample marginal log hazard ratio used as the simulation benchmark.
  ms.cox <- NULL
  psi <- NULL
  
  if (constant.HR) {
    
    ms.cox <- coxph(
      Surv(time, Delta) ~ A,
      data = dat
    )
    
    ## Extract the estimated marginal log hazard ratio for treatment
    psi <- coef(ms.cox)["A"]
    
    # --------------------------------------------------------------------------
    # Time-Varying Treatment Effect
    # --------------------------------------------------------------------------
    
  } else {
    
    ## Under the time-varying specification, allow the treatment effect to
    ## vary with follow-up time through the A-by-time interaction
    ms.cox <- coxph(
      Surv(time, Delta) ~ A * time,
      data = dat
    )
    
    ## No single scalar marginal log hazard ratio is extracted when the
    ## treatment effect varies over time
    psi <- NULL
  }
  
  # ----------------------------------------------------------------------------
  # Counterfactual Survival Curves
  # ----------------------------------------------------------------------------
  
  ## Estimate the marginal counterfactual survival curves under no treatment
  ## and treatment by evaluating the fitted model at A = 0 and A = 1
  sf.cox <- survfit(
    ms.cox,
    newdata = data.frame(A = c(0, 1))
  )
  
  ## Interpolate the counterfactual survival curve under A = 0 onto the common
  ## simulation time grid
  s0.cox <- approx(
    sf.cox$time,
    sf.cox$surv[, 1],
    xout   = t.grid,
    method = "constant",
    rule   = 2
  )$y
  
  ## Interpolate the counterfactual survival curve under A = 1 onto the common
  ## simulation time grid
  s1.cox <- approx(
    sf.cox$time,
    sf.cox$surv[, 2],
    xout   = t.grid,
    method = "constant",
    rule   = 2
  )$y
  
  # ----------------------------------------------------------------------------
  # Marginal RMST Contrast
  # ----------------------------------------------------------------------------
  
  ## Calculate the marginal RMST contrast by integrating the difference between
  ## the two counterfactual survival curves over the simulation time grid
  output <- list(
    psi  = psi,
    rmst = rmst.fn(t.grid, s1.cox) - rmst.fn(t.grid, s0.cox)
  )
  
  return(output)
}