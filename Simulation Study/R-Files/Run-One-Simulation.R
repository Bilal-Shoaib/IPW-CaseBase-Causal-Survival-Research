# ==============================================================================
# Single Simulation Replicate
# ==============================================================================

## Run one simulation replicate by generating a dataset, estimating stabilized
## inverse probability weights, fitting the MS-Cox and IPW case-base models,
## and calculating the corresponding treatment effects and RMST contrasts.
run.one.sim <- function(
    n, r, l,
    treatment.predictor.coefs,
    intercept,
    outcome.predictor.coefs,
    beta.A,
    censoring.limit,
    constant.HR,
    i
) {
  
  # ----------------------------------------------------------------------------
  # Generate Simulated Dataset
  # ----------------------------------------------------------------------------
  
  ## Generate one dataset under the specified data-generating mechanism
  data.set <- generate.data(
    n                         = n,
    r                         = r,
    l                         = l,
    treatment.predictor.coefs = treatment.predictor.coefs,
    intercept                 = intercept,
    outcome.predictor.coefs   = outcome.predictor.coefs,
    beta.A                    = beta.A,
    censoring.limit           = censoring.limit,
    constant.HR               = constant.HR
  )
  
  # ----------------------------------------------------------------------------
  # Propensity Score Estimation
  # ----------------------------------------------------------------------------
  
  ## Estimate the probability of treatment conditional on all baseline
  ## covariates using logistic regression
  p.score <- fitted(
    glm(
      A ~ X1 + X2 + X3 + X4 + X5 +
        X6 + X7 + X8 + X9 + X10,
      data = data.set,
      family = binomial(link = "logit")
    )
  )
  
  ## Extract the observed treatment indicator and estimate the marginal
  ## probability of receiving treatment
  A <- data.set$A
  p.A1 <- mean(A)
  
  # ----------------------------------------------------------------------------
  # Stabilized Inverse Probability Weights
  # ----------------------------------------------------------------------------
  
  ## Construct stabilized inverse probability weights using the estimated
  ## propensity scores and the marginal probability of treatment
  ipw.stabilized <- ifelse(
    A == 1,
    p.A1 / p.score,
    (1 - p.A1) / (1 - p.score)
  )
  
  # ==============================================================================
  # Marginal Structural Cox Model
  # ==============================================================================
  
  ## Initialize objects used to store the Cox model estimates
  ms.cox <- NULL
  cox.psi <- NULL
  cox.psi.se <- NULL
  
  if (constant.HR) {
    
    # --------------------------------------------------------------------------
    # Constant Hazard Ratio
    # --------------------------------------------------------------------------
    
    ## Fit the marginal structural Cox model using the stabilized IPWs
    ms.cox <- coxph(
      Surv(time, Delta) ~ A,
      data = data.set,
      weights = ipw.stabilized,
      robust = TRUE
    )
    
    ## Extract the estimated marginal log hazard ratio and its robust standard
    ## error
    cox.psi <- coef(ms.cox)["A"]
    cox.psi.se <- sqrt(diag(vcov(ms.cox))["A"])
    
    
  } else {
    
    # --------------------------------------------------------------------------
    # Time-Varying Hazard Ratio
    # --------------------------------------------------------------------------
    
    ## Allow the treatment effect to vary over time through the treatment-by-time
    ## interaction
    ms.cox <- coxph(
      Surv(time, Delta) ~ A * time,
      data = data.set,
      weights = ipw.stabilized,
      robust = TRUE
    )
    
    ## A single scalar treatment-effect estimate is not extracted under the
    ## time-varying specification
    cox.psi <- NULL
    cox.psi.se <- NULL
  }
  
  # ----------------------------------------------------------------------------
  # MS-Cox Counterfactual Survival Curves and RMST
  # ----------------------------------------------------------------------------
  
  ## Define the two treatment profiles used to obtain counterfactual survival
  ## curves under A = 0 and A = 1
  survival.curves.data <- data.frame(
    A = c(0, 1)
  )
  
  ## Calculate the RMST contrast and its standard error from the fitted
  ## marginal structural Cox model
  cox.rmst.list <- cox.rmst.list(
    ms.cox = ms.cox,
    t.grid = t.grid,
    cox.psi.se = cox.psi.se,
    survival.curves.data = survival.curves.data
  )
  
  ## Extract the RMST contrast and its standard error
  cox.rmst <- cox.rmst.list$cox.rmst
  cox.rmst.se <- cox.rmst.list$cox.rmst.se
  
  # ==============================================================================
  # IPW Case-Base Model
  # ==============================================================================
  
  ## Set the number of controls sampled per case in the case-base design
  fsh.ratio <- 100
  
  ## Initialize objects used to store the case-base model estimates
  cb.ipw <- NULL
  cb.psi <- NULL
  cb.psi.se <- NULL
  cb.var <- NULL
  
  if (constant.HR) {
    
    # --------------------------------------------------------------------------
    # Constant Hazard Ratio
    # --------------------------------------------------------------------------
    
    ## Fit the IPW case-base model using the stabilized inverse probability
    ## weights
    cb.ipw <- FSH(
      Delta ~ A,
      data = data.set,
      ratio = fsh.ratio,
      w = ipw.stabilized
    )
    
    ## Extract the estimated marginal log hazard ratio and its
    ## heteroskedasticity-consistent standard error
    cb.psi <- coef(cb.ipw)["A"]
    cb.var <- vcovHC(
      cb.ipw,
      type = "HC0"
    )
    cb.psi.se <- sqrt(diag(cb.var)["A"])
    
    
  } else {
    
    # --------------------------------------------------------------------------
    # Time-Varying Hazard Ratio
    # --------------------------------------------------------------------------
    
    ## Fit the IPW case-base model with a treatment-by-time interaction to
    ## allow the treatment effect to vary over follow-up
    cb.ipw <- FSH(
      Delta ~ A * pspline(time, df = 4),
      data = data.set,
      ratio = fsh.ratio,
      w = ipw.stabilized
    )
    
    ## A single scalar treatment-effect estimate and standard error are not
    ## extracted under the time-varying specification
    cb.psi <- NULL
    cb.psi.se <- NULL
    cb.var <- NULL
  }
  
  # ----------------------------------------------------------------------------
  # IPW Case-Base RMST
  # ----------------------------------------------------------------------------
  
  ## Calculate the RMST contrast and its standard error from the fitted IPW
  ## case-base model
  cb.rmst.list <- cb.rmst.list(
    cb.ipw,
    cb.var,
    t.grid,
    constant.HR
  )
  
  ## Extract the RMST contrast and its standard error
  cb.rmst <- cb.rmst.list$cb.rmst
  cb.rmst.se <- cb.rmst.list$cb.rmst.se
  
  # ==============================================================================
  # Store Simulation Results
  # ==============================================================================
  
  ## Return all treatment-effect and RMST estimates from the current simulation
  ## replicate, together with the replicate identifier
  output <- list(
    simulation.number = i,
    cox.psi = cox.psi,
    cox.psi.se = cox.psi.se,
    cb.psi = cb.psi,
    cb.psi.se = cb.psi.se,
    cox.rmst = cox.rmst,
    cox.rmst.se = cox.rmst.se,
    cb.rmst = cb.rmst,
    cb.rmst.se = cb.rmst.se
  )
  
  return(output)
}