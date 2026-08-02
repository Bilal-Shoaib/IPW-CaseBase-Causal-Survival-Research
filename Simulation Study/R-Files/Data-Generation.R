# ==============================================================================
# Data-Generating Function
# ==============================================================================

## Generate a simulated survival dataset with baseline confounding, a
## treatment effect, and independent random censoring.
generate.data <- function(
    n,
    r,
    l,
    treatment.predictor.coefs,
    intercept,
    outcome.predictor.coefs,
    beta.A,
    censoring.limit,
    constant.HR = TRUE,
    shape = 1 
) {
  
  # ----------------------------------------------------------------------------
  # Baseline Covariates
  # ----------------------------------------------------------------------------
  
  ## Generate 10 baseline covariates from a range of continuous and discrete
  ## distributions. These covariates are subsequently used to generate both
  ## treatment assignment and the event-time distribution.
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  X3 <- runif(n)
  X4 <- runif(n)
  X5 <- rlnorm(n)
  X6 <- rlnorm(n)
  X7 <- rexp(n, rate = r)
  X8 <- rbinom(n, 1, 0.6)
  X9 <- rbinom(n, 1, 0.4)
  X10 <- rpois(n, lambda = l)
  
  ## Combine the baseline covariates into a single design matrix for the
  ## treatment and outcome models
  X <- cbind(
    X1, X2, X3, X4, X5,
    X6, X7, X8, X9, X10
  )
  
  # ----------------------------------------------------------------------------
  # Input Validation
  # ----------------------------------------------------------------------------
  
  ## Ensure that the treatment and outcome coefficient vectors contain one
  ## coefficient for each simulated baseline covariate
  stopifnot(length(treatment.predictor.coefs) == ncol(X))
  stopifnot(length(outcome.predictor.coefs) == ncol(X))
  
  ## Ensure that all supplied model coefficients are finite
  stopifnot(all(is.finite(treatment.predictor.coefs)))
  stopifnot(all(is.finite(outcome.predictor.coefs)))
  
  # ----------------------------------------------------------------------------
  # Treatment Assignment
  # ----------------------------------------------------------------------------
  
  ## Generate treatment assignment from a logistic model in which all baseline
  ## covariates contribute to treatment assignment, thereby inducing
  ## confounding between treatment and the outcome
  weighted.covariates <- as.vector(
    X %*% treatment.predictor.coefs
  )
  
  linear.predictor.A <- intercept + weighted.covariates
  
  ## Convert the treatment linear predictor to a propensity score and generate
  ## the binary treatment indicator
  A <- rbinom(
    n,
    1,
    plogis(linear.predictor.A)
  )
  
  # ----------------------------------------------------------------------------
  # Event-Time Model and Censoring
  # ----------------------------------------------------------------------------
  
  ## Construct the covariate component of the event-time linear predictor.
  ## The treatment coefficient beta.A represents the true causal treatment
  ## effect in the constant-hazard-ratio setting.
  weighted.covariates <- as.vector(
    X %*% outcome.predictor.coefs
  )
  
  ## Generate independent random censoring times from a uniform distribution
  ## bounded by the specified censoring limit
  C <- runif(
    n,
    0,
    censoring.limit
  )
  
  ## Set the baseline event-time linear predictor. The intercept controls the
  ## baseline event-time distribution, while the covariates induce prognostic
  ## heterogeneity across individuals.
  linear.predictor.Y <- -1 + weighted.covariates
  
  T.star <- numeric(n)
  
  # ----------------------------------------------------------------------------
  # Constant Hazard-Ratio Outcome Model
  # ----------------------------------------------------------------------------
  
  if (constant.HR) {
    
    ## Add the treatment effect to the event-time linear predictor. Under this
    ## specification, beta.A determines the constant log hazard ratio
    ## associated with treatment.
    linear.predictor.Y <- linear.predictor.Y + (beta.A * A)
    
    ## Generate event times from a Weibull distribution. The scale parameter
    ## depends on the individual-specific linear predictor, allowing both
    ## treatment and baseline covariates to affect event-time risk.
    T.star <- rweibull(
      n,
      shape = shape,
      scale = exp(-linear.predictor.Y / shape)
    )
    
    # --------------------------------------------------------------------------
    # Time-Varying Hazard-Ratio Outcome Model
    # --------------------------------------------------------------------------
    
  } else {
    
    ## Placeholder for the time-varying treatment-effect specification
    ## (implemented separately when constant.HR = FALSE)
    linear.predictor.Y <- NULL
    
    T.star <- NULL
  }
  
  # ----------------------------------------------------------------------------
  # Observed Survival Data
  # ----------------------------------------------------------------------------
  
  ## Determine the observed follow-up time as the minimum of the event and
  ## censoring times
  time <- pmin(
    T.star,
    C
  )
  
  ## Define the event indicator: 1 if the event occurs before or at the
  ## censoring time, and 0 otherwise
  Delta <- as.integer(
    T.star <= C
  )
  
  # ----------------------------------------------------------------------------
  # Return Simulated Dataset
  # ----------------------------------------------------------------------------
  
  ## Combine treatment, observed survival outcomes, and baseline covariates
  ## into the final simulated analysis dataset
  data_set <- data.frame(
    time = as.numeric(time),
    Delta = Delta,
    A = A,
    X
  )
  
  return(data_set)
}