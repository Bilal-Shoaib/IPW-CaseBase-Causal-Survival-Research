# ==============================================================================
# Restricted Mean Survival Time (RMST)
# ==============================================================================

## Numerically integrate a survival curve over time using the trapezoidal rule.
## The returned vector gives the cumulative RMST evaluated at each time point.
rmst.fn <- function(time, surv) {
  dt <- diff(time)
  mids <- (head(surv, -1) + tail(surv, -1)) / 2
  c(0, cumsum(mids * dt))
}

# ==============================================================================
# MS-Cox RMST and Standard Error
# ==============================================================================

## Calculate the RMST contrast and its standard error from the marginal
## structural Cox model. The contrast compares restricted mean survival under
## treatment (A = 1) with restricted mean survival under control (A = 0).
cox.rmst.list <- function(
    ms.cox, t.grid, cox.psi.se, survival.curves.data 
) {
  
  ## Estimate the counterfactual survival curves under A = 0 and A = 1
  sf <- survfit(
    ms.cox,
    newdata = survival.curves.data
  )
  
  ## Interpolate the counterfactual survival curve under A = 0 onto the
  ## specified time grid
  S0 <- approx(
    sf$time,
    sf$surv[, 1],
    t.grid,
    method = "constant",
    rule = 2
  )$y
  
  ## Interpolate the counterfactual survival curve under A = 1 onto the
  ## specified time grid
  S1 <- approx(
    sf$time,
    sf$surv[, 2],
    t.grid,
    method = "constant",
    rule = 2
  )$y
  
  # ----------------------------------------------------------------------------
  # RMST Contrast
  # ----------------------------------------------------------------------------
  
  ## Calculate the RMST contrast as the difference between the restricted mean
  ## survival under treatment and control
  cox.rmst <- rmst.fn(t.grid, S1) - rmst.fn(t.grid, S0)
  
  # ----------------------------------------------------------------------------
  # Delta-Method Standard Error
  # ----------------------------------------------------------------------------
  
  ## Construct the gradient of the RMST functional with respect to the Cox
  ## treatment coefficient. The gradient involves the integral of
  ## S1(t) * log{S1(t)} over the restriction period.
  integrand <- ifelse(
    S1 > 0,
    S1 * log(S1),
    0
  )
  
  g.cox <- rmst.fn(t.grid, integrand)
  
  ## Propagate uncertainty in the estimated treatment coefficient to the RMST
  ## contrast using the delta method
  cox.rmst.se <- abs(g.cox) * cox.psi.se
  
  return(list(
    cox.rmst = cox.rmst,
    cox.rmst.se = cox.rmst.se
  ))
}

# ==============================================================================
# IPW Case-Base RMST and Standard Error
# ==============================================================================

## Calculate the RMST contrast and its standard error from the IPW case-base
## model. The RMST standard error is obtained using a numerical delta method.
cb.rmst.list <- function(
    cb.ipw, cb.var, t.grid, constant.HR
) {
  
  ## Extract the fitted model coefficients and the model terms used to construct
  ## the corresponding design matrix
  beta.hat <- coef(cb.ipw)
  term.form <- delete.response(terms(cb.ipw))
  
  # ----------------------------------------------------------------------------
  # Counterfactual Survival Function
  # ----------------------------------------------------------------------------
  
  ## Construct the counterfactual survival function under treatment level a
  ## directly from the fitted case-base model
  cb.surv <- function(beta, a, grid) {
    
    ## Construct prediction data over the specified time grid for treatment
    ## level a
    nd <- data.frame(
      time = grid,
      A = a,
      offset = 0
    )
    
    ## Generate the model matrix using the fitted model's own terms and align
    ## its columns with the estimated coefficient vector
    X <- model.matrix(term.form, data = nd)
    X <- X[, names(beta), drop = FALSE]
    
    ## Obtain the estimated hazard function from the model's linear predictor
    lambda <- exp(as.vector(X %*% beta))
    
    ## Numerically integrate the hazard function using the trapezoidal rule to
    ## obtain the cumulative hazard
    Lambda <- c(
      0,
      cumsum(
        diff(grid) *
          (head(lambda, -1) + tail(lambda, -1)) / 2
      )
    )
    
    ## Convert the cumulative hazard to the corresponding survival function
    exp(-Lambda)
  }
  
  # ----------------------------------------------------------------------------
  # RMST Contrast
  # ----------------------------------------------------------------------------
  
  ## Define the RMST contrast as a function of the model coefficient vector.
  ## This function is evaluated at the fitted coefficients and is also used
  ## when calculating the numerical gradient.
  rmst.diff <- function(beta) {
    rmst.fn(t.grid, cb.surv(beta, 1, t.grid)) -
      rmst.fn(t.grid, cb.surv(beta, 0, t.grid))
  }
  
  ## Evaluate the RMST contrast at the fitted model coefficients
  cb.rmst <- rmst.diff(beta.hat)
  
  # ----------------------------------------------------------------------------
  # Numerical Gradient
  # ----------------------------------------------------------------------------
  
  ## Calculate the numerical gradient of the RMST contrast with respect to each
  ## fitted model coefficient using a central finite-difference approximation
  p <- length(beta.hat)
  m <- length(t.grid)
  
  grad <- matrix(
    NA_real_,
    nrow = m,
    ncol = p
  )
  
  colnames(grad) <- names(beta.hat)
  
  ## Finite-difference step size
  eps <- 1e-6
  
  ## Perturb each coefficient in turn to approximate its contribution to the
  ## RMST contrast
  for (j in seq_len(p)) {
    bp <- bm <- beta.hat
    
    bp[j] <- bp[j] + eps
    bm[j] <- bm[j] - eps
    
    grad[, j] <- (
      rmst.diff(bp) - rmst.diff(bm)
    ) / (2 * eps)
  }
  
  # ----------------------------------------------------------------------------
  # Delta-Method Standard Error
  # ----------------------------------------------------------------------------
  
  ## Align the covariance matrix with the coefficient ordering used by the
  ## numerical gradient before applying the delta-method variance formula
  cb.var <- cb.var[
    colnames(grad),
    colnames(grad),
    drop = FALSE
  ]
  
  ## Calculate the RMST standard error at each restriction time using the
  ## delta-method variance g' V g
  cb.rmst.se <- vapply(
    seq_len(m),
    function(i) {
      g <- grad[i, ]
      
      sqrt(
        drop(
          t(g) %*% cb.var %*% g
        )
      )
    },
    numeric(1)
  )
  
  return(list(
    cb.rmst = cb.rmst,
    cb.rmst.se = cb.rmst.se
  ))
}