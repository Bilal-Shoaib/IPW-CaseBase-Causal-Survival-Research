# ==============================================================================
# Restricted Mean Survival Time (RMST)
# ==============================================================================

## Numerically integrate a survival curve over the specified follow-up period
## using the trapezoidal rule. The cumulative integral at each time point
## represents the RMST up to that point.
rmst.fn <- function(time, surv) {
  dt <- diff(time)
  mids <- (head(surv, -1) + tail(surv, -1)) / 2
  c(0, cumsum(mids * dt))
}

# ==============================================================================
# MS-Cox RMST and Standard Error
# ==============================================================================

## Calculate the counterfactual RMST contrast and its standard error from the
## marginal structural Cox model. The RMST contrast is defined as the
## difference in restricted mean survival between RHC and no RHC.
cox.rmst.list <- function(
    ms.cox, t.grid, cox.psi.se, survival.curves.data 
) {
  
  # ----------------------------------------------------------------------------
  # Counterfactual Survival Curves
  # ----------------------------------------------------------------------------
  
  ## Estimate the counterfactual survival curves under both treatment levels
  ## from the fitted marginal structural Cox model
  sf <- survfit(
    ms.cox,
    newdata = survival.curves.data
  )
  
  ## Interpolate the counterfactual survival curve under no RHC onto the
  ## specified time grid
  S0 <- approx(
    sf$time,
    sf$surv[, 1],
    t.grid,
    method = "constant",
    rule = 2
  )$y
  
  ## Interpolate the counterfactual survival curve under RHC onto the
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
  ## survival under RHC and no RHC
  cox.rmst <- rmst.fn(t.grid, S1) - rmst.fn(t.grid, S0)
  
  # ----------------------------------------------------------------------------
  # Delta-Method Standard Error
  # ----------------------------------------------------------------------------
  
  ## Construct the gradient of the RMST functional with respect to the
  ## treatment log-hazard coefficient. For the Cox model, the derivative
  ## involves the integral of S1(t) log{S1(t)} over the restriction period.
  integrand <- ifelse(S1 > 0, S1 * log(S1), 0)
  g.cox <- rmst.fn(t.grid, integrand)
  
  ## Propagate uncertainty from the estimated Cox treatment coefficient to the
  ## RMST contrast using the delta method
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
## model. The survival curves are reconstructed from the fitted model and the
## variance of the RMST contrast is obtained using a numerical delta method.
cb.rmst.list <- function(cb.ipw, cb.var, t.grid, constant.HR) {
  
  ## Extract the estimated model coefficients and the corresponding model
  ## terms used to construct the design matrix
  beta.hat <- coef(cb.ipw)
  term.form <- delete.response(terms(cb.ipw))
  
  # ----------------------------------------------------------------------------
  # Counterfactual Survival Function
  # ----------------------------------------------------------------------------
  
  ## Construct the counterfactual survival function under treatment level a
  ## directly from the case-base model's design matrix
  cb.surv <- function(beta, a, grid) {
    
    ## Create prediction data over the specified time grid under treatment
    ## level a
    nd <- data.frame(
      time = grid,
      A = a,
      offset = 0
    )
    
    ## Construct the model matrix using the same terms as the fitted model and
    ## align its columns with the estimated coefficient vector
    X <- model.matrix(term.form, data = nd)
    X <- X[, names(beta), drop = FALSE]
    
    ## Convert the linear predictor to the estimated hazard function
    lambda <- exp(as.vector(X %*% beta))
    
    ## Numerically integrate the estimated hazard using the trapezoidal rule
    ## to obtain the cumulative hazard
    Lambda <- c(
      0,
      cumsum(
        diff(grid) * (head(lambda, -1) + tail(lambda, -1)) / 2
      )
    )
    
    ## Convert cumulative hazard to survival probability
    exp(-Lambda)
  }
  
  # ----------------------------------------------------------------------------
  # RMST Contrast
  # ----------------------------------------------------------------------------
  
  ## Define the RMST contrast as a function of the model coefficient vector.
  ## This function is subsequently used to obtain the numerical gradient
  ## required for the delta-method standard error.
  rmst.diff <- function(beta) {
    rmst.fn(t.grid, cb.surv(beta, 1, t.grid)) -
      rmst.fn(t.grid, cb.surv(beta, 0, t.grid))
  }
  
  ## Evaluate the RMST contrast at the estimated model coefficients
  cb.rmst <- rmst.diff(beta.hat)
  
  # ----------------------------------------------------------------------------
  # Numerical Gradient
  # ----------------------------------------------------------------------------
  
  ## Calculate the numerical gradient of the RMST contrast with respect to
  ## each model coefficient using a central finite-difference approximation
  p <- length(beta.hat)
  m <- length(t.grid)
  
  grad <- matrix(NA_real_, nrow = m, ncol = p)
  colnames(grad) <- names(beta.hat)
  
  ## Step size used for the finite-difference approximation
  eps <- 1e-6
  
  ## Perturb each coefficient in turn and approximate the corresponding
  ## derivative of the RMST contrast
  for (j in seq_len(p)) {
    bp <- bm <- beta.hat
    bp[j] <- bp[j] + eps
    bm[j] <- bm[j] - eps
    
    grad[, j] <- (rmst.diff(bp) - rmst.diff(bm)) / (2 * eps)
  }
  
  # ----------------------------------------------------------------------------
  # Delta-Method Standard Error
  # ----------------------------------------------------------------------------
  
  ## Align the model covariance matrix with the coefficient ordering used by
  ## the numerical gradient before applying the delta-method variance formula
  cb.var <- cb.var[colnames(grad), colnames(grad), drop = FALSE]
  
  ## Calculate the RMST standard error at each restriction time using
  ## Var{g(beta)} = grad' Var(beta) grad
  cb.rmst.se <- vapply(
    seq_len(m),
    function(i) {
      g <- grad[i, ]
      sqrt(drop(t(g) %*% cb.var %*% g))
    },
    numeric(1)
  )
  
  return(list(
    cb.rmst = cb.rmst,
    cb.rmst.se = cb.rmst.se
  ))
}