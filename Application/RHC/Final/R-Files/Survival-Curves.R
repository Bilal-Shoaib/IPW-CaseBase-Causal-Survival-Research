# ==============================================================================
# Counterfactual Survival Curves: MS-Cox and IPW Case-Base Models
# ==============================================================================

# ------------------------------------------------------------------------------
# MS-Cox Counterfactual Survival Curves
# ------------------------------------------------------------------------------

## Estimate counterfactual survival under both treatment levels using the
## marginal structural Cox model. A = 0 corresponds to no RHC, while A = 1
## corresponds to RHC.
sf.cox <- survfit(
  model.5,
  newdata = survival.curves.data
)

## Interpolate the estimated survival curve under no RHC onto the common
## time grid used for subsequent comparison and RMST calculation
S0.cox <- approx(
  sf.cox$time,
  sf.cox$surv[, 1],
  xout = tau.grid,
  method = "constant",
  rule = 2
)$y

## Interpolate the corresponding lower 95% confidence limit under no RHC
S0.lower.cox <- approx(
  sf.cox$time,
  sf.cox$lower[, 1],
  xout = tau.grid,
  method = "constant",
  rule = 2
)$y

## Interpolate the corresponding upper 95% confidence limit under no RHC
S0.upper.cox <- approx(
  sf.cox$time,
  sf.cox$upper[, 1],
  xout = tau.grid,
  method = "constant",
  rule = 2
)$y

## Interpolate the estimated survival curve under RHC onto the common time grid
S1.cox <- approx(
  sf.cox$time,
  sf.cox$surv[, 2],
  xout = tau.grid,
  method = "constant",
  rule = 2
)$y

## Interpolate the corresponding lower 95% confidence limit under RHC
S1.lower.cox <- approx(
  sf.cox$time,
  sf.cox$lower[, 2],
  xout = tau.grid,
  method = "constant",
  rule = 2
)$y

## Interpolate the corresponding upper 95% confidence limit under RHC
S1.upper.cox <- approx(
  sf.cox$time,
  sf.cox$upper[, 2],
  xout = tau.grid,
  method = "constant",
  rule = 2
)$y

# ==============================================================================
# IPW Case-Base Counterfactual Survival Curves
# ==============================================================================

## Estimate counterfactual survival under no RHC and RHC using the IPW
## case-base model over the common time grid
Sf.cb <- absoluteRisk(
  object = model.8,
  time = tau.grid,
  newdata = survival.curves.data,
  type = "survival"
)

## Extract the estimated survival curves for no RHC and RHC, respectively
S0.cb <- Sf.cb[, 2]
S1.cb <- Sf.cb[, 3]

# ------------------------------------------------------------------------------
# Bootstrap Confidence Intervals for the IPW Case-Base Curves
# ------------------------------------------------------------------------------

## Use nonparametric bootstrap resampling to obtain pointwise confidence
## intervals for the counterfactual survival curves estimated by the IPW
## case-base model
set.seed(s)
n.boot <- 200

## For each bootstrap sample, refit the IPW case-base model and obtain
## counterfactual survival estimates under both treatment levels
boot.surv <- replicate(n.boot, {
  idx <- sample(nrow(data.set), replace = TRUE)
  boot.data <- data.set[idx, ]
  
  boot.fit <- FSH(
    Delta ~ pspline(time, df = 2) * A,
    data = boot.data,
    ratio = 10,
    w = ipw.stabilized[idx]
  )
  
  ar <- absoluteRisk(
    boot.fit,
    time = tau.grid,
    newdata = survival.curves.data,
    type = "survival"
  )
  
  list(
    S0 = ar[, 2],
    S1 = ar[, 3]
  )
})

## Organize the bootstrap estimates into separate matrices for the no-RHC
## and RHC counterfactual survival curves
S0.boot <- sapply(1:n.boot, function(i) boot.surv[["S0", i]])
S1.boot <- sapply(1:n.boot, function(i) boot.surv[["S1", i]])

## Calculate the pointwise lower 95% bootstrap confidence limit for the
## no-RHC counterfactual survival curve
S0.lower.cb <- apply(
  S0.boot,
  1,
  quantile,
  probs = 0.025,
  na.rm = TRUE
)

## Calculate the pointwise upper 95% bootstrap confidence limit for the
## no-RHC counterfactual survival curve
S0.upper.cb <- apply(
  S0.boot,
  1,
  quantile,
  probs = 0.975,
  na.rm = TRUE
)

## Calculate the pointwise lower 95% bootstrap confidence limit for the
## RHC counterfactual survival curve
S1.lower.cb <- apply(
  S1.boot,
  1,
  quantile,
  probs = 0.025,
  na.rm = TRUE
)

## Calculate the pointwise upper 95% bootstrap confidence limit for the
## RHC counterfactual survival curve
S1.upper.cb <- apply(
  S1.boot,
  1,
  quantile,
  probs = 0.975,
  na.rm = TRUE
)

# ==============================================================================
# Head-to-Head Comparison of Counterfactual Survival Curves
# ==============================================================================

## Compare the MS-Cox and IPW case-base estimates separately under each
## treatment condition. The confidence bands are plotted first so that the
## survival curves remain visible on top of the shaded regions.

# ------------------------------------------------------------------------------
# Counterfactual Survival Under No RHC
# ------------------------------------------------------------------------------

## Initialize the plotting region for the counterfactual survival curves
## under no RHC (A = 0)
plot(
  tau.grid, S0.cox,
  type = "n",
  lwd = 2,
  col = col.cox,
  lty = 1,
  ylim = c(0, 1),
  xlab = "Time",
  ylab = "Survival Probability",
  main = "Counterfactual Survival Curves (No RHC)"
)

## Add the 95% confidence band for the MS-Cox survival estimate
polygon(
  c(tau.grid, rev(tau.grid)),
  c(S0.lower.cox, rev(S0.upper.cox)),
  col = fill.cox,
  border = NA
)

## Add the 95% bootstrap confidence band for the IPW case-base survival
## estimate
polygon(
  c(tau.grid, rev(tau.grid)),
  c(S0.lower.cb, rev(S0.upper.cb)),
  col = fill.cb,
  border = NA
)

## Overlay the estimated survival curves. The MS-Cox curve is a step function,
## whereas the case-base model produces a smooth survival curve.
lines(S0.cox, type = "s", lwd = 2, col = col.cox, lty = 1)
lines(S0.cb, type = "l", lwd = 2, col = col.cb, lty = 1)

## Identify the two estimation approaches in the comparison plot
legend(
  "topright",
  legend = c(
    "Semi-Parametric: MS-Cox",
    "Parametric: Case Base"
  ),
  col = c(
    col.cox, col.cb
  ),
  lty = c(1, 1),
  lwd = 2,
  bty = "n"
)

# ------------------------------------------------------------------------------
# Counterfactual Survival Under RHC
# ------------------------------------------------------------------------------

## Initialize the plotting region for the counterfactual survival curves
## under RHC (A = 1)
plot(
  tau.grid, S1.cox,
  type = "n",
  lwd = 2,
  col = col.cox,
  lty = 1,
  ylim = c(0, 1),
  xlab = "Time",
  ylab = "Survival Probability",
  main = "Counterfactual Survival Curves (RHC)"
)

## Add the 95% confidence band for the MS-Cox survival estimate
polygon(
  c(tau.grid, rev(tau.grid)),
  c(S1.lower.cox, rev(S1.upper.cox)),
  col = fill.cox,
  border = NA
)

## Add the 95% bootstrap confidence band for the IPW case-base survival
## estimate
polygon(
  c(tau.grid, rev(tau.grid)),
  c(S1.lower.cb, rev(S1.upper.cb)),
  col = fill.cb,
  border = NA
)

## Overlay the estimated survival curves under RHC
lines(S1.cox, type = "s", lwd = 2, col = col.cox, lty = 1)
lines(S1.cb, type = "l", lwd = 2, col = col.cb, lty = 1)

## Identify the two estimation approaches in the comparison plot
legend(
  "topright",
  legend = c(
    "Semi-Parametric: MS-Cox",
    "Parametric: Case Base"
  ),
  col = c(
    col.cox, col.cb
  ),
  lty = c(1, 1),
  lwd = 2,
  bty = "n"
)