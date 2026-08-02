# ==============================================================================
# RMST Comparison: MS-Cox vs. IPW Case-Base
# ==============================================================================

# ------------------------------------------------------------------------------
# RMST Estimates and Standard Errors
# ------------------------------------------------------------------------------

## Extract the robust standard error of the treatment log-hazard coefficient
## from the marginal structural Cox model
model.5.psi.se <- sqrt(diag(vcov(model.5))["A"])

## Obtain the heteroskedasticity-consistent covariance matrix for the IPW
## case-base model. This covariance matrix is used in the numerical
## delta-method calculation of the RMST standard error.
model.8.var <- cb.var <- vcovHC(model.8, type = "HC0")


## Calculate the RMST contrast and its delta-method standard error for the
## marginal structural Cox model
cox.rmst.list <- cox.rmst.list(
  model.5, tau.grid, model.5.psi.se, survival.curves.data
)

## Calculate the RMST contrast and its delta-method standard error for the
## IPW case-base model
cb.rmst.list <- cb.rmst.list(
  model.8, model.8.var, tau.grid, constant.HR = TRUE
)

# ------------------------------------------------------------------------------
# Extract RMST Estimates and Standard Errors
# ------------------------------------------------------------------------------

## Extract the estimated RMST contrasts and their corresponding standard errors
## from both estimation approaches
cox.rmst <- cox.rmst.list$cox.rmst
cox.rmst.se <- cox.rmst.list$cox.rmst.se

cb.rmst <- cb.rmst.list$cb.rmst
cb.rmst.se <- cb.rmst.list$cb.rmst.se

# ------------------------------------------------------------------------------
# Pointwise 95% Confidence Intervals
# ------------------------------------------------------------------------------

## Construct pointwise 95% Wald confidence intervals for the MS-Cox RMST
## contrast over the range of restriction times
cox.lower <- cox.rmst - 1.96 * cox.rmst.se
cox.upper <- cox.rmst + 1.96 * cox.rmst.se

## Construct pointwise 95% Wald confidence intervals for the IPW case-base
## RMST contrast over the range of restriction times
cb.lower <- cb.rmst - 1.96 * cb.rmst.se
cb.upper <- cb.rmst + 1.96 * cb.rmst.se

# ==============================================================================
# RMST Comparison Plot
# ==============================================================================

## Determine a common y-axis range that accommodates the point estimates and
## confidence intervals from both estimation approaches
ylim <- range(
  c(
    cox.lower, cox.upper,
    cb.lower, cb.upper
  ),
  na.rm = TRUE
)

## Initialize the plotting region for comparing the estimated causal RMST
## contrasts over the full range of restriction times
plot(
  tau.grid, cox.rmst,
  type = "n",
  ylim = ylim,
  xlab = expression(tau),
  ylab = "RMST Difference",
  main = "Causal RMST up to Tau"
)

# ------------------------------------------------------------------------------
# 95% Confidence Bands
# ------------------------------------------------------------------------------

## Add the pointwise 95% confidence band for the MS-Cox RMST contrast
polygon(
  c(tau.grid, rev(tau.grid)),
  c(cox.lower, rev(cox.upper)),
  col = adjustcolor(col.cox, alpha.f = 0.2),
  border = NA
)

## Add the pointwise 95% confidence band for the IPW case-base RMST contrast
polygon(
  c(tau.grid, rev(tau.grid)),
  c(cb.lower, rev(cb.upper)),
  col = adjustcolor(col.cb, alpha.f = 0.2),
  border = NA
)

# ------------------------------------------------------------------------------
# RMST Contrast Estimates
# ------------------------------------------------------------------------------

## Overlay the estimated causal RMST contrast from the MS-Cox model
lines(tau.grid, cox.rmst, col = col.cox, lwd = 2)

## Overlay the estimated causal RMST contrast from the IPW case-base model
lines(tau.grid, cb.rmst, col = col.cb, lwd = 2)

## Add a reference line at an RMST difference of zero, corresponding to equal
## restricted mean survival under the two treatment conditions
abline(h = 0, lty = 3)

# ------------------------------------------------------------------------------
# Model Legend
# ------------------------------------------------------------------------------

## Identify the two estimation approaches in the comparison plot
legend(
  "bottomleft",
  legend = c("MS-Cox", "IPW Case-Base"),
  col = c(col.cox, col.cb),
  lwd = 2,
  bty = "n"
)