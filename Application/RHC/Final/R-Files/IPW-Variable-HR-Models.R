# ==============================================================================
# Model 7: Unadjusted Case-Base Model with Time-Varying Hazard Ratio
# ==============================================================================

## Set seed to ensure reproducibility of the case-base sampling procedure
set.seed(s)

## Fit an unadjusted case-base smooth hazard model allowing the treatment
## effect to vary over time. The spline term models the baseline time effect,
## while its interaction with treatment permits a non-constant hazard ratio.
model.7 <- FSH(
  Delta ~ pspline(time, df = 2) * A,
  data = data.set,
  ratio = fsh.ratio
)

# ------------------------------------------------------------------------------
# Survival Curve Prediction Setup
# ------------------------------------------------------------------------------

## Define the two treatment levels used to generate counterfactual survival
## curves under no RHC (A = 0) and RHC (A = 1)
trt.types <- c(0, 1)

## Define the treatment profile used for individual-level prediction
individual.profile.data <- data.frame(A = 1)

## Create a prediction dataset containing both treatment levels for estimation
## of the corresponding counterfactual survival curves
survival.curves.data <- data.frame(A = trt.types)

## Determine the maximum observed follow-up time and define the restricted
## four-month follow-up period used for focused visualization
tau.max <- max(data.set$time, na.rm = TRUE)
t.max <- 150

## Define the full follow-up grid and the restricted four-month grid
tau.grid <- seq(0, tau.max, by = 1)
t.grid <- seq(0, t.max, by = 1)

# ------------------------------------------------------------------------------
# Time-Varying Hazard Ratio: Full Follow-Up
# ------------------------------------------------------------------------------

## Display the estimated treatment hazard ratio as a function of follow-up
## time, including its 95% confidence interval
plot(
  model.7,
  type = "hr",
  newdata = expand.grid(
    time = tau.grid,
    A = trt.types
  ),
  var = "A",
  xvar = "time",
  ci = TRUE,
  ylim = c(0,2),
  xlab="Time",
  ylab="Hazard Ratio",
  main="Model 7: Variable Hazard Ratio"
)

## Add a reference line corresponding to a hazard ratio of 1, representing
## no instantaneous treatment effect
abline(h = 1, col = "red", lty = 2)

# ------------------------------------------------------------------------------
# Time-Varying Hazard Ratio: First Four Months
# ------------------------------------------------------------------------------

## Restrict the hazard-ratio plot to the first 150 days (approximately four
## months) to provide a more detailed view of the early follow-up period
plot(
  model.7,
  type = "hr",
  newdata = expand.grid(
    time = t.grid,
    A = trt.types
  ),
  var = "A",
  xvar = "time",
  ci = TRUE,
  ylim = c(1,1.5),
  xlab="Time",
  ylab="Hazard Ratio",
  main="Model 7: Variable Hazard Ratio - First 4 Months"
)

## Add the no-effect reference line
abline(h = 1, col = "red", lty = 2)

# ==============================================================================
# Model 8: IPW Case-Base Model with Time-Varying Hazard Ratio
# ==============================================================================

## Set seed to ensure reproducibility of the case-base sampling procedure
set.seed(s)

## Fit the IPW case-base smooth hazard model, allowing the treatment effect
## to vary over time. Stabilized IPWs are incorporated through the modified
## FSH function to estimate the treatment effect in the weighted
## pseudo-population.
model.8 <- FSH(
  Delta ~ pspline(time, df = 2) * A,
  data = data.set,
  ratio = fsh.ratio,
  w = ipw.stabilized
)

# ------------------------------------------------------------------------------
# Time-Varying Hazard Ratio: Full Follow-Up
# ------------------------------------------------------------------------------

## Display the estimated time-varying treatment hazard ratio and its 95%
## confidence interval over the full follow-up period
plot(
  model.8,
  type = "hr",
  newdata = expand.grid(
    time = tau.grid,
    A = trt.types
  ),
  var = "A",
  xvar = "time",
  ci = TRUE,
  ylim = c(0,2),
  xlab="Time",
  ylab="Hazard Ratio",
  main="Model 8: Variable Hazard Ratio"
)

## Add the no-effect reference line
abline(h = 1, col = "red", lty = 2)

# ------------------------------------------------------------------------------
# Time-Varying Hazard Ratio: First Four Months
# ------------------------------------------------------------------------------

## Restrict the hazard-ratio plot to the first 150 days to provide a more
## detailed visualization of the early treatment effect
plot(
  model.8,
  type = "hr",
  newdata = expand.grid(
    time = t.grid,
    A = trt.types
  ),
  var = "A",
  xvar = "time",
  ci = TRUE,
  ylim = c(1,1.5),
  xlab="Time",
  ylab="Hazard Ratio",
  main="Model 8: Variable Hazard Ratio - First 4 Months"
)

## Add the no-effect reference line
abline(h = 1, col = "red", lty = 2)

# ==============================================================================
# Head-to-Head Comparison: Models 7 and 8
# ==============================================================================

## Extract the estimated time-varying treatment hazard ratios and their 95%
## confidence intervals from the unadjusted and IPW-adjusted case-base models.
## The resulting estimates are used to construct a common comparison plot.
get.hr <- function(model) {
  invisible(pdf(NULL))
  
  out <- plot(
    model,
    type = "hr",
    newdata = data.frame(
      time = tau.grid,
      A = 0
    ),
    var = "A",
    xvar = "time",
    ci = TRUE
  )
  
  invisible(dev.off())
  
  out
}

## Extract the time-varying hazard-ratio estimates from Models 7 and 8
hr.7 <- get.hr(model.7)
hr.8 <- get.hr(model.8)

# ------------------------------------------------------------------------------
# Combined Hazard-Ratio Plot
# ------------------------------------------------------------------------------

## Initialize an empty plotting region using the time and hazard-ratio ranges
## from the extracted estimates. The curves and confidence intervals for both
## models are added to this common plotting region below.
plot(
  hr.7$time,
  hr.7$hazard_ratio,
  type = "n",
  ylim = c(0, 2),
  xlab = "Time",
  ylab = "Hazard Ratio",
  main = "Head-to-head Comparison of Variable HR Models"
)

# ------------------------------------------------------------------------------
# 95% Confidence Intervals
# ------------------------------------------------------------------------------

## Add the 95% confidence interval band for the unadjusted case-base model
polygon(
  c(hr.7$time, rev(hr.7$time)),
  c(hr.7$lowerbound, rev(hr.7$upperbound)),
  col = adjustcolor("red", alpha.f = 0.20),
  border = NA
)

## Add the 95% confidence interval band for the IPW-adjusted case-base model
polygon(
  c(hr.8$time, rev(hr.8$time)),
  c(hr.8$lowerbound, rev(hr.8$upperbound)),
  col = adjustcolor("blue", alpha.f = 0.20),
  border = NA
)

# ------------------------------------------------------------------------------
# Estimated Hazard-Ratio Curves
# ------------------------------------------------------------------------------

## Overlay the estimated time-varying hazard-ratio curve from the unadjusted
## case-base model
lines(
  hr.7$time,
  hr.7$hazard_ratio,
  col = "red",
  lwd = 2
)

## Overlay the estimated time-varying hazard-ratio curve from the IPW-adjusted
## case-base model
lines(
  hr.8$time,
  hr.8$hazard_ratio,
  col = "blue",
  lwd = 2
)

# ------------------------------------------------------------------------------
# Reference Line and Legend
# ------------------------------------------------------------------------------

## Add a horizontal reference line at HR = 1, corresponding to no difference
## in the instantaneous hazard between treatment groups
abline(
  h = 1,
  col = "black",
  lty = 2
)

## Identify the two estimation approaches shown in the combined plot
legend(
  "bottomleft",
  legend = c("Unadjusted", "IPW-Adjusted"),
  col = c("red", "blue"),
  lwd = 2,
  bty = "n"
)