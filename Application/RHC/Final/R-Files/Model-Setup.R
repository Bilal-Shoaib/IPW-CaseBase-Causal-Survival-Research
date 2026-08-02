# ==============================================================================
# Case-Base Model Setup
# ==============================================================================

## Set the case-base sampling ratio used across all case-base model fits.
## A ratio of 100 samples 100 background person-moments per observed case.
fsh.ratio <- 100

# ------------------------------------------------------------------------------
# Plotting Parameters
# ------------------------------------------------------------------------------

## Define colours used to distinguish the marginal structural Cox model
## and the case-base model in survival and hazard ratio plots
col.cox <- "red"
col.cb <- "blue"

## Define semi-transparent fill colours for confidence interval bands
fill.cox <- rgb(1,0,0,0.20)
fill.cb  <- rgb(0,0,1,0.20)

## Define semi-transparent colours for the exposed and control treatment groups
col.exposed <- rgb(0,0,1,0.5)
col.control <- rgb(1,0,0,0.5)

# ==============================================================================
# Counterfactual Survival Curve Setup
# ==============================================================================

## Define the two treatment levels used to generate counterfactual survival
## curves: A = 0 represents no RHC and A = 1 represents RHC
trt.types <- c(0, 1)

## Create a single-treatment profile for prediction under the treated condition
individual.profile.data <- data.frame(A = 1)

## Create a prediction dataset containing both treatment levels for estimation
## of counterfactual survival under each treatment condition
survival.curves.data <- data.frame(A = trt.types)

# ------------------------------------------------------------------------------
# Time Grid Definition
# ------------------------------------------------------------------------------

## Determine the maximum observed follow-up time in the analysis dataset
tau.max <- max(data.set$time, na.rm = TRUE)

## Restrict the primary survival-curve display to the first 150 days
## (approximately four months) to focus on the relevant follow-up period
t.max <- 150

## Define the time grid used for RMST integration, spanning the full observed
## follow-up period at one-day intervals
tau.grid <- seq(0, tau.max, by = 1)

## Define the time grid used for plotting and estimating survival curves over
## a restricted four-month follow-up period
t.grid <- seq(0, t.max, by = 1)