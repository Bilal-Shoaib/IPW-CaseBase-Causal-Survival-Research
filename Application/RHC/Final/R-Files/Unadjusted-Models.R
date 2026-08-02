# ==============================================================================
# Model 1: Unadjusted Marginal Structural Cox Model
# ==============================================================================

## Set seed to ensure reproducibility of any stochastic procedures used
## during model fitting
set.seed(s)

## Fit an unadjusted Cox proportional hazards model with treatment as the
## sole predictor. The resulting coefficient estimates the marginal
## log-hazard ratio comparing RHC with no RHC.
model.1 <- coxph(
  Surv(time, Delta) ~ A,
  data = data.set,
  robust = TRUE
)

## Extract the treatment coefficient and its standard error
beta.1 <- coef(model.1)["A"]
se.1 <- sqrt(diag(vcov(model.1))["A"])

## Transform the treatment coefficient from the log-hazard ratio scale to the
## hazard ratio scale and calculate the corresponding 95% Wald confidence
## interval
cat("HR-hat:", exp(beta.1), "\n")
cat("SE of HR-hat:", se.1, "\n")
cat("95% CI: (", exp(beta.1 + c(-1,1)*1.96*se.1), ")\n")

# ==============================================================================
# Model 2: Unadjusted Case-Base Model
# ==============================================================================

## Set seed to ensure reproducibility of the case-base sampling procedure
set.seed(s)

## Fit an unadjusted smooth hazard model using case-base sampling, with
## treatment as the sole predictor. The sampling ratio determines the number
## of background person-moments sampled for each observed case.
model.2 <- fitSmoothHazard(
  Delta ~ A,
  data = data.set,
  ratio = fsh.ratio
)

## Extract the treatment coefficient and calculate its robust standard error
## using the heteroskedasticity-consistent covariance estimator
beta.2 <- coef(model.2)["A"]
se.2 <- sqrt(diag(vcovHC(model.2))["A"])

## Transform the treatment coefficient to the hazard ratio scale and calculate
## the corresponding 95% Wald confidence interval
cat("HR-hat:", exp(beta.2), "\n")
cat("SE of HR-hat:", se.2, "\n")
cat("95% CI: (", exp(beta.2 + c(-1,1)*1.96*se.2), ")\n")