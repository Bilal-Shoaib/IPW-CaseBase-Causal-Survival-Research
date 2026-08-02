# ==============================================================================
# Model 5: Marginal Structural Cox Model
# ==============================================================================

## Set seed to ensure reproducibility of any stochastic procedures used
## during model fitting
set.seed(s)

## Fit the marginal structural Cox model using stabilized inverse probability
## weights (IPWs). The treatment-only model estimates the marginal treatment
## effect in the IPW-defined pseudo-population.
model.5 <- coxph(
  Surv(time, Delta) ~ A,
  data = data.set,
  weights = ipw.stabilized,
  robust = TRUE
)

## Extract the treatment coefficient and its robust standard error
beta.5 <- coef(model.5)["A"]
se.5 <- sqrt(diag(vcov(model.5))["A"])

## Transform the treatment coefficient to the hazard ratio scale and calculate
## the corresponding 95% Wald confidence interval
cat("HR-hat:", exp(beta.5), "\n")
cat("SE of HR-hat:", se.5, "\n")
cat("95% CI: (", exp(beta.5 + c(-1,1)*1.96*se.5), ")\n")

# ==============================================================================
# Model 6: IPW Case-Base Proportional Hazards Model
# ==============================================================================

## Set seed to ensure reproducibility of the case-base sampling procedure
set.seed(s)

## Fit the modified smooth hazard model using case-base sampling and stabilized
## IPWs. The custom FSH function incorporates the IPWs as observation-level
## weights while accounting for the case-base sampling mechanism.
model.6 <- FSH(
  Delta ~ A,
  data = data.set,
  ratio = fsh.ratio,
  w = ipw.stabilized
)

## Extract the treatment coefficient and calculate its heteroskedasticity-
## consistent standard error
beta.6 <- coef(model.6)["A"]
se.6 <- sqrt(diag(vcovHC(model.6))["A"])

## Transform the treatment coefficient to the hazard ratio scale and calculate
## the corresponding 95% Wald confidence interval
cat("HR-hat:", exp(beta.6), "\n")
cat("SE of HR-hat:", se.6, "\n")
cat("95% CI: (", exp(beta.6 + c(-1,1)*1.96*se.6), ")\n")