# ==============================================================================
# Model 3: Covariate-Adjusted Cox Model
# ==============================================================================

## Set seed to ensure reproducibility of any stochastic procedures used
## during model fitting
set.seed(s)

## Fit a covariate-adjusted Cox proportional hazards model. Treatment is
## adjusted for the full set of observed baseline demographic, socioeconomic,
## comorbidity, and disease-severity covariates.
model.3 <- coxph(
  Surv(time, Delta) ~ A +
    Gender +
    Age +
    Weight +
    Race +
    Education +
    Income +
    InsuranceClass +
    CHFhx +
    CARDIOhx +
    DEMENThx +
    PSYCHhx +
    CHRPULhx +
    RENALhx +
    LIVERhx +
    GIBLEDhx +
    MALIGhx +
    IMMUNhx +
    TRANShx +
    AMIhx + 
    Cancer + 
    Category + 
    Severity,
  data = data.set,
  robust = TRUE
)

## Extract the treatment coefficient and its robust standard error
beta.3 <- coef(model.3)["A"]
se.3 <- sqrt(diag(vcov(model.3))["A"])

## Transform the treatment coefficient to the hazard ratio scale and calculate
## the corresponding 95% Wald confidence interval
cat("HR-hat:", exp(beta.3), "\n")
cat("SE of HR-hat:", se.3, "\n")
cat("95% CI: (", exp(beta.3 + c(-1,1)*1.96*se.3), ")\n")

# ==============================================================================
# Model 4: Covariate-Adjusted Case-Base Model
# ==============================================================================

## Set seed to ensure reproducibility of the case-base sampling procedure
set.seed(s)

## Fit a covariate-adjusted smooth hazard model using case-base sampling.
## The model uses the same treatment and baseline covariates as the adjusted
## Cox model to provide a directly comparable alternative estimator.
model.4 <- fitSmoothHazard(
  Delta ~ A +
    Gender +
    Age +
    Weight +
    Race +
    Education +
    Income +
    InsuranceClass +
    CHFhx +
    CARDIOhx +
    DEMENThx +
    PSYCHhx +
    CHRPULhx +
    RENALhx +
    LIVERhx +
    GIBLEDhx +
    MALIGhx +
    IMMUNhx +
    TRANShx +
    AMIhx + 
    Cancer + 
    Category + 
    Severity,
  data = data.set,
  ratio = fsh.ratio
)

## Extract the treatment coefficient and calculate its heteroskedasticity-
## consistent standard error
beta.4 <- coef(model.4)["A"]
se.4 <- sqrt(diag(vcovHC(model.4))["A"])

## Transform the treatment coefficient to the hazard ratio scale and calculate
## the corresponding 95% Wald confidence interval
cat("HR-hat:", exp(beta.4), "\n")
cat("SE of HR-hat:", se.4, "\n")
cat("95% CI: (", exp(beta.4 + c(-1,1)*1.96*se.4), ")\n")