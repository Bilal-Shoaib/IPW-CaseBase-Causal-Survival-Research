# ==============================================================================
# Extract Marginal Estimates for All Simulation Settings
# ==============================================================================

## Moderate Treatment Effect
true.estimates.moderate <- intermediate.step(
  r, l, 
  intercept,
  outcome.predictor.coefs,
  beta.A3,
  censoring.limit.3,
  constant.HR
)

## Null Treatment Effect
true.estimates.null <- intermediate.step(
  r, l, 
  intercept,
  outcome.predictor.coefs,
  beta.A1,
  censoring.limit.1,
  constant.HR
)

## Strong Treatment Effect
true.estimates.strong <- intermediate.step(
  r, l, 
  intercept,
  outcome.predictor.coefs,
  beta.A2,
  censoring.limit.2,
  constant.HR
)