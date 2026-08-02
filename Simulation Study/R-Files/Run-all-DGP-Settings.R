# ==============================================================================
# Run Simulation Study: All Settings
# ==============================================================================

## Setting 1: Null Treatment Effect
run.complete.ss(
  n1, beta.A1,
  censoring.limit.1,
  true.estimates.null,
  generate.plots = TRUE
)

## Setting 2: Strong Treatment Effect
run.complete.ss(
  n2, beta.A2,
  censoring.limit.2,
  true.estimates.strong,
  generate.plots = TRUE
)

## Setting 3: Moderate Treatment Effect
run.complete.ss(
  n3, beta.A3,
  censoring.limit.3,
  true.estimates.moderate,
  generate.plots = TRUE
)

## Setting 4: Moderate Treatment Effect with Smaller Sample Size
run.complete.ss(
  n4, beta.A4,
  censoring.limit.4,
  true.estimates.moderate
)

## Setting 5: Moderate Treatment Effect with Small Sample Size
run.complete.ss(
  n5, beta.A5,
  censoring.limit.5,
  true.estimates.moderate
)

## Setting 6: Moderate Treatment Effect with Heavy Censoring
run.complete.ss(
  n6, beta.A6,
  censoring.limit.6,
  true.estimates.moderate
)