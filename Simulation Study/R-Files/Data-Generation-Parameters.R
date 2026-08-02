# ==============================================================================
# Simulation Study Parameters
# ==============================================================================

## General Simulation Parameters
seed <- 123

tau.max <- 4
t.grid <- seq(0, tau.max, by = 0.1)

# Covariate distribution parameters
r <- 2
l <- 1.5

# Treatment assignment model coefficients
treatment.predictor.coefs <- c(
  0.7,  -0.4,  0.5,  0.2,  0.6,
  -0.3,   0.1, -0.2,  0.4, -0.5
)

# Outcome model covariate coefficients
outcome.predictor.coefs <- c(
  0.6,  -0.3,  0.4,  0.1,  0.5,
  -0.2,  0.05, -0.1,  0.3, -0.4
)

# Treatment assignment model intercept
intercept <- 0.3

# Specify a constant marginal hazard ratio
constant.HR <- TRUE

# Weibull baseline hazard shape parameter
shape <- 1

# ==============================================================================
# Simulation Setting 1: Null Treatment Effect
# ==============================================================================

n1 <- 1000
beta.A1 <- 0
censoring.limit.1 <- 5

# ==============================================================================
# Simulation Setting 2: Strong Treatment Effect
# ==============================================================================

n2 <- 1000
beta.A2 <- 1
censoring.limit.2 <- 5

# ==============================================================================
# Simulation Setting 3: Moderate Treatment Effect
# ==============================================================================

n3 <- 1000
beta.A3 <- 0.5
censoring.limit.3 <- 5

# ==============================================================================
# Simulation Setting 4: Moderate Treatment Effect with Smaller Sample Size
# ==============================================================================

n4 <- 500
beta.A4 <- 0.5
censoring.limit.4 <- 5

# ==============================================================================
# Simulation Setting 5: Moderate Treatment Effect with Small Sample Size
# ==============================================================================

n5 <- 250
beta.A5 <- 0.5
censoring.limit.5 <- 5

# ==============================================================================
# Simulation Setting 6: Moderate Treatment Effect with Heavy Censoring
# ==============================================================================

n6 <- 1000
beta.A6 <- 0.5
censoring.limit.6 <- 3