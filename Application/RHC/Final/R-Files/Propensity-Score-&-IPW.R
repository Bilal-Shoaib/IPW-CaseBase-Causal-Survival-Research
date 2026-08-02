# ==============================================================================
# Propensity Score Estimation
# ==============================================================================

## Set seed to ensure reproducibility of any downstream procedures involving
## randomization or resampling
set.seed(s)

## Estimate the propensity score, defined as the probability of receiving
## treatment conditional on the observed baseline covariates
propensity.score <- fitted(
  glm(
    A ~ Gender +
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
    family = binomial(link = "logit")
  )
)

## Separate estimated propensity scores by treatment group to assess the
## degree of overlap in treatment assignment probabilities
propensity.score.0 <- propensity.score[A == 0]
propensity.score.1 <- propensity.score[A == 1]

# ------------------------------------------------------------------------------
# Propensity Score Overlap
# ------------------------------------------------------------------------------

## Plot the propensity score distributions for treated and untreated patients
## to visually assess the positivity/overlap condition
par(
  mar=c(2,2,2,0)
)

hist(
  propensity.score.0,
  col=col.control,
  breaks=seq(0,1,by=0.02),
  ylim=c(0,200),
  main="Propensity Score Overlap",
  xlab="Propensity Scores"
)

hist(
  propensity.score.1,
  col=col.exposed,
  breaks=seq(0,1,by=0.02),
  add=T
)

box()

legend(
  "topright",
  legend=c('A = 0', 'A = 1'),
  fill=c(col.control,col.exposed),
  border="black",
  bty="n"
)

# ==============================================================================
# Stabilized Inverse Probability Weights
# ==============================================================================

## Estimate the marginal probability of treatment, used to construct the
## numerator of the stabilized treatment weights
p.A1 <- mean(A)

## Calculate unstabilized inverse probability of treatment weights (IPWs)
## by taking the inverse of the estimated probability of the treatment
## actually received
ipw.unstabilized <- ifelse(
  A == 1,
  1 / propensity.score,
  1 / (1 - propensity.score)
)

## Calculate stabilized IPWs by incorporating the marginal probability of
## treatment in the numerator. Stabilization reduces the variability of the
## weights while preserving their role in creating a weighted pseudo-population
ipw.stabilized <- ifelse(
  A == 1,
  p.A1 / propensity.score,
  (1 - p.A1) / (1 - propensity.score)
)

# ------------------------------------------------------------------------------
# IPW Distribution
# ------------------------------------------------------------------------------

## Compare the distributions of the unstabilized and stabilized IPWs to
## assess the magnitude and variability of the resulting weights
par(
  mar=c(2,2,2,0)
)

hist(
  ipw.unstabilized,
  col=col.control,
  breaks=seq(0, max(ipw.unstabilized) + 0.5,by=0.1),
  ylim=c(0,1000),
  xlim=c(0, 10),
  main="Unstablized vs. Stabilized IPWs",
  ylab="Frequency",
  xlab="IPWs"
)

hist(
  ipw.stabilized,
  col=col.exposed,
  breaks=seq(0, max(ipw.stabilized) + 0.5,by=0.1),
  add=T
)

box()

legend(
  "topright",
  legend=c('Unstabilized', 'Stabilized'),
  fill=c(col.control,col.exposed),
  border="black", bty="n"
)