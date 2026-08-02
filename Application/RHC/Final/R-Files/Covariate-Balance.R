# ==============================================================================
# Descriptive Statistics: Table One
# ==============================================================================

## Specify the baseline covariates to be included in the descriptive table
covariates = c(
  "Gender",
  "Age",
  "Weight",
  "Race",
  "Education",
  "Income",
  "InsuranceClass",
  
  "CHFhx",
  "CARDIOhx",
  "DEMENThx",
  "PSYCHhx",
  "CHRPULhx",
  "RENALhx",
  "LIVERhx",
  "GIBLEDhx",
  "MALIGhx",
  "IMMUNhx",
  "TRANShx",
  "AMIhx",
  
  "Cancer",
  "Category",
  "Severity"
)

## Generate the unweighted baseline characteristics table, stratified by
## treatment group. Statistical significance tests are omitted because the
## primary objective is descriptive comparison of baseline characteristics.
table.one.original <- CreateTableOne(
  vars = covariates,
  strata = "A",
  data = data.set,
  test = FALSE
)

## Display standardized mean differences (SMDs) to quantify baseline
## differences between the treatment groups
print(table.one.original, smd = TRUE)

# ==============================================================================
# Weighted Descriptive Statistics: Table Two
# ==============================================================================

## Construct a survey design object using the stabilized IPWs. Setting ids = ~0
## specifies that observations are treated as independent sampling units.
data.set.ipw <- svydesign(
  ids = ~ 0,
  data = data.set,
  weights = ipw.stabilized
)

## Generate the weighted baseline characteristics table using the IPW-defined
## pseudo-population, stratified by treatment group
table.one.ipw <- svyCreateTableOne(
  vars = covariates,
  strata = "A",
  data = data.set.ipw,
  test = FALSE
)

## Display standardized mean differences for the weighted covariate
## distributions to assess residual imbalance after IPW
print(table.one.ipw, smd = TRUE)

# ==============================================================================
# Covariate Balance Before and After IPW
# ==============================================================================

## Assess covariate balance before and after application of the stabilized
## inverse probability weights. The unweighted and weighted SMDs are retained
## to facilitate direct comparison of balance across treatment groups.
balance <- bal.tab(
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
  weights = ipw.stabilized,
  method = "weighting",
  estimand = "Average Treatment Effect",
  un = TRUE
)

## Visualize the standardized mean differences before and after IPW.
## An absolute SMD below 0.1 is used as the criterion for acceptable
## covariate balance.
love.plot(
  balance,
  stats = "mean.diffs",
  threshold = 0.1,
  abs = TRUE,
  line = TRUE,
  var.order = "unadjusted",
  colors = c(col.control, col.exposed),
  shapes = c(17, 16),
  limits = c(m = c(0, 0.3)),
  title = "Covariate Balance Before and After IPW",
  labels = TRUE
)

## Covariate balance improved substantially following application of the
## stabilized inverse probability weights, indicating improved comparability
## between treatment groups in the weighted pseudo-population.