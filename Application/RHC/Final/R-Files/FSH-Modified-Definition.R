# ==============================================================================
# Modified Case-Base Smooth Hazard Model
# ==============================================================================

## Modified implementation of fitSmoothHazard() from the casebase package.
## The function extends the standard case-base fitting procedure to allow
## observation-level weights, including inverse probability weights (IPWs).
FSH <- function(
    formula,
    data,
    time,
    w = rep(1, nrow(data)),
    ratio = 100, ...  
) {
  
  ## Preserve the original function call for reproducibility and method
  ## dispatch in downstream model summaries
  cl <- match.call()
  
  ## Explicitly expand "." in the model formula so that all covariates in the
  ## supplied dataset are represented in the model specification
  formula <- expand.dot.formula(
    formula,
    data = data
  )
  
  ## Extract the event variable from the left-hand side of the model formula
  eventVar <- all.vars(formula[[2]])
  
  ## If the follow-up time variable is not explicitly supplied, infer it from
  ## the data using the event variable and the standard case-base structure
  if (missing(time)) {
    varNames <- checkArgsTimeEvent(
      data = data,
      event = eventVar
    )
    
    timeVar <- varNames$time
    
  } else {
    ## Otherwise, use the time variable explicitly supplied by the user
    timeVar <- time
  }
  
  ## Identify the distinct event types present in the data to determine whether
  ## the analysis involves a single event or competing risks
  typeEvents <- sort(unique(data[[eventVar]]))
  
  
  # ----------------------------------------------------------------------------
  # Case-Base Sampling
  # ----------------------------------------------------------------------------
  
  ## Generate the case-base sample when the input data have not already been
  ## sampled. The original dataset is retained so that the supplied observation
  ## weights can be incorporated into the case-base sampling procedure.
  if (!inherits(data, "cbData")) {
    
    originalData <- as.data.frame(data)
    originalData$w <- w
    
    ## Sample cases and controls from the person-time data. The sampling ratio
    ## controls the number of person-moments sampled relative to observed cases.
    sampleData <- sampleCaseBase(
      originalData,
      timeVar,
      eventVar,
      comprisk = (length(typeEvents) > 2),
      ratio
    )
    
  } else {
    ## If the input is already a case-base dataset, retain it directly.
    ## The original unsampled data are unavailable at this stage.
    originalData <- data.frame()
    sampleData <- data
  }
  
  
  # ----------------------------------------------------------------------------
  # Weighted Case-Base Model
  # ----------------------------------------------------------------------------
  
  ## Include the case-base sampling offset in the model formula. This offset
  ## accounts for the sampling mechanism when estimating the hazard function.
  formula <- update(formula, ~ . + offset(offset))
  
  ## Fit a weighted binomial regression model to the case-base sample.
  ## Observation-level weights allow the model to incorporate IPWs in addition
  ## to the case-base sampling structure.
  fittingFunction <- function(formula) glm(
    formula,
    data = sampleData,
    family = binomial,
    weights = w
  )
  
  ## Estimate the smooth hazard model using the weighted case-base sample
  out <- fittingFunction(formula)
  
  
  # ----------------------------------------------------------------------------
  # Model Object and Metadata
  # ----------------------------------------------------------------------------
  
  ## Preserve the lower-level glm call while replacing the displayed call with
  ## the original FSH call for consistency with the custom model interface
  out$lower_call <- out$call 
  out$call <- cl
  
  ## Store the original data and key variable names for downstream methods
  ## such as prediction and plotting
  out$originalData <- originalData
  out$typeEvents <- typeEvents
  out$timeVar <- timeVar
  out$eventVar <- eventVar
  
  ## Record the number of case and background person-moments in the sampled
  ## case-base dataset
  num_pm <- table(sampleData[[eventVar]])
  out$num_cm <- num_pm[2]
  out$num_bm <- num_pm[1]
  out$formula <- formula
  
  ## Store the case-base sampling offset separately, then reset the model's
  ## offset to zero so that absolute risk estimation can be performed without
  ## retaining the sampling correction in the prediction stage
  out$offset <- out$data$offset
  out$data$offset <- 0
  
  ## Assign a custom class so that methods defined for single-event case-base
  ## models can be applied to the fitted object
  class(out) <- c("singleEventCB", class(out))
  
  return(out)
}


# ==============================================================================
# Formula Expansion Helper
# ==============================================================================

## Expand "." in a model formula into the complete set of covariates available
## in the supplied dataset. This ensures that the modified FSH function can
## explicitly evaluate formulas containing "." before fitting the model.
expand.dot.formula <- function(formula, data = NULL) {
  
  ## Check whether the formula contains the "." shorthand for all available
  ## predictors
  if (isTRUE("." %in% all.vars(formula))) {
    
    ## Preserve the original formula attributes while expanding the terms
    att <- attributes(formula)
    
    try.terms <- try(
      stats::terms(formula, data = data),
      silent = TRUE
    )
    
    ## Replace the shorthand formula with its explicitly expanded version when
    ## the terms object can be constructed successfully
    if (!is(try.terms, "try-error")) {
      formula <- formula(try.terms)
    }
    
    ## Restore the original formula attributes after expansion
    attributes(formula) <- att
  }
  
  ## Return the expanded formula
  return(formula)
}