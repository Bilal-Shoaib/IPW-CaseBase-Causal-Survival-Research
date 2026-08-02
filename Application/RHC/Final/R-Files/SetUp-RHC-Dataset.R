# ==============================================================================
# Package Imports
# ==============================================================================

### Exploratory Data Analysis and Covariate Balance
library("Hmisc")
library("tableone")
library("xtable")
library("ggplot2")
library("survey")
library("cobalt")

### Survival Analysis and Causal Effect Estimation
library("survival")
library("casebase")
library("sandwich")
library("foreach")
library("doParallel")


# ==============================================================================
# Data Preparation
# ==============================================================================

## Set random seed for reproducibility
s = 123

## Load the Right Heart Catheterization (RHC) dataset
getHdata(rhc)


# ------------------------------------------------------------------------------
# Primary Analysis Variables
# ------------------------------------------------------------------------------

## Define follow-up time as the time from hospital admission to the last
## available contact or death date
time <- as.numeric(
  pmax(rhc$dthdte, rhc$lstctdte, na.rm = TRUE) - rhc$sadmdte
)

## Define the binary event indicator:
## 1 = death, 0 = censored
Delta <- ifelse(
  rhc$death == "Yes",
  1,
  0
)

## Define the binary treatment indicator:
## 1 = right heart catheterization (RHC) performed, 0 = not performed
A <- ifelse(
  rhc$swang1 == "RHC",
  1,
  0
)


# ------------------------------------------------------------------------------
# Baseline Covariates
# ------------------------------------------------------------------------------

### Demographic and Socioeconomic Characteristics

## Patient sex
Gender <- rhc$sex

## Patient age at baseline
Age <- as.numeric(rhc$age)

## Patient weight in kilograms
## Note: The original dataset encodes some missing values as 0. These values
## are retained because replacing them with NA substantially alters the
## estimated propensity scores and covariate balance following weighting.
Weight <- as.numeric(rhc$wtkilo1)

## Patient race
Race <- droplevels(rhc$race)

## Educational attainment
Education <- as.numeric(rhc$edu)

## Income category
Income <- rhc$income

## Insurance classification
InsuranceClass <- rhc$ninsclas


### Pre-existing Comorbidities

## Names of baseline comorbidity variables to be recoded
hx.vars <- c(
  "chfhx",
  "cardiohx",
  "dementhx",
  "psychhx",
  "chrpulhx",
  "renalhx",
  "liverhx",
  "gibledhx",
  "malighx",
  "immunhx",
  "transhx",
  "amihx"
)

## Convert labelled comorbidity variables to binary factors:
## "No" = absence of comorbidity, "Yes" = presence of comorbidity
rhc[hx.vars] <- lapply(rhc[hx.vars], function(x) {
  x <- as.numeric(x)
  
  factor(
    ifelse(x == 1, "Yes", "No"),
    levels = c("No", "Yes")
  )
})

## Extract recoded comorbidity variables for the analysis dataset
CHFhx <- rhc$chfhx
CARDIOhx <- rhc$cardiohx
DEMENThx <- rhc$dementhx
PSYCHhx <- rhc$psychhx
CHRPULhx <- rhc$chrpulhx
RENALhx <- rhc$renalhx
LIVERhx <- rhc$liverhx
GIBLEDhx <- rhc$gibledhx
MALIGhx <- rhc$malighx
IMMUNhx <- rhc$immunhx
TRANShx <- rhc$transhx
AMIhx <- rhc$amihx


### Disease Severity and Clinical Characteristics

## Indicator of cancer diagnosis
Cancer <- rhc$ca

## Primary disease category
Category <- droplevels(rhc$cat1)

## Baseline disease severity score
Severity <- as.numeric(rhc$das2d3pc)


# ------------------------------------------------------------------------------
# Construct Analysis Dataset
# ------------------------------------------------------------------------------

## Combine the outcome, treatment, and baseline covariates into a single
## analysis-ready dataset
data.set <- data.frame(
  time = time,
  Delta = Delta,
  A = A,
  
  Gender = Gender,
  Age = Age,
  Weight = Weight,
  Race = Race,
  Education = Education,
  Income = Income,
  InsuranceClass = InsuranceClass,
  
  CHFhx = CHFhx,
  CARDIOhx = CARDIOhx,
  DEMENThx = DEMENThx,
  PSYCHhx = PSYCHhx,
  CHRPULhx = CHRPULhx,
  RENALhx = RENALhx,
  LIVERhx = LIVERhx,
  GIBLEDhx = GIBLEDhx,
  MALIGhx = MALIGhx,
  IMMUNhx = IMMUNhx,
  TRANShx = TRANShx,
  AMIhx = AMIhx,
  
  Cancer = Cancer,
  Category = Category,
  Severity = Severity
)