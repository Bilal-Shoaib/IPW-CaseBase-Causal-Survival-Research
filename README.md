# Causal Effect Estimation for Restricted Mean Survival Time under Case-Base Sampling with Inverse Probability Weights

## Project Overview

This project contains two R Markdown analyses for causal restricted mean survival time (RMST) under case-base sampling with inverse probability weighting (IPW):

1. a simulation study evaluating the performance of the estimators across multiple data-generating scenarios; and
2. an application to the Right Heart Catheterization (RHC) dataset.

Both analyses are authored by Bilal Shaikh and Sumeet Kalia and are configured to produce PDF output.

---

# Repository Structure

The repository is organized into two primary directories:

```text
.
├── Application/
│   ├── RHC-Application.Rmd
│   ├── Figures/
│   └── Tables/
└── Simulation Study/
    ├── Simulation-Study.Rmd
    ├── Figures/
    └── Tables/
```

Each directory contains the R Markdown analysis together with the figures and tables produced by that analysis.

---

## 1. Application

The `Application` directory contains the R Markdown analysis titled:

> **Causal RMST under case-base sampling with IPW: Right Heart Catheterization**

The analysis uses the Right Heart Catheterization dataset and considers both an intention-to-treat analogue and a per-protocol analogue.

### `RHC-Application.Rmd`

The application proceeds through the following main components:

- configuration of the analysis;
- a common plotting theme;
- data preparation;
- construction of follow-up under the intention-to-treat and per-protocol principles;
- propensity score estimation and treatment weights;
- censoring weights;
- exploratory data analysis and covariate balance;
- model fitting;
- hazard-ratio results;
- RMST results;
- sensitivity analyses;
- case-base resampling variability; and
- session information.

The model-fitting section evaluates ten configurations. These include unadjusted and covariate-adjusted MS-Cox models, IPW and IPW-plus-IPCW MS-Cox models, corresponding case-base models, and time-varying hazard-ratio configurations.

The application calculates counterfactual survival curves and RMST contrasts for RHC and no RHC. RMST estimates are evaluated under both intention-to-treat and per-protocol analyses.

### `Figures/`

The application R Markdown file writes the following figures to this directory:

```text
Figures/
├── covariate-balance.pdf
├── propensity-overlap.pdf
├── rmst-contrast.pdf
├── counterfactual-survival.pdf
└── hr-time-varying.pdf
```

These figures cover covariate balance, propensity-score overlap, RMST contrasts, counterfactual survival, and time-varying hazard ratios.

### `Tables/`

The application R Markdown file writes the following tables to this directory:

```text
Tables/
├── weights.tex
├── hr.tex
├── rmst.tex
├── sensitivity.tex
└── resampling.tex
```

These tables contain summaries of the treatment and censoring weights, hazard-ratio results, RMST results, sensitivity analyses, and case-base resampling variability.

---

## 2. Simulation Study

The `Simulation Study` directory contains the R Markdown analysis titled:

> **Causal RMST under case-base sampling with IPW: Simulation Study**

The simulation evaluates the estimators across eight specified scenarios.

### `Simulation-Study.Rmd`

The simulation study contains the following main components:

- global configuration;
- data-generating mechanism;
- reference-value calculation;
- weighted case-base fitting;
- RMST contrast and delta-method standard-error calculation;
- a single simulation replicate;
- execution of the simulation study;
- performance-measure calculation;
- LaTeX table generation;
- reference values based on integration against a model-based alternative;
- figure generation;
- computation-time benchmarking; and
- session information.

The data-generating mechanism uses ten covariates. Treatment assignment is generated from a treatment linear predictor, and event and censoring times are generated according to the specified simulation mechanism.

The eight simulation scenarios are:

| Scenario | Description |
| --- | --- |
| S1 | Null treatment effect |
| S2 | Moderate treatment effect |
| S3 | Reference configuration |
| S4 | Reduced sample size |
| S5 | Small sample size |
| S6 | Heavy censoring |
| S7 | Strong confounding |
| S8 | Weak prognostic covariate effect |

The simulation is configured for 2,000 replicates per scenario. Reference values are obtained through Monte Carlo integration over the covariate distribution.

The estimators evaluated in an individual replicate include:

- `M1`: IPW MS-Cox model;
- `M2`: IPW case-base model with an exponential time basis;
- `M3`: IPW case-base model with a Gompertz time basis;
- `M4`: IPW case-base model with a spline time basis;
- `M5`: IPW MS-Cox model with a time-varying treatment effect; and
- `M6`: IPW case-base model with a time-varying treatment effect.

The performance summaries include bias, Monte Carlo standard error, empirical standard deviation, average estimated standard error, mean squared error, confidence-interval coverage, and rejection rate. Convergence-related quantities are also summarized.

The simulation is run in parallel and saves the resulting simulation object as:

```text
simulation-results.rds
```

### `Figures/`

The simulation R Markdown file writes the following figures to this directory:

```text
Figures/
├── rmst-curves.pdf
├── hazard-ratio-time.pdf
└── coverage.pdf
```

These figures summarize RMST curves, time-varying hazard ratios, and empirical coverage across the simulation scenarios.

### `Tables/`

The simulation R Markdown file writes the following tables to this directory:

```text
Tables/
├── reference.tex
├── hr.tex
├── rmst.tex
├── convergence.tex
└── timing.tex
```

These tables contain reference values, hazard-ratio results, RMST results, convergence summaries, and computation-time results.

---

# Running the Analyses

Both analyses are contained in their respective `.Rmd` files and specify:

```yaml
output: pdf_document
```

The analysis-specific directories are created by the R Markdown code when required.

### Application

The application can be rendered from:

```text
Application/RHC-Application.Rmd
```

The analysis generates its figures under:

```text
Application/Figures/
```

and its tables under:

```text
Application/Tables/
```

### Simulation Study

The simulation study can be rendered from:

```text
Simulation Study/Simulation-Study.Rmd
```

The analysis generates its figures under:

```text
Simulation Study/Figures/
```

and its tables under:

```text
Simulation Study/Tables/
```

The simulation study also performs a computation-time benchmark separately from the main simulation study.

---

# Computational Workflow

The two analyses follow complementary roles within the project:

```text
                    Causal RMST under
                 Case-Base Sampling with IPW
                            │
             ┌──────────────┴──────────────┐
             │                             │
             ▼                             ▼
     Simulation Study                 RHC Application
             │                             │
             ▼                             ▼
     Data-generating                 Data preparation
        mechanism                    and follow-up
             │                             │
             ▼                             ▼
     Reference values                 Weights and
             │                       covariate balance
             ▼                             │
     Repeated estimation                   ▼
             │                       Model fitting
             ▼                             │
     Performance measures                  ▼
             │                       HR and RMST
             ▼                           results
     Figures and tables                    │
             │                             ▼
             └──────────────┬──────────────┘
                            ▼
                     Figures and tables
```

The simulation study evaluates estimator behavior under controlled scenarios, while the RHC application applies the analysis to the Right Heart Catheterization data under two analytic principles.

---

# Output Organization

The project keeps analysis outputs separate from the R Markdown source files:

| Directory | Contents |
| --- | --- |
| `Application/Figures/` | Figures generated by the RHC application |
| `Application/Tables/` | Tables generated by the RHC application |
| `Simulation Study/Figures/` | Figures generated by the simulation study |
| `Simulation Study/Tables/` | Tables generated by the simulation study |

The `.Rmd` files contain the analysis code and the code used to generate these outputs.

---

# Session Information

Both R Markdown files include a final `Session information` section using:

```r
sessionInfo()
```

This records the R session information associated with each analysis.
