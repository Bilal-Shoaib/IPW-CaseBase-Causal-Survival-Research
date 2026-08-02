# Causal Effect Estimation for Restricted Mean Survival Time under Case-Base Sampling with Inverse Probability Weights

## Project Overview

This repository contains the complete computational materials for the research project “Causal effect estimation for restricted mean survival time under Case-Base sampling with inverse probability weights.”

The project develops and evaluates a computationally efficient framework for causal survival analysis by extending the Case-Base modeling framework to accommodate inverse probability weights (IPW). The proposed approach enables estimation of marginal treatment effects while substantially reducing computational burden relative to the weighted marginal structural Cox (MS-Cox) model, particularly in analyses involving extremely large datasets.

The study considers two principal methodological contributions. First, the Case-Base framework is extended to incorporate inverse probability weighting, providing a computationally efficient alternative for estimating causal treatment effects in survival settings. Second, a closed-form expression is derived for estimating the standard errors of restricted mean survival time (RMST) curves using the delta method, avoiding the computational burden associated with repeated bootstrap model fitting.

The proposed methodology is evaluated through an extensive simulation study and subsequently demonstrated using a real-world clinical application based on the Right Heart Catheterization (RHC) dataset. Across the simulation scenarios considered, the proposed Case-Base approach demonstrated comparable statistical performance to the weighted MS-Cox model while offering substantial computational advantages.

---

## Key Contributions

The findings of this study highlight two methodological contributions to the causal survival analysis literature:

1. Extension of Case-Base modeling to causal survival analysis.
   The Case-Base modeling framework is extended to accommodate inverse probability weights, providing a computationally efficient approach for estimating marginal treatment effects while preserving the causal interpretation of the estimands.

2. Closed-form standard errors for RMST curves.
   A closed-form expression is derived for estimating the standard errors of restricted mean survival time curves using the delta method. This provides an analytical alternative to computationally intensive bootstrap procedures and facilitates efficient statistical inference for RMST.

The proposed methodology was evaluated through an extensive simulation study and subsequently demonstrated in a real-world clinical application. The simulation results showed comparable bias, precision, confidence interval coverage, and overall accuracy between the proposed Case-Base approach and the weighted MS-Cox model across the scenarios considered.

The real-world application further demonstrates the flexibility of the proposed framework by accommodating a time-dependent treatment effect in the RHC dataset. Collectively, the results demonstrate that substantial computational gains can be achieved without compromising statistical performance, making the proposed approach particularly attractive for large-scale survival analyses.

---

# Repository Structure

The repository is organized into three primary directories:

```text
.
├── Application/
├── Simulation Study/
└── Data Sets/
```

Each directory serves a distinct purpose within the research project.

---

## 1. Application

The `Application` directory contains the complete computational materials for the real-world application of the proposed methodology to the Right Heart Catheterization (RHC) dataset.

```text
Application/
├── Other Attempts/
└── RHC/
    ├── Final/
    │   ├── Rmd-Files/
    │   ├── R-Files/
    │   ├── Misc/
    │   └── Plots and Diagrams/
    └── Other versions/
```

### `Other Attempts/`

Contains exploratory and alternative implementations considered during the development of the application. These materials document earlier analytical approaches and are retained for completeness.

### `RHC/`

Contains the analyses conducted using the RHC dataset.

#### `Final/`

Contains the finalized version of the RHC analysis and all associated computational materials.

##### `Rmd-Files/`

Contains the complete `.Rmd` file in which all code used for the final RHC analysis is combined into a single document.

This provides the full analysis as a single, sequential workflow.

##### `R-Files/`

Contains the same code used in the `.Rmd` file, separated into multiple `.R` scripts according to their respective analytical tasks.

This version is provided to improve organization, readability, and ease of navigation through the analysis.

##### `Misc/`

Contains miscellaneous supporting materials associated with the final RHC analysis.

##### `Plots and Diagrams/`

Contains plots, figures, and diagrams generated throughout the RHC analysis.

#### `Other versions/`

Contains earlier versions and alternative implementations of the RHC analysis that preceded the finalized version.

---

## 2. Simulation Study

The `Simulation Study` directory contains the complete computational materials used to evaluate the statistical and computational performance of the proposed methodology.

The directory follows the same organizational structure used for the finalized application:

```text
Simulation Study/
├── Rmd-Files/
├── R-Files/
├── Misc/
└── Plots and Diagrams/
```

### `Rmd-Files/`

Contains the complete `.Rmd` file containing the simulation study in its entirety, including the simulation setup, data generation, model fitting, performance evaluation, and results.

### `R-Files/`

Contains the same simulation code divided into separate `.R` files according to their respective analytical components.

The segmented scripts provide a more organized representation of the simulation workflow and facilitate inspection of individual components.

### `Misc/`

Contains miscellaneous supporting materials associated with the simulation study.

### `Plots and Diagrams/`

Contains plots, figures, and diagrams generated from the simulation study.

---

## 3. Data Sets

The `Data Sets` directory contains additional datasets that were considered during the development of the research project.

These datasets are not required to reproduce either the Application or the Simulation Study. They are retained as supplementary materials documenting datasets that were initially considered during the development of the project.

---

# Running the Code

## Requirements

The computational analyses are implemented in R and are intended to be run using RStudio.

To reproduce the analyses, users should have:

* R installed;
* RStudio installed; and
* all required R packages installed.

The repository should be downloaded or cloned in its entirety so that the directory structure and relative file paths remain intact.

## R Markdown Files

The primary method for reproducing the analyses is through the `.Rmd` files contained in the `Rmd-Files` directories.

### Application

To reproduce the real-world application:

1. Open the repository in RStudio.
2. Navigate to:

```text
Application/RHC/Final/Rmd-Files/
```

3. Open the `.Rmd` file.
4. Execute the code chunks sequentially from top to bottom.

### Simulation Study

To reproduce the simulation study:

1. Open the repository in RStudio.
2. Navigate to:

```text
Simulation Study/Rmd-Files/
```

3. Open the `.Rmd` file.
4. Execute the code chunks sequentially from top to bottom.

> Important: The code should be executed in its existing order. Later code chunks may depend on objects, functions, or results generated by earlier chunks.

---

# Code Organization

Both the Application and Simulation Study are provided in two complementary formats.

| Format       | Description                                                                 |
| ------------ | --------------------------------------------------------------------------- |
| `.Rmd` files | Complete analysis contained within a single, sequential R Markdown document |
| `.R` files   | The same code separated into logically organized scripts                    |

The `.Rmd` files provide the most direct route for reproducing the complete analyses. The segmented `.R` files provide a more convenient format for reviewing, navigating, and modifying individual components of the computational workflow.

---

# Research Outputs

The repository contains the computational materials underlying the methodological development, simulation evaluation, and real-world application presented in this research.

The overall workflow can be summarized as:

```text
Methodological Development
          │
          ▼
    Simulation Study
          │
          ▼
  Statistical Evaluation
          │
          ▼
    RHC Application
          │
          ▼
      Final Results
```

The simulation study evaluates the statistical performance and computational efficiency of the proposed methodology, while the RHC application demonstrates its practical implementation in a real-world clinical setting.

---

# Citation

*Citation information will be added following publication and archival of the research materials.*

[Archive link to be added]