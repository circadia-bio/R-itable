# 🧬 R-itable

**Profile-likelihood heritability estimation for family cohort studies — no SOLAR required.**

[![R](https://img.shields.io/badge/R-%3E%3D4.1-276DC3)](https://cran.r-project.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE.md)
[![R CMD CHECK](https://github.com/circadia-bio/R-itable/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/circadia-bio/R-itable/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![Status](https://img.shields.io/badge/status-early%20development-orange)](https://github.com/circadia-bio/R-itable)

---

> ⚠️ **R-itable is in early development.** The API may change without notice, and the package has not undergone peer review. `herit_ace()`/`herit_bivar()`/`herit_bivar_batch()` have been validated against real SOLAR Eclipse output (28 of 31 real-cohort trait pairs match to 4 decimal places -- see `NEWS.md` for methodology and the two known exceptions); `herit_vc()`/`herit_batch()` have not yet undergone the same formal comparison. Verify outputs independently before using in any research context.

---

## 📖 What is R-itable?

**R-itable** (`library(Ritable)`) estimates narrow-sense heritability (h²) for
quantitative traits in pedigree-based family cohort studies. It implements a
profile-likelihood variance-components approach equivalent to SOLAR Eclipse —
without any proprietary dependencies, compiled code, or external binaries.

Built for neuroimaging and biomedical cohorts where you need to run heritability
over dozens of traits across multiple covariate models and get results you can
trace back to first principles.

---

## ✨ Features

- 🔬 **Profile-likelihood VC estimator** — eigendecomposition of the GRM, 1-D
  optimisation, exact LRT with one-sided χ²(1) boundary correction (matching SOLAR).
- 👯 **Twin cohorts** — `build_grm(mz_col = ...)` gives monozygotic pairs (or
  larger groups) the correct kinship of 1.0, instead of the ordinary
  full-sibling 0.5.
- 🏠 **Household/common-environment effects** — `herit_ace()` fits nested
  sporadic/AE/CE/ACE models and tests whether A and C are jointly
  identifiable in your sample, with optional SOLAR-`-screen`-style
  automatic covariate screening (`covs = ...`).
- 🔗 **Bivariate genetic/environmental correlations** — `herit_bivar()` /
  `herit_bivar_batch()` estimate rhoG, rhoE, and derived rhoP between trait
  pairs, with LRTs and Benjamini-Hochberg FDR correction across a batch.
- 📐 **Profile-likelihood CIs** — not Wald ±1.96 SE; proper asymmetric intervals
  via `uniroot()`.
- 🔄 **INT transformation** — inverse-normal transform applied internally;
  also exported as `int_transform()` for use in other pipelines.
- 📦 **Batch mode** — `herit_batch()` iterates over traits × covariate models
  and returns a tidy data frame, ready for tables and figures.
- 🌲 **Forest plots** — `plot_forest()` for immediate visualisation of batch
  output (requires ggplot2).
- 🧩 **Minimal dependencies** — core functions require only base R, `kinship2`,
  `rlang`, and `cli`.

---

## 🗂️ Project Structure

```
R-itable/
├── R/
│   ├── Ritable-package.R      # package-level docs and colour palette
│   ├── build_grm.R            # build additive GRM from pedigree (incl. MZ twins)
│   ├── build_household.R      # build household/shared-environment matrix
│   ├── int_transform.R        # rank-based inverse-normal transform
│   ├── herit_vc.R             # single-trait VC estimator
│   ├── herit_batch.R          # batch wrapper
│   ├── herit_ace.R            # A vs C identifiability (sporadic/AE/CE/ACE)
│   ├── herit_bivar.R          # bivariate genetic/environmental correlation
│   ├── vc_utils.R             # internal general-ML helpers (ACE, bivariate)
│   └── plot_forest.R          # ggplot2 forest plot
├── tests/testthat/
│   ├── helper-fixtures.R      # shared synthetic pedigree/twin/data fixtures
│   ├── test-build_grm.R
│   ├── test-build_household.R
│   ├── test-int_transform.R
│   ├── test-herit_vc.R
│   ├── test-herit_batch.R
│   ├── test-herit_ace.R
│   └── test-herit_bivar.R
├── vignettes/
│   ├── getting-started.Rmd
│   └── twin-cohorts.Rmd       # MZ relatedness, household effects, bivariate correlation
├── data-raw/
│   └── prepare_data.R
├── DESCRIPTION
├── NAMESPACE
└── NEWS.md
```

---

## 🚀 Getting Started

### Prerequisites

```r
install.packages(c("kinship2", "rlang", "cli"))

# For plotting:
install.packages("ggplot2")
```

### Installation

```r
# From GitHub (recommended while pre-CRAN)
remotes::install_github("circadia-bio/R-itable")
```

### Basic usage

```r
library(Ritable)

# 1. Build GRM from pedigree
A <- build_grm(my_pedigree, study_ids = my_data$IID)

# 2. Single trait
herit_vc("bmi", grm = A, data = my_data, covs = c("age", "sex"))

# 3. Many traits x models -> tidy data frame with columns:
#    label, trait, covariates, n, h2, se, ci_lo, ci_hi, pval,
#    var_covariates, sigma2_a, sigma2_e
res <- herit_batch(
  traits    = c("bmi", "hdl", "systolic_bp"),
  grm       = A,
  data      = my_data,
  covs_list = list(
    unadj = NULL,
    cov1  = c("age", "sex"),
    cov2  = c("age", "sex", "age2")
  )
)

# 4. Forest plot
plot_forest(res, model_filter = "cov2")
```

### Twin cohorts: MZ relatedness, household effects, bivariate correlations

```r
# GRM with MZ-twin relatedness override (mztwin_col groups MZ pairs/triplets)
A <- build_grm(twin_pedigree, study_ids = twin_data$IID, mz_col = "mztwin")
C <- build_household(twin_pedigree$hhid, twin_pedigree$id)

# A vs C identifiability: sporadic / AE / CE / full ACE
herit_ace("anxiety_score", grm = A, household = C, data = twin_data, id_col = "IID")

# ...with automatic covariate screening (SOLAR's polygenic -screen algorithm:
# single-pass LRT test of each candidate, drop if p >= prob_level)
twin_data$age_sex <- twin_data$age * twin_data$sex
res <- herit_ace("anxiety_score", grm = A, household = C, data = twin_data,
                covs = c("age", "sex", "age_sex"), id_col = "IID")
res$covs_kept      # covariates that survived screening
res$covariate_lrt_p # each candidate's individual LRT p-value

# Bivariate genetic/environmental correlation between two traits
herit_bivar("anxiety_score", "depression_score", grm = A, data = twin_data, id_col = "IID")

# Many trait pairs at once, with BH-FDR correction per correlation type
pairs <- data.frame(
  trait1 = c("anxiety_score", "depression_score"),
  trait2 = c("sleep_quality", "sleep_quality")
)
herit_bivar_batch(pairs, grm = A, data = twin_data, id_col = "IID")
```

For a full walkthrough see `vignette("getting-started", package = "Ritable")`, and
`vignette("twin-cohorts", package = "Ritable")` for the workflow above in more detail.

---

## 📦 Dependencies

| Package | Role |
|---------|------|
| `kinship2` | Pedigree object and kinship matrix |
| `rlang` | Error/warning handling, tidy eval |
| `cli` | Progress bars and formatted messages |
| `ggplot2` *(Suggests)* | `plot_forest()` |

---

## 👥 Authors

| Role | Name | Affiliation |
|------|------|-------------|
| Author & maintainer | Lucas França | Northumbria University / Circadia Lab |
| Author | Mario Leocadio-Miguel | Northumbria University / Circadia Lab |

---

## 🤝 Related Tools

- 🧪 [**ptestR**](https://github.com/circadia-bio/ptestR) — permutation tests for R
- 🌙 [**SleepDiaries**](https://github.com/circadia-bio/SleepDiaries) — sleep diary PWA
- ⚡ [**ACTT_validation_study**](https://github.com/circadia-bio/ACTT_validation_study) — actigraphy validation
- 🔬 [**circadia-bio**](https://github.com/circadia-bio) — the Circadia Lab GitHub organisation

---

## 📄 Licence

Released under the [MIT License](./LICENSE.md).

Copyright © Circadia Lab, 2026
