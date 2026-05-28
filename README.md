# 🧬 R-itable

**Profile-likelihood heritability estimation for family cohort studies — no SOLAR required.**

[![R](https://img.shields.io/badge/R-%3E%3D4.1-276DC3)](https://cran.r-project.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE.md)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![Status](https://img.shields.io/badge/status-early%20development-orange)](https://github.com/circadia-bio/R-itable)

---

> ⚠️ **R-itable is in early development and has not been formally tested.** The API may change without notice, estimation results have not yet been validated against a reference implementation such as SOLAR Eclipse, and the package has not undergone peer review. Use with caution and verify outputs independently before using in any research context.

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
│   ├── build_grm.R            # build additive GRM from pedigree
│   ├── int_transform.R        # rank-based inverse-normal transform
│   ├── herit_vc.R             # single-trait VC estimator
│   ├── herit_batch.R          # batch wrapper
│   └── plot_forest.R          # ggplot2 forest plot
├── tests/testthat/
│   ├── helper-fixtures.R      # shared synthetic pedigree/data
│   ├── test-build_grm.R
│   ├── test-int_transform.R
│   ├── test-herit_vc.R
│   └── test-herit_batch.R
├── vignettes/
│   └── getting-started.Rmd
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

# 3. Many traits x models -> tidy data frame
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

For a full walkthrough see `vignette("getting-started", package = "Ritable")`.

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
