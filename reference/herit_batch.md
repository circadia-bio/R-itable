# Batch heritability estimation over multiple traits

A convenience wrapper around
[`herit_vc()`](https://r-itable.circadia-lab.uk/reference/herit_vc.md)
that iterates over a vector of trait names, optionally across multiple
covariate models, and returns a tidy data frame.

## Usage

``` r
herit_batch(
  traits,
  grm,
  data,
  covs_list = NULL,
  id_col = "IID",
  min_n = 80L,
  ci_level = 0.95,
  .progress = TRUE
)
```

## Arguments

- traits:

  Character vector of trait column names in `data`.

- grm:

  Numeric matrix: additive genetic relationship matrix as returned by
  [`build_grm()`](https://r-itable.circadia-lab.uk/reference/build_grm.md).

- data:

  Data frame containing ID, trait, and covariate columns.

- covs_list:

  A named list of covariate vectors, where each element defines one
  covariate model. If `NULL`, a single unadjusted model is run. Example:
  `list(unadj = NULL, cov1 = c("age", "sex"), cov2 = c("age", "sex", "age2"))`.

- id_col:

  Name of the individual ID column. Default `"IID"`.

- min_n:

  Minimum sample size to attempt estimation. Default `80`.

- ci_level:

  Profile-likelihood CI level. Default `0.95`.

- .progress:

  Logical. Show a cli progress bar. Default `TRUE`.

## Value

A data frame (tibble-compatible) with one row per successfully fitted
model and columns: `label`, `trait`, `covariates`, `n`, `h2`, `se`,
`ci_lo`, `ci_hi`, `pval`, `sigma2_a`, `sigma2_e`. Failed / skipped
models are silently omitted.

## See also

[`herit_vc()`](https://r-itable.circadia-lab.uk/reference/herit_vc.md),
[`build_grm()`](https://r-itable.circadia-lab.uk/reference/build_grm.md),
[`plot_forest()`](https://r-itable.circadia-lab.uk/reference/plot_forest.md)

## Examples

``` r
if (FALSE) { # \dontrun{
A    <- build_grm(my_pedigree, study_ids = my_data$IID)

res  <- herit_batch(
  traits     = c("bmi", "systolic_bp", "hdl"),
  grm        = A,
  data       = my_data,
  covs_list  = list(
    unadj = NULL,
    cov1  = c("age", "sex"),
    cov2  = c("age", "sex", "age2")
  )
)

# Significant adjusted models
subset(res, grepl("cov2", label) & pval < 0.05)
} # }
```
