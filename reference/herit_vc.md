# Profile-likelihood variance-components heritability estimation

Estimates narrow-sense heritability (h²) for a single quantitative trait
using a profile-likelihood variance-components approach equivalent to
SOLAR Eclipse. The phenotype is inverse-normal transformed internally.

## Usage

``` r
herit_vc(
  trait,
  grm,
  data,
  covs = NULL,
  id_col = "IID",
  label = NULL,
  min_n = 80L,
  ci_level = 0.95,
  verbose = TRUE
)
```

## Arguments

- trait:

  Character string: name of the trait column in `data`.

- grm:

  Numeric matrix: the additive genetic relationship matrix for all
  individuals in `data`, as returned by
  [`build_grm()`](https://r-itable.circadia-lab.uk/reference/build_grm.md).
  Row and column names must match the ID column of `data`.

- data:

  Data frame containing `id_col`, `trait`, and any covariate columns.

- covs:

  Character vector of covariate column names, or `NULL` for an
  intercept-only model. Covariates are mean-centred and scaled to unit
  variance before fitting.

- id_col:

  Name of the individual ID column in `data`. Default `"IID"`.

- label:

  Optional string label for this model (used in batch output). If
  `NULL`, defaults to `"<trait>_adj"` or `"<trait>_unadj"`.

- min_n:

  Minimum number of complete observations required to attempt
  estimation. Models with fewer observations return `NULL` silently.
  Default `80`.

- ci_level:

  Confidence level for the profile-likelihood interval. Default `0.95`.

- verbose:

  Logical. Print progress to the console. Default `TRUE`.

## Value

A named list with elements:

- `label`:

  Model label.

- `trait`:

  Trait name.

- `covariates`:

  Covariate names joined by `"+"`, or `""`.

- `n`:

  Sample size after dropping missing values.

- `h2`:

  MLE of narrow-sense heritability.

- `se`:

  Standard error from profile-likelihood curvature (Wald).

- `ci_lo`, `ci_hi`:

  Profile-likelihood confidence interval bounds.

- `pval`:

  One-sided LRT p-value with chi-squared(1) boundary correction.

- `var_covariates`:

  Proportion of phenotypic variance explained by fixed-effect covariates
  (R² on INT-transformed phenotype). `NA` for unadjusted models.
  Corresponds to the "variance explained" column in Leocadio-Miguel et
  al. (2025).

- `sigma2_a`:

  Additive genetic variance (sigma²_g in SOLAR notation).

- `sigma2_e`:

  Residual environmental variance (sigma²_e).

Returns `NULL` if `n < min_n` or if the GRM subset is degenerate.

## Details

**Model:** Omega = sigma2_p \[h2 A + (1 - h2) I_n\]

**Optimisation:** eigendecomposition of A followed by 1-D
profile-likelihood optimisation over h2 in (0, 1) using a coarse grid to
seed [`stats::optimize()`](https://rdrr.io/r/stats/optimize.html).

**LRT:** one-sided chi-squared(1) boundary correction – the null is on
the boundary of the parameter space (h2 = 0), so the p-value is halved
relative to a standard chi-squared test. This matches the SOLAR Eclipse
convention.

**CIs:** derived by [`uniroot()`](https://rdrr.io/r/stats/uniroot.html)
on the profile log-likelihood, falling back to a Wald interval if the
root-finding fails.

## See also

[`build_grm()`](https://r-itable.circadia-lab.uk/reference/build_grm.md),
[`herit_batch()`](https://r-itable.circadia-lab.uk/reference/herit_batch.md),
[`int_transform()`](https://r-itable.circadia-lab.uk/reference/int_transform.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Build GRM then estimate heritability
A   <- build_grm(my_pedigree, study_ids = my_data$IID)
res <- herit_vc("bmi", grm = A, data = my_data,
                covs = c("age", "sex", "age2"))
str(res)
} # }
```
