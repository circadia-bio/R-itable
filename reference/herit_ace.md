# A vs C identifiability check (ACE model comparison)

Fits four nested variance-component models for a single quantitative
trait – sporadic (no familial effect), polygenic/AE (additive genetic
only), household/CE (common environment only), and full ACE (A and C
estimated jointly) – and returns their log-likelihoods together with the
two likelihood-ratio tests needed to judge whether A and C are
separately identifiable in the sample: CE vs sporadic (is there a
familial signal at all?) and ACE vs AE (does C add anything once A is in
the model?). This mirrors SOLAR Eclipse's `house` / `polygenic -screen`
model-comparison workflow (see `solar_full_analysis.tcl`, Part 3).

## Usage

``` r
herit_ace(
  trait,
  grm,
  household,
  data,
  id_col = "IID",
  transform = TRUE,
  min_n = 80L,
  verbose = TRUE
)
```

## Arguments

- trait:

  Character string: name of the trait column in `data`.

- grm:

  Numeric matrix: additive genetic relationship matrix, as returned by
  [`build_grm()`](https://r-itable.circadia-lab.uk/reference/build_grm.md).

- household:

  Numeric matrix: household indicator matrix, as returned by
  [`build_household()`](https://r-itable.circadia-lab.uk/reference/build_household.md).
  Must have the same row/column names as `grm`.

- data:

  Data frame containing `id_col` and `trait`.

- id_col:

  Name of the individual ID column in `data`. Default `"IID"`.

- transform:

  Logical. Apply
  [`int_transform()`](https://r-itable.circadia-lab.uk/reference/int_transform.md)
  to `trait` before fitting. Default `TRUE`, matching SOLAR Eclipse
  convention for skewed phenotypes.

- min_n:

  Minimum number of complete observations required. Default `80L`.

- verbose:

  Logical. Print progress via `cli`. Default `TRUE`.

## Value

A named list with elements:

- `trait`:

  Trait name.

- `n`:

  Sample size after dropping missing values.

- `loglik_sporadic`, `loglik_ae`, `loglik_ce`, `loglik_ace`:

  Log-likelihoods of the four nested models.

- `h2_ae`:

  Heritability under the AE model.

- `c2_ce`:

  Common-environment proportion under the CE model.

- `h2_ace`, `c2_ace`:

  Joint estimates under the full ACE model.

- `chisq_c_vs_sporadic`, `p_c_vs_sporadic`:

  One-sided boundary LRT of CE vs sporadic (evidence of a familial
  effect at all).

- `chisq_ace_vs_ae`, `p_ace_vs_ae`:

  One-sided boundary LRT of ACE vs AE (evidence C adds anything once A
  is estimated). A near-zero statistic with `c2_ace` at or near the 0
  boundary indicates A and C are not jointly identifiable in this sample
  – report the simpler AE model (see package vignette / Leocadio-Miguel
  et al. companion analyses).

Returns `NULL` if `n < min_n`.

## Details

Unlike
[`herit_vc()`](https://r-itable.circadia-lab.uk/reference/herit_vc.md),
which diagonalises a single variance component via
[`eigen()`](https://rdrr.io/r/base/eigen.html), this function estimates
one or two variance components jointly by direct numerical maximum
likelihood over the full `n x n` covariance matrix (Cholesky-based GLS),
since `A` and household indicator `C` do not in general share
eigenvectors.

**Model:** `Omega = sigma2_p * [h2*A + c2*C + (1 - h2 - c2)*I]`, jointly
estimated for the ACE model via multi-start Nelder-Mead on a
logit-reparametrised `(h2, c2)` (guaranteeing `h2, c2 >= 0` and
`h2 + c2 < 1`); AE and CE are 1-D special cases.

**LRTs:** one-sided chi-squared(1) boundary correction, as in
[`herit_vc()`](https://r-itable.circadia-lab.uk/reference/herit_vc.md),
since both null hypotheses (`c2 = 0`, or `h2 = 0`) sit on the boundary
of the parameter space.

## See also

[`build_grm()`](https://r-itable.circadia-lab.uk/reference/build_grm.md),
[`build_household()`](https://r-itable.circadia-lab.uk/reference/build_household.md),
[`herit_vc()`](https://r-itable.circadia-lab.uk/reference/herit_vc.md),
[`herit_bivar()`](https://r-itable.circadia-lab.uk/reference/herit_bivar.md)
