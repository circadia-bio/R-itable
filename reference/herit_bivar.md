# Bivariate genetic and environmental correlation between two traits

Fits a bivariate AE variance-components model for two quantitative
traits and estimates the additive-genetic correlation (rhoG), the
non-shared environmental correlation (rhoE), and the derived phenotypic
correlation (rhoP = rhoG*sqrt(h2_1*h2_2) +
rhoE*sqrt((1-h2_1)*(1-h2_2))), together with likelihood-ratio tests of
rhoG = 0, rhoE = 0, and the joint test of rhoP = 0 (rhoG = rhoE = 0).
This mirrors SOLAR Eclipse's `polygenic -testrhog` / `-testrhoe`
bivariate output (see `solar_full_analysis.tcl`, Part 5).

## Usage

``` r
herit_bivar(
  trait1,
  trait2,
  grm,
  data,
  id_col = "IID",
  transform = TRUE,
  min_n = 80L,
  verbose = TRUE
)
```

## Arguments

- trait1, trait2:

  Character strings: names of the two trait columns in `data`.

- grm:

  Numeric matrix: additive genetic relationship matrix, as returned by
  [`build_grm()`](https://r-itable.circadia-lab.uk/reference/build_grm.md).

- data:

  Data frame containing `id_col`, `trait1`, and `trait2`.

- id_col:

  Name of the individual ID column in `data`. Default `"IID"`.

- transform:

  Logical. Apply
  [`int_transform()`](https://r-itable.circadia-lab.uk/reference/int_transform.md)
  to both traits before fitting. Default `TRUE`.

- min_n:

  Minimum number of individuals contributing to the analysis sample
  required (an individual counts if at least one of `trait1`, `trait2`
  is observed – see Details). Default `80L`.

- verbose:

  Logical. Print progress via `cli`. Default `TRUE`.

## Value

A named list with elements:

- `trait1`, `trait2`:

  Trait names.

- `n`:

  Number of individuals contributing to the analysis (union of those
  with `trait1` and/or `trait2` observed – see Details).

- `n_complete`:

  Number of individuals with both traits observed. Equal to `n` when
  neither trait has missing data.

- `h2_1`, `h2_2`:

  Heritability of each trait under the joint bivariate model.

- `rhoG`, `rhoE`, `rhoP`:

  Genetic, environmental, and derived phenotypic correlations.

- `p_rhoG`, `p_rhoE`:

  LRT p-values for rhoG = 0 and rhoE = 0 (standard two-sided
  chi-squared(1); these correlations are not boundary parameters).

- `p_rhoP`:

  LRT p-value for the joint null rhoG = rhoE = 0 (chi-squared(2)).

Returns `NULL` if `n < min_n`.

## Details

**Model:** traits are stacked into a `2n`-vector with block covariance
`Sigma = [[sigma2_g1*A + sigma2_e1*I, rhoG*sqrt(sigma2_g1*sigma2_g2)*A + rhoE*sqrt(sigma2_e1*sigma2_e2)*I], [..., sigma2_g2*A + sigma2_e2*I]]`,
where `sigma2_g_i = h2_i * sigma2_p_i` and
`sigma2_e_i = (1-h2_i) * sigma2_p_i`. All six free parameters
(`h2_1, h2_2, rhoG, rhoE, sigma2_p1, sigma2_p2`) are estimated jointly
by multi-start Nelder-Mead ML on a reparametrised scale (logit for `h2`,
`atanh` for the correlations, `log` for the variances) so that all
constraints hold automatically.

**Missing data:** if either trait has missing values, individuals
missing exactly one of the two traits are still included in the analysis
sample, contributing their observed trait's marginal likelihood,
matching SOLAR Eclipse's default `UnbalancedTraits` behaviour
(individuals are only excluded if *both* traits, or the ID, are
missing). This is both more statistically appropriate than dropping them
(more information used under the usual MAR assumption) and matches SOLAR
exactly once combined with the SOLAR-exact `inormal` tie handling in
[`int_transform()`](https://r-itable.circadia-lab.uk/reference/int_transform.md)
– see NEWS.md for validation. When neither trait has missing data (the
common case), this reduces to the same fast eigenbasis computation as
before – the unbalanced path is only used when needed, since it cannot
exploit the eigenbasis shortcut (the missingness pattern is
individual-specific, so the joint covariance no longer shares a common
eigenbasis with `A` across all n individuals).

**Numerical robustness:** both traits are standardised to unit variance
before fitting. h2 and the correlations are scale-invariant, so this
changes nothing about the true ML answer, but avoids a real failure mode
for raw (non-INT-transformed) traits with an extreme variance – e.g. a
data-entry problem producing values in the billions – where the
unstandardised covariance matrix loses all floating-point precision.
rhoG/rhoE are also bounded to `[-0.9, 0.9]`, matching SOLAR's actual
default parameter bounds (confirmed from source, not just documentation
– see NEWS.md) rather than the mathematical `[-1, 1]`.

As with
[`herit_ace()`](https://r-itable.circadia-lab.uk/reference/herit_ace.md),
the unbalanced path uses direct numerical ML over the full covariance
matrix rather than an eigendecomposition shortcut.

Validated on the Gomes et al. twin dataset (31 trait pairs), combining
correct per-trait inverse-normal transform selection, SOLAR-exact
`inormal` tie handling (see
[`int_transform()`](https://r-itable.circadia-lab.uk/reference/int_transform.md)),
unbalanced-trait handling, and the numerical robustness fixes above: 28
of 31 pairs match SOLAR Eclipse to 4 decimal places exactly on both rhoG
and rhoE, and all 31 pairs agree on nominal significance for both. The 2
remaining pairs both involve a single phenotype (`CRONOTIPO_MCTQ_m` in
the validation dataset) that produces wildly different results from
SOLAR's reported values even under a fully numerically-robust fit –
strong evidence that the underlying data values differ between this
package's copy and whatever SOLAR was actually run on, not a modelling
or numerical gap in this package. See NEWS.md.

## See also

[`build_grm()`](https://r-itable.circadia-lab.uk/reference/build_grm.md),
[`herit_vc()`](https://r-itable.circadia-lab.uk/reference/herit_vc.md),
[`herit_ace()`](https://r-itable.circadia-lab.uk/reference/herit_ace.md),
[`herit_bivar_batch()`](https://r-itable.circadia-lab.uk/reference/herit_bivar_batch.md)
