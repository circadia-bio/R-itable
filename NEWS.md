## R-itable 0.2.0  (2026-08)

* `build_grm()` gains an `mz_col` argument for twin cohorts: pairs (or larger
  groups, e.g. triplets) sharing a monozygotic-twin label get the correct
  MZ kinship (A = 1.0) instead of the ordinary full-sibling default
  (A = 0.5), via `kinship2::pedigree(relation = ...)`. Also fixes a latent
  bug where cohorts with character founder IDs (e.g. synthetic twin-pair
  parent labels) would error, because the previous `NA -> 0` coercion for
  missing parents assumed numeric IDs; missing parents are now coded `NA`
  uniformly, which is accepted by `kinship2::pedigree()` regardless of ID
  type and is verified bit-identical to the old behaviour on numeric-ID
  pedigrees.

* `build_household()` -- new. Builds the block-diagonal household/shared-
  environment indicator matrix from an `hhid` grouping vector, for use with
  `herit_ace()`.

* `herit_ace()` -- new. Fits sporadic / AE / CE / full-ACE nested
  variance-component models for a single trait and returns the two
  likelihood-ratio tests needed to judge A/C identifiability: CE vs
  sporadic, and ACE vs AE. Uses direct numerical ML over the full
  covariance matrix (Cholesky-based GLS, multi-start Nelder-Mead), since A
  and a household indicator do not in general share eigenvectors the way
  `herit_vc()`'s single-component model does.

* `herit_bivar()` / `herit_bivar_batch()` -- new. Bivariate AE model
  estimating the genetic correlation (rhoG), non-shared environmental
  correlation (rhoE), and derived phenotypic correlation (rhoP) between two
  traits, with LRTs for rhoG = 0, rhoE = 0, and the joint rhoG = rhoE = 0.
  `herit_bivar_batch()` iterates over a set of trait pairs and applies
  Benjamini-Hochberg FDR correction separately within each of the three
  p-value families, matching the `parse_and_fdr.py` SOLAR post-processing
  convention. Likelihood evaluations are reduced from an O(n^3) 2n x 2n
  Cholesky to n independent 2x2 blocks by diagonalising the GRM once
  (both covariance blocks are linear combinations of A and I, so they share
  A's eigenvectors) -- this is what makes batches over ~30 pairs tractable
  (~0.4s/pair vs ~64s/pair for the naive implementation). The optimiser
  seeds its variance-scale parameters from the empirical variance of each
  trait rather than assuming unit variance, which matters for any trait
  that isn't inverse-normal transformed (assuming unit variance previously
  caused catastrophic convergence failures -- correlations pinned to the
  +/-1 boundary -- for raw-scale traits with large or small variance).

  Validated against SOLAR Eclipse output on the Gomes et al. twin dataset
  (n=161, 31 trait pairs), matching SOLAR's exact choice of which
  phenotypes were inverse-normal transformed (9 of 17) vs left raw: mean
  |delta rhoG| = 0.019, mean |delta rhoE| = 0.013, with 6/31 pairs matching
  SOLAR to 4 decimal places on both rhoG and rhoE. Not bit-identical to
  SOLAR -- ML vs SOLAR's own optimiser and variance parametrisation will
  always leave some numerical daylight, especially for weakly-identified
  (non-significant) correlations where the likelihood surface is flat. The
  largest residual disagreement is `CRONOTIPO_MCTQ_m`, which has an
  anomalous raw variance (~2.3e18, values in the billions) suggestive of a
  units/encoding problem upstream in `pheno_novo.csv` -- flagged for
  Mario/Lucas to check independently of this package.

* `int_transform()` gains a `method` argument. `"vdw"` is now the
  **default**, and it exactly replicates SOLAR Eclipse's `inormal` command
  -- confirmed by reading SOLAR's actual source (`lib/solar.tcl`, `proc
  inormal`, cloned from `github.com/kochunov/solar-eclipse`), not just its
  documentation. The previous docstring's claim that R-itable's old
  Blom-style formula "is the standard transformation used in SOLAR
  Eclipse" was incorrect, and an intermediate fix within this same release
  (a plain Van der Waerden `qnorm(r/(n+1))`, correct only for untied data)
  was still subtly wrong under ties. SOLAR's actual tie-handling, now
  replicated exactly: for a block of tied values spanning sorted rank
  positions `i` to `j`, average `qnorm(i/(n+1)), ..., qnorm(j/(n+1))` --
  i.e. average the *z-scores* of the individual rank positions, not the
  z-score of the *averaged* rank. These differ whenever a nonlinear
  `qnorm` is applied to a tied block that isn't symmetric around the
  middle of the sample (e.g. a large block of tied zeros in a
  zero-inflated phenotype). `"blom"` (`qnorm((r-0.5)/n)`, conventional
  rank-average tie handling) is retained as an explicit option for
  continuity with R-itable <= 0.2.0. **This changes the default numeric
  output of `int_transform()`, `herit_vc()`, `herit_batch()`,
  `herit_ace()`, and `herit_bivar()`** for any call using the default
  `transform`/`method` settings, particularly for heavily-tied traits.
  Pass `method = "blom"` (for `int_transform()`) to reproduce pre-0.2.0
  numeric output exactly.

* `herit_bivar()` gains proper support for single-trait missing data.
  Previously, any individual missing either trait was dropped via a
  complete-case restriction. Individuals missing exactly one of the two
  traits are now retained and contribute their observed trait's marginal
  likelihood, matching SOLAR Eclipse's default `UnbalancedTraits` behaviour
  (only individuals missing *both* traits, or the ID, are excluded). This
  is both more statistically appropriate (uses more of the available data
  under the usual MAR assumption) and matches SOLAR much more closely on
  affected pairs. The returned list gains `n` (all contributing
  individuals) and `n_complete` (individuals with both traits observed);
  `n` is now the primary sample-size figure, matching SOLAR's convention.
  When neither trait has missing data (the common case) this is a no-op --
  the fast eigenbasis path is used exactly as before; the slower general
  Cholesky path only engages when there is missingness in one of the two
  traits.

  **Final validation against SOLAR Eclipse** on the Gomes et al. twin
  dataset (n=161, 31 trait pairs), combining every fix in this release
  (correct per-trait inverse-normal transform selection, the SOLAR-exact
  `inormal` tie handling above, data-driven optimiser starting values, and
  unbalanced-trait handling): **29 of 31 pairs now match SOLAR to 4
  decimal places exactly** on both rhoG and rhoE (mean |delta rhoG| =
  0.010, mean |delta rhoE| = 0.009, both driven entirely by the 2
  remaining pairs below; median |delta rhoG| = |delta rhoE| = 0.0000).
  All 31 pairs agree with SOLAR on nominal significance for both rhoG and
  rhoE. The bivariate model itself was also independently confirmed
  algebraically identical to SOLAR's own documented omega formula for
  multivariate polygenic models (SOLAR Manual Chapter 9). The 2 remaining
  pairs both involve the same phenotype:

  `CRONOTIPO_MCTQ_m` has an anomalous raw variance (~2.3e18, values in the
  billions) strongly suggestive of a units/encoding problem upstream in
  `pheno_novo.csv`, independent of this package -- flagged to Mario/Lucas
  separately, not something this package can or should paper over.

  (An earlier hypothesis in this release's development -- that the
  remaining gap on `DR_ATV_m`/`DR_ATM_m`/`DR_ATL_m` reflected an
  inherently hard, flat likelihood surface under heavy zero-inflation --
  turned out to be wrong; those pairs are now resolved exactly by the
  `inormal` tie-handling fix above, once the actual SOLAR source was read
  rather than assumed.)

## R-itable 0.1.0  (2026-05)

* `herit_vc()` and `herit_batch()` now return `var_covariates`: the proportion
  of phenotypic variance explained by fixed-effect covariates (R² on the
  INT-transformed phenotype). This corresponds to the "variance explained"
  column reported in Leocadio-Miguel et al. (2025, *J Sleep Res*). Returns
  `NA` for unadjusted models. `sigma2_a` and `sigma2_e` are now documented
  using SOLAR's sigma²_g / sigma²_e notation for clarity.

* `build_grm()` — build an additive genetic relationship matrix from a
  pedigree data frame via `kinship2::kinship()`. Supports custom column names,
  graceful handling of missing parents, and informative errors for common
  mistakes.

* `herit_vc()` — profile-likelihood variance-components heritability estimator
  for a single quantitative trait. Features: inverse-normal transformation,
  one-sided LRT with chi-squared(1) boundary correction, profile-likelihood
  95% CIs, and zero-variance covariate detection.

* `herit_batch()` — iterate `herit_vc()` over many traits x covariate model
  combinations; returns a tidy data frame. Includes a cli progress bar.

* `int_transform()` — exported rank-based inverse-normal transformation
  (Blom-style).

* `plot_forest()` — ggplot2 forest plot method for `herit_batch()` output,
  with optional model filtering and significance shading.
