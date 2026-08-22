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
  (~0.4s/pair vs ~64s/pair for the naive implementation).

  Validated against SOLAR Eclipse output on the Gomes et al. twin dataset
  (n=161, 31 trait pairs): mean |delta rhoG| ~ 0.05, mean |delta rhoE| ~
  0.02, with matching nominal-significance calls on 29/31 (rhoG) and 31/31
  (rhoE) pairs. Not bit-identical to SOLAR -- differences arise from ML vs
  SOLAR's own optimiser and variance parametrisation -- but the substantive
  findings (including the depression-anxiety and sleep-quality signals)
  replicate.

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
