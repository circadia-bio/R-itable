## R-itable 0.2.0  (2026-08)

### ✨ New features

* `build_grm()` gains an `mz_col` argument for twin cohorts: pairs (or larger
  groups, e.g. triplets) sharing a monozygotic-twin label get the correct
  MZ kinship (A = 1.0) instead of the ordinary full-sibling default
  (A = 0.5), via `kinship2::pedigree(relation = ...)`.

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
  Cholesky to n independent 2x2 blocks by diagonalising the GRM once (both
  covariance blocks are linear combinations of A and I, so they share A's
  eigenvectors) -- this is what makes batches over ~30 pairs tractable
  (~0.4s/pair vs ~64s/pair for the naive implementation).

  **Validated against SOLAR Eclipse** on the Gomes et al. twin dataset
  (n=161, 31 trait pairs), combining every fix in this release (correct
  per-trait inverse-normal transform selection, SOLAR-exact `inormal` tie
  handling, SOLAR's actual `[-0.9,0.9]` rho bounds, and unbalanced-trait
  handling): **28 of 31 pairs match SOLAR to 4 decimal places exactly** on
  both rhoG and rhoE, and all 31 pairs agree with SOLAR on nominal
  significance for both. The bivariate model itself was also independently
  confirmed algebraically identical to SOLAR's own documented omega
  formula for multivariate polygenic models (SOLAR Manual Chapter 9). The
  2 remaining pairs both involve the same phenotype, `CRONOTIPO_MCTQ_m`
  (see Bug fixes below for why, and why it isn't one).

### 🐛 Bug fixes

* **`build_grm()`**: cohorts with character founder IDs (e.g. synthetic
  twin-pair parent labels like `"FA12"`) would error, because the previous
  `NA -> 0` coercion for missing parents assumed numeric IDs. Missing
  parents are now coded `NA` uniformly, which `kinship2::pedigree()`
  accepts regardless of ID type -- verified bit-identical to the old
  behaviour on numeric-ID pedigrees, so this is a pure bug fix, not a
  behaviour change for existing callers.

* **`int_transform()`**: the default formula didn't actually match SOLAR
  Eclipse's `inormal`, despite the docstring's claim that it did. Fixed in
  two passes once the claim was actually checked against SOLAR's source
  (`lib/solar.tcl`, `proc inormal`, cloned from
  `github.com/kochunov/solar-eclipse`) rather than assumed:
  1. Switched the untied-case formula from Blom (`qnorm((r-0.5)/n)`) to
     Van der Waerden (`qnorm(r/(n+1))`), which is what SOLAR actually uses.
  2. That alone was still subtly wrong under ties. SOLAR's real tie
     handling: for a block of tied values spanning sorted rank positions
     `i` to `j`, average `qnorm(i/(n+1)), ..., qnorm(j/(n+1))` -- i.e.
     average the *z-scores* of the individual rank positions, not the
     z-score of the *averaged* rank. These differ whenever a nonlinear
     `qnorm` is applied to a tied block that isn't symmetric around the
     middle of the sample (e.g. a large block of tied zeros in a
     zero-inflated phenotype). Now replicated exactly.

  `"blom"` is kept as an explicit `method` option for continuity with
  R-itable <= 0.2.0. **This changes the default numeric output of
  `int_transform()`, `herit_vc()`, `herit_batch()`, `herit_ace()`, and
  `herit_bivar()`** for any call using default settings, particularly for
  heavily-tied traits. Pass `method = "blom"` to reproduce pre-0.2.0
  output exactly. Confirmed to fully resolve the `DR_ATV_m`/`DR_ATM_m`/
  `DR_ATL_m` discrepancy against SOLAR described in an earlier draft of
  this changelog as an "inherently hard, flat likelihood surface" -- that
  diagnosis was wrong; it was this tie-handling bug.

* **`herit_bivar()`**: individuals missing exactly one of the two traits
  were dropped via a complete-case restriction. They're now retained and
  contribute their observed trait's marginal likelihood, matching SOLAR
  Eclipse's default `UnbalancedTraits` behaviour (only individuals missing
  *both* traits, or the ID, are excluded). The returned list gains `n`
  (all contributing individuals) and `n_complete` (individuals with both
  traits observed); `n` is now the primary sample-size figure, matching
  SOLAR's convention. When neither trait has missing data (the common
  case) this is a no-op -- the fast eigenbasis path is used exactly as
  before; the slower general Cholesky path only engages when needed.

* **`herit_bivar()`**: rhoG/rhoE were unbounded (mathematical `[-1, 1]`)
  instead of SOLAR's actual default `[-0.9, 0.9]` -- confirmed by reading
  SOLAR's source (`lib/solar.tcl`, `proc polymod`:
  `set rholow -0.9; set rhoup 0.9`), not assumed. Also, both traits are
  now standardised to unit variance before fitting: h2 and the
  correlations are scale-invariant so this changes nothing about the true
  ML answer, but it fixes a real numerical failure mode for raw
  (non-INT-transformed) traits with an extreme variance, where the
  unstandardised covariance matrix loses all floating-point precision
  during fitting (observed as correlations pinning to a boundary, or the
  optimiser failing to move from its starting point at all).

  Root-caused against `CRONOTIPO_MCTQ_m` (raw variance ~2.3e18 in
  `pheno_novo.csv`), which was the last pair not resolved by the fixes
  above. Checking this against SOLAR's own source explains why SOLAR
  never hits the failure this package did: its `sd`/`mean` parameters are
  "initialized... in c++" directly from the raw data rather than
  numerically searched, so SOLAR is architecturally immune to
  extreme-scale numerical failure regardless of what the values represent.
  Reproducing that robustness here (the standardisation fix above) does
  make the fit numerically well-behaved -- but the resulting
  maximum-likelihood estimate for this phenotype (rhoG in the 0.7-0.9
  range) is now *more* different from SOLAR's reported values (0.04-0.13)
  than the original broken fit was, not less. Since the same code
  reproduces SOLAR exactly on the other 29 pairs, the most likely
  explanation is that this package's copy of `pheno_novo.csv` has
  different `CRONOTIPO_MCTQ_m` values than whatever SOLAR was actually run
  on -- a data provenance question, not something fixable in this package.
  Flagged to Mario/Lucas to check against the original SOLAR input file.

### 🧪 Tests

* New `test-build_household.R`, `test-herit_ace.R`, `test-herit_bivar.R`,
  and a `make_twin_data()` fixture (MZ + DZ pairs with simulated shared
  genetic and household components) in `helper-fixtures.R`.
* `test-build_grm.R`: new MZ-twin cases, including a dedicated same-mother
  triplet pedigree (kinship2 requires MZ-marked individuals to share both
  sex and mother, which a naive fixture reuse doesn't guarantee).
* `test-herit_bivar.R`: unbalanced-trait (single-missing-trait) cases, and
  `n` vs `n_complete` behaviour.
* `test-int_transform.R`: `"vdw"` matches Van der Waerden exactly when
  untied; a dedicated tied-block case confirming z-scores are averaged
  (not ranks); `"blom"` matches the conventional rank-average formula.

### 🔧 Other changes

* `.Rbuildignore`: `.zenodo.json` and `CITATION.cff` excluded from
  `R CMD check` (both are legitimate top-level GitHub/Zenodo files, not
  standard R package structure).

### 📚 Documentation

* `_pkgdown.yml`: added `build_household`, `herit_ace`, `herit_bivar`, and
  `herit_bivar_batch` to the reference index (the pkgdown build was
  failing without this).
* README: covers the new twin/ACE/bivariate workflow, updated feature list
  and project structure tree.

## R-itable 0.1.0  (2026-05)

### 🌱 Initial release

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

### ✨ New features

* `herit_vc()` and `herit_batch()` now return `var_covariates`: the proportion
  of phenotypic variance explained by fixed-effect covariates (R² on the
  INT-transformed phenotype). This corresponds to the "variance explained"
  column reported in Leocadio-Miguel et al. (2025, *J Sleep Res*). Returns
  `NA` for unadjusted models.

### 📚 Documentation

* `sigma2_a` and `sigma2_e` are now documented using SOLAR's sigma²_g /
  sigma²_e notation for clarity.
