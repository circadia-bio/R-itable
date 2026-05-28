## itable 0.1.0  (2026-05)

### New features

* `build_grm()` — build an additive genetic relationship matrix from a
  pedigree data frame via `kinship2::kinship()`. Supports custom column names,
  graceful handling of missing parents, and informative errors for common
  mistakes.

* `herit_vc()` — profile-likelihood variance-components heritability estimator
  for a single quantitative trait. Features: inverse-normal transformation,
  one-sided LRT with chi²(1) boundary correction, profile-likelihood 95% CIs,
  and zero-variance covariate detection.

* `herit_batch()` — iterate `herit_vc()` over many traits × covariate model
  combinations; returns a tidy data frame. Includes a cli progress bar.

* `int_transform()` — exported rank-based inverse-normal transformation
  (Blom-style).

* `plot_forest()` — ggplot2 forest plot method for `herit_batch()` output,
  with optional model filtering and significance shading.
