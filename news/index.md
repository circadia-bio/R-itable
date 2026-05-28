# Changelog

## itable 0.1.0 (2026-05)

### New features

- [`build_grm()`](https://r-itable.circadia-lab.uk/reference/build_grm.md)
  — build an additive genetic relationship matrix from a pedigree data
  frame via
  [`kinship2::kinship()`](https://rdrr.io/pkg/kinship2/man/kinship.html).
  Supports custom column names, graceful handling of missing parents,
  and informative errors for common mistakes.

- [`herit_vc()`](https://r-itable.circadia-lab.uk/reference/herit_vc.md)
  — profile-likelihood variance-components heritability estimator for a
  single quantitative trait. Features: inverse-normal transformation,
  one-sided LRT with chi²(1) boundary correction, profile-likelihood 95%
  CIs, and zero-variance covariate detection.

- [`herit_batch()`](https://r-itable.circadia-lab.uk/reference/herit_batch.md)
  — iterate
  [`herit_vc()`](https://r-itable.circadia-lab.uk/reference/herit_vc.md)
  over many traits × covariate model combinations; returns a tidy data
  frame. Includes a cli progress bar.

- [`int_transform()`](https://r-itable.circadia-lab.uk/reference/int_transform.md)
  — exported rank-based inverse-normal transformation (Blom-style).

- [`plot_forest()`](https://r-itable.circadia-lab.uk/reference/plot_forest.md)
  — ggplot2 forest plot method for
  [`herit_batch()`](https://r-itable.circadia-lab.uk/reference/herit_batch.md)
  output, with optional model filtering and significance shading.
