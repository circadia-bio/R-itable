# Package index

## Genetic relationship matrix

Build the additive GRM (and household matrix) from a pedigree data
frame.

- [`build_grm()`](https://r-itable.circadia-lab.uk/reference/build_grm.md)
  : Build an additive genetic relationship matrix from a pedigree
- [`build_household()`](https://r-itable.circadia-lab.uk/reference/build_household.md)
  : Build a household (shared-environment) indicator matrix

## Heritability estimation

Single-trait and batch profile-likelihood VC estimation.

- [`herit_vc()`](https://r-itable.circadia-lab.uk/reference/herit_vc.md)
  : Profile-likelihood variance-components heritability estimation
- [`herit_batch()`](https://r-itable.circadia-lab.uk/reference/herit_batch.md)
  : Batch heritability estimation over multiple traits

## Household effects & twin ACE models

A vs C identifiability for family and twin cohorts.

- [`herit_ace()`](https://r-itable.circadia-lab.uk/reference/herit_ace.md)
  : A vs C identifiability check (ACE model comparison)

## Bivariate genetic/environmental correlation

Genetic and environmental correlation between trait pairs.

- [`herit_bivar()`](https://r-itable.circadia-lab.uk/reference/herit_bivar.md)
  : Bivariate genetic and environmental correlation between two traits
- [`herit_bivar_batch()`](https://r-itable.circadia-lab.uk/reference/herit_bivar_batch.md)
  : Batch bivariate genetic/environmental correlations with FDR
  correction

## Utilities

Helper functions used internally and exported for pipelines.

- [`int_transform()`](https://r-itable.circadia-lab.uk/reference/int_transform.md)
  : Inverse-normal (rank-based) transformation
- [`Ritable_colours`](https://r-itable.circadia-lab.uk/reference/Ritable_colours.md)
  : R-itable colour palette

## Visualisation

Forest plots for heritability results.

- [`plot_forest()`](https://r-itable.circadia-lab.uk/reference/plot_forest.md)
  : Forest plot of heritability estimates

## Package

Package-level documentation.

- [`Ritable`](https://r-itable.circadia-lab.uk/reference/Ritable-package.md)
  [`Ritable-package`](https://r-itable.circadia-lab.uk/reference/Ritable-package.md)
  : Ritable: Pedigree-Based Heritability Estimation for Family Cohort
  Studies
