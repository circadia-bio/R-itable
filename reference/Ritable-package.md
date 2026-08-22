# Ritable: Pedigree-Based Heritability Estimation for Family Cohort Studies

Provides profile-likelihood and general maximum-likelihood
variance-components estimation of narrow-sense heritability (h2),
household/common-environment effects (c2), and bivariate genetic and
environmental correlations for quantitative traits in family cohort and
twin studies. Additive genetic relationship matrices are built from
pedigrees via 'kinship2', including monozygotic-twin relatedness
overrides. Phenotypes are inverse-normal transformed internally.
Likelihood-ratio tests use a one-sided chi-squared boundary correction
equivalent to SOLAR Eclipse. Ninety-five percent confidence intervals
are derived from the profile likelihood rather than Wald approximations.
Batch estimation over many traits or trait pairs returns tidy data
frames ready for downstream visualisation (forest plots, heatmaps) and
includes Benjamini-Hochberg FDR correction for bivariate correlation
batches.

## See also

Useful links:

- <https://r-itable.circadia-lab.uk>

- <https://github.com/circadia-bio/R-itable>

- Report bugs at <https://github.com/circadia-bio/R-itable/issues>

## Author

**Maintainer**: Lucas França <lucas.franca@northumbria.ac.uk>
([ORCID](https://orcid.org/0000-0003-0853-1319))

Authors:

- Lucas França <lucas.franca@northumbria.ac.uk>
  ([ORCID](https://orcid.org/0000-0003-0853-1319))

- Mario Leocadio-Miguel <mario.miguel@northumbria.ac.uk>
  ([ORCID](https://orcid.org/0000-0002-7248-3529))
