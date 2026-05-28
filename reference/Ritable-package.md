# Ritable: Pedigree-Based Heritability Estimation for Family Cohort Studies

Provides profile-likelihood variance-components estimation of
narrow-sense heritability (h2) for quantitative traits in family cohort
studies. Additive genetic relationship matrices are built from pedigrees
via 'kinship2'. Phenotypes are inverse-normal transformed internally.
Likelihood-ratio tests use a one-sided chi-squared boundary correction
equivalent to SOLAR Eclipse. Ninety-five percent confidence intervals
are derived from the profile likelihood rather than Wald approximations.
Batch estimation over many traits returns tidy data frames ready for
downstream visualisation (forest plots, heatmaps).

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
