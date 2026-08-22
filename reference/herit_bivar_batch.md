# Batch bivariate genetic/environmental correlations with FDR correction

Iterates
[`herit_bivar()`](https://r-itable.circadia-lab.uk/reference/herit_bivar.md)
over a set of trait pairs and applies Benjamini-Hochberg FDR correction
separately within each of the three p-value families (rhoP, rhoG, rhoE),
matching `parse_and_fdr.py`'s post-processing of SOLAR output.

## Usage

``` r
herit_bivar_batch(
  pairs,
  grm,
  data,
  id_col = "IID",
  transform = TRUE,
  min_n = 80L,
  .progress = TRUE
)
```

## Arguments

- pairs:

  A two-column character matrix or data frame of trait name pairs
  (column 1 = trait1, column 2 = trait2), or a list of length-2
  character vectors.

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

- .progress:

  Logical. Show a cli progress bar. Default `TRUE`.

## Value

A data frame with one row per successfully fitted pair and columns
`trait1`, `trait2`, `n`, `n_complete`, `h2_1`, `h2_2`, `rhoG`, `rhoE`,
`rhoP`, `p_rhoG`, `p_rhoE`, `p_rhoP`, `q_rhoG`, `q_rhoE`, `q_rhoP`
(BH-adjusted). Failed / skipped pairs are silently omitted.

## See also

[`herit_bivar()`](https://r-itable.circadia-lab.uk/reference/herit_bivar.md)
