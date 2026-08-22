# Build a household (shared-environment) indicator matrix

Constructs the block-diagonal indicator matrix **C** used as the shared
/common-environment ("household") random-effect design in
[`herit_ace()`](https://r-itable.circadia-lab.uk/reference/herit_ace.md):
`C[i, j] = 1` if individuals `i` and `j` belong to the same household
group, `0` otherwise (including the diagonal convention `C[i, i] = 1`).

## Usage

``` r
build_household(hhid, ids)
```

## Arguments

- hhid:

  Vector of household group labels, one per individual, in the same
  order as `ids`. Individuals sharing an `hhid` value are treated as
  co-resident (e.g. twin pairs reared together).

- ids:

  Vector of individual IDs, the same length as `hhid`. Used to name the
  rows/columns of the returned matrix so it aligns with a GRM from
  [`build_grm()`](https://r-itable.circadia-lab.uk/reference/build_grm.md).

## Value

A symmetric numeric matrix of dimension `length(ids) x length(ids)`,
with row/column names `as.character(ids)`.

## Examples

``` r
hh  <- c("A", "A", "B", "B", "C")
ids <- 1:5
build_household(hh, ids)
#>   1 2 3 4 5
#> 1 1 1 0 0 0
#> 2 1 1 0 0 0
#> 3 0 0 1 1 0
#> 4 0 0 1 1 0
#> 5 0 0 0 0 1
```
