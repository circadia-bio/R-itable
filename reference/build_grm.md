# Build an additive genetic relationship matrix from a pedigree

Constructs the additive genetic relationship matrix **A** (= 2 x kinship
matrix) for a set of study subjects, using a pedigree that may include
additional founder parents not in the study sample.

## Usage

``` r
build_grm(
  ped_df,
  study_ids = NULL,
  id_col = "id",
  pat_col = "pat",
  mom_col = "mom",
  sex_col = "sex",
  mz_col = NULL
)
```

## Arguments

- ped_df:

  A data frame with at least the columns `id`, `pat`, `mom`, and `sex`.
  Missing parents should be `NA` or `0`. `sex` should be numeric: `1` =
  male, `2` = female (any other value is recoded to `1` with a warning).

- study_ids:

  Character or integer vector of IDs for the study subjects whose
  sub-matrix should be extracted. Defaults to all IDs in `ped_df`.

- id_col, pat_col, mom_col, sex_col:

  Column names for the four required pedigree fields. Defaults are
  `"id"`, `"pat"`, `"mom"`, `"sex"`.

- mz_col:

  Optional column name in `ped_df` identifying monozygotic twin/multiple
  groups (e.g. a shared family label present only for MZ sets,
  blank/`NA` for everyone else, as in a `MZTWIN` column). All pairwise
  combinations of individuals sharing the same non-blank, non-`NA` value
  of `mz_col` are passed to
  [`kinship2::pedigree()`](https://rdrr.io/pkg/kinship2/man/pedigree.html)
  as monozygotic (`code = 1`), overriding the default full-sibling
  kinship of 0.25 with the correct MZ kinship of 0.5 (A = 1.0). Groups
  of size \\\geq\\ 3 (e.g. MZ triplets) are supported: every pair within
  the group is marked. `NULL` (default) disables this and reproduces
  prior behaviour exactly.

## Value

A symmetric numeric matrix of dimension
`length(study_ids) x length(study_ids)`, with diagonal entries equal to
1 for non-inbred individuals and row/column names matching `study_ids`.

## Details

The function calls
[`kinship2::kinship()`](https://rdrr.io/pkg/kinship2/man/kinship.html)
on the full pedigree (including founders), then multiplies by 2 to
obtain the additive relationship matrix, and finally subsets to
`study_ids`. Including founders in the pedigree ensures that kinship
coefficients between study subjects connected only through founders are
estimated correctly.

For twin cohorts, `pat`/`mom` are typically synthetic per-pair founder
labels rather than genotyped individuals; such founders should still be
present as their own rows in `ped_df` (with `pat`/`mom` set to `NA`) so
that
[`kinship2::pedigree()`](https://rdrr.io/pkg/kinship2/man/pedigree.html)
can resolve the family structure. Without an `mz_col`, twins are treated
as ordinary full siblings (A = 0.5); supplying `mz_col` corrects MZ
pairs to A = 1.0 while leaving DZ pairs (which already have the correct
full-sibling kinship) unchanged.

## Examples

``` r
# Minimal two-generation pedigree: two couples, four offspring
ped <- data.frame(
  id  = 1:8,
  pat = c(0, 0, 0, 0, 1, 1, 3, 3),
  mom = c(0, 0, 0, 0, 2, 2, 4, 4),
  sex = c(1, 2, 1, 2, 1, 2, 1, 2)
)
A <- build_grm(ped, study_ids = 5:8)
round(A, 3)
#>     5   6   7   8
#> 5 1.0 0.5 0.0 0.0
#> 6 0.5 1.0 0.0 0.0
#> 7 0.0 0.0 1.0 0.5
#> 8 0.0 0.0 0.5 1.0

# With an MZ-twin column: subjects 7 and 8 are MZ (forced to the same sex,
# since kinship2 requires MZ pairs to share sex -- ordinarily they would
# already match)
ped$sex[ped$id %in% c(7, 8)] <- 1
ped$mztwin <- c(NA, NA, NA, NA, NA, NA, "MZ1", "MZ1")
A_mz <- build_grm(ped, study_ids = 5:8, mz_col = "mztwin")
round(A_mz, 3)  # A[7,8] = 1.0 instead of 0.5
#>     5   6 7 8
#> 5 1.0 0.5 0 0
#> 6 0.5 1.0 0 0
#> 7 0.0 0.0 1 1
#> 8 0.0 0.0 1 1
```
