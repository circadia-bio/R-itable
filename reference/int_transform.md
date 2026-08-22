# Inverse-normal (rank-based) transformation

Applies a rank-based inverse-normal transformation (INT) to a numeric
vector, placing empirical quantiles onto a standard normal scale.

## Usage

``` r
int_transform(x, method = c("vdw", "blom"), ties = "average")
```

## Arguments

- x:

  Numeric vector. `NA`s are preserved.

- method:

  Which INT formula to use. `"vdw"` (default) is the Van der
  Waerden (1952) transformation, `\Phi^{-1}(r_i / (n+1))`, which is what
  SOLAR Eclipse's `inormal` command actually uses. `"blom"` is the
  Blom (1958) transformation, `\Phi^{-1}((r_i - 0.5) / n)`, used by
  R-itable prior to this fix under the (incorrect) belief that it
  matched SOLAR.

- ties:

  Method passed to [`base::rank()`](https://rdrr.io/r/base/rank.html).
  Default `"average"`.

## Value

A numeric vector of the same length as `x`, with `NA`s in the same
positions and non-missing values transformed to approximate normality.

## Details

For `method = "vdw"` (default):
\$\$\Phi^{-1}\\\left(\frac{r_i}{n+1}\right)\$\$ For `method = "blom"`:
\$\$\Phi^{-1}\\\left(\frac{r_i - 0.5}{n}\right)\$\$ where \\r_i\\ is the
rank of observation \\i\\ and \\n\\ is the number of non-missing
observations.

Van der Waerden is SOLAR Eclipse's actual `inormal` formula (confirmed
against SOLAR's documentation and independent technical references); it
is the default here so that
[`herit_vc()`](https://r-itable.circadia-lab.uk/reference/herit_vc.md),
[`herit_ace()`](https://r-itable.circadia-lab.uk/reference/herit_ace.md),
and
[`herit_bivar()`](https://r-itable.circadia-lab.uk/reference/herit_bivar.md)
match SOLAR's convention out of the box. `"blom"` is retained as an
option for continuity with earlier R-itable versions and for users who
specifically want that formula for other reasons. The two formulas agree
closely for moderate-to-large n and differ more near the extremes of
small samples or heavily tied data.

## Examples

``` r
set.seed(1)
x <- rexp(200, rate = 0.5)   # right-skewed
hist(x,              main = "Raw")

hist(int_transform(x), main = "INT (Van der Waerden)")

```
