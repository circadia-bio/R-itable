# Inverse-normal (rank-based) transformation

Applies a rank-based inverse-normal transformation (INT) to a numeric
vector, placing empirical quantiles onto a standard normal scale.

## Usage

``` r
int_transform(x, ties = "average")
```

## Arguments

- x:

  Numeric vector. `NA`s are preserved.

- ties:

  Method passed to [`base::rank()`](https://rdrr.io/r/base/rank.html).
  Default `"average"`.

## Value

A numeric vector of the same length as `x`, with `NA`s in the same
positions and non-missing values transformed to approximate normality.

## Details

The Blom-style formula used is: \$\$\Phi^{-1}\\\left(\frac{r_i -
0.5}{n}\right)\$\$ where \\r_i\\ is the rank of observation \\i\\ and
\\n\\ is the number of non-missing observations. This is the standard
transformation used in SOLAR Eclipse and most variance-components
heritability software.

## Examples

``` r
set.seed(1)
x <- rexp(200, rate = 0.5)   # right-skewed
hist(x,              main = "Raw")

hist(int_transform(x), main = "INT")

```
