# Forest plot of heritability estimates

Produces a ggplot2 forest plot from the output of
[`herit_batch()`](https://r-itable.circadia-lab.uk/reference/herit_batch.md)
or a list of
[`herit_vc()`](https://r-itable.circadia-lab.uk/reference/herit_vc.md)
results coerced to a data frame.

## Usage

``` r
plot_forest(
  results,
  model_filter = NULL,
  colour_by = "trait",
  sig_threshold = 0.05,
  title = NULL,
  x_limits = c(0, 1)
)
```

## Arguments

- results:

  Data frame as returned by
  [`herit_batch()`](https://r-itable.circadia-lab.uk/reference/herit_batch.md),
  containing at least columns `label`, `trait`, `h2`, `ci_lo`, `ci_hi`,
  `pval`.

- model_filter:

  Optional character vector of model name substrings to keep (matched
  against `label`). E.g. `"cov2"` to show only the age + sex +
  age-squared model.

- colour_by:

  Column to colour points by. Default `"trait"`. Set to `NULL` for a
  monochrome plot.

- sig_threshold:

  Numeric. Traits with `pval` below this threshold are shown with a
  filled point; others with an open point. Default `0.05`.

- title:

  Optional plot title string.

- x_limits:

  Numeric vector of length 2 for the x-axis range. Default `c(0, 1)`.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## Details

Requires **ggplot2** (listed in `Suggests`). An informative error is
thrown if it is not installed.

## See also

[`herit_batch()`](https://r-itable.circadia-lab.uk/reference/herit_batch.md)

## Examples

``` r
if (FALSE) { # \dontrun{
res <- herit_batch(c("bmi", "hdl", "systolic_bp"),
                   grm = A, data = my_data,
                   covs_list = list(unadj = NULL,
                                    cov2  = c("age", "sex", "age2")))

plot_forest(res, model_filter = "cov2", title = "Adjusted heritability")
} # }
```
