#' Inverse-normal (rank-based) transformation
#'
#' Applies a rank-based inverse-normal transformation (INT) to a numeric
#' vector, placing empirical quantiles onto a standard normal scale.
#'
#' @param x Numeric vector. `NA`s are preserved.
#' @param method Which INT formula to use. `"vdw"` (default) exactly
#'   replicates SOLAR Eclipse's `inormal` command, confirmed against SOLAR's
#'   actual source (`lib/solar.tcl`, `proc inormal`): sort ascending, assign
#'   `z_i = qnorm(i / (n+1))` to each individual rank position `i`, and for
#'   tied values average the resulting *z-scores* of their individual rank
#'   positions (not the z-score of the averaged rank -- these differ for
#'   nonlinear qnorm whenever a tied block is not symmetric around n/2,
#'   e.g. a large block of tied values at one tail, such as zero-inflated
#'   phenotypes). `"blom"` is the Blom (1958) transformation,
#'   `qnorm((r-0.5)/n)` with conventional averaged-rank tie handling, used
#'   by R-itable prior to this fix under the (incorrect) belief that it
#'   matched SOLAR.
#' @param ties Method passed to [base::rank()] when `method = "blom"`.
#'   Ignored for `method = "vdw"`, which always uses SOLAR's own tie
#'   handling (see Details). Default `"average"`.
#'
#' @return A numeric vector of the same length as `x`, with `NA`s in the same
#'   positions and non-missing values transformed to approximate normality.
#'
#' @details
#' For `method = "vdw"` with no ties, this is the Van der Waerden (1952)
#' transformation, `qnorm(r_i / (n+1))`, and is bit-identical to it. With
#' ties, tied individuals get the mean of `qnorm(i/(n+1))` over the tied
#' rank positions `i`, exactly matching SOLAR's `inormal`. This distinction
#' only matters for heavily-tied data; it was confirmed to fully resolve
#' the remaining discrepancy between R-itable and SOLAR Eclipse on
#' zero-inflated phenotypes in the Gomes et al. validation dataset (\eqn{>}
#' 25% tied values) -- see NEWS.md.
#'
#' For `method = "blom"`:
#' \deqn{\Phi^{-1}\!\left(\frac{r_i - 0.5}{n}\right)}
#' using [base::rank()]'s `ties` argument for tie handling (default
#' `"average"`, i.e. the z-score of the averaged rank -- the conventional
#' approach for this formula, not SOLAR's).
#'
#' @examples
#' set.seed(1)
#' x <- rexp(200, rate = 0.5)   # right-skewed
#' hist(x,              main = "Raw")
#' hist(int_transform(x), main = "INT (SOLAR-matching Van der Waerden)")
#'
#' # Heavy ties (e.g. zero-inflated activity duration): SOLAR's tie handling
#' # (default) differs meaningfully from naive rank-average tie handling.
#' y <- c(rep(0, 40), round(runif(60, 1, 100)))
#' int_transform(y)[1]                        # SOLAR-exact
#' qnorm(rank(y, ties.method = "average")[1] / (length(y) + 1))  # naive
#'
#' @importFrom stats qnorm
#' @export
int_transform <- function(x, method = c("vdw", "blom"), ties = "average") {
  method <- match.arg(method)

  if (method == "blom") {
    n <- sum(!is.na(x))
    r <- rank(x, ties.method = ties, na.last = "keep")
    return(qnorm((r - 0.5) / n))
  }

  # method == "vdw": exact replication of SOLAR's inormal tie handling.
  n_total <- length(x)
  na_idx  <- is.na(x)
  n       <- sum(!na_idx)
  out     <- rep(NA_real_, n_total)
  if (n == 0) return(out)

  ord   <- order(x[!na_idx])
  xs    <- x[!na_idx][ord]
  z_all <- qnorm(seq_len(n) / (n + 1))

  run    <- rle(xs)
  ends   <- cumsum(run$lengths)
  starts <- ends - run$lengths + 1
  z_by_sorted_pos <- numeric(n)
  for (k in seq_along(run$lengths)) {
    z_by_sorted_pos[starts[k]:ends[k]] <- mean(z_all[starts[k]:ends[k]])
  }
  out[!na_idx][ord] <- z_by_sorted_pos
  out
}
