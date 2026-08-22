#' Inverse-normal (rank-based) transformation
#'
#' Applies a rank-based inverse-normal transformation (INT) to a numeric
#' vector, placing empirical quantiles onto a standard normal scale.
#'
#' @param x Numeric vector. `NA`s are preserved.
#' @param method Which INT formula to use. `"vdw"` (default) is the
#'   Van der Waerden (1952) transformation, `\Phi^{-1}(r_i / (n+1))`, which
#'   is what SOLAR Eclipse's `inormal` command actually uses. `"blom"` is
#'   the Blom (1958) transformation, `\Phi^{-1}((r_i - 0.5) / n)`, used by
#'   R-itable prior to this fix under the (incorrect) belief that it
#'   matched SOLAR.
#' @param ties Method passed to [base::rank()]. Default `"average"`.
#'
#' @return A numeric vector of the same length as `x`, with `NA`s in the same
#'   positions and non-missing values transformed to approximate normality.
#'
#' @details
#' For `method = "vdw"` (default):
#' \deqn{\Phi^{-1}\!\left(\frac{r_i}{n+1}\right)}
#' For `method = "blom"`:
#' \deqn{\Phi^{-1}\!\left(\frac{r_i - 0.5}{n}\right)}
#' where \eqn{r_i} is the rank of observation \eqn{i} and \eqn{n} is the
#' number of non-missing observations.
#'
#' Van der Waerden is SOLAR Eclipse's actual `inormal` formula (confirmed
#' against SOLAR's documentation and independent technical references); it
#' is the default here so that [herit_vc()], [herit_ace()], and
#' [herit_bivar()] match SOLAR's convention out of the box. `"blom"` is
#' retained as an option for continuity with earlier R-itable versions and
#' for users who specifically want that formula for other reasons. The two
#' formulas agree closely for moderate-to-large n and differ more near the
#' extremes of small samples or heavily tied data.
#'
#' @examples
#' set.seed(1)
#' x <- rexp(200, rate = 0.5)   # right-skewed
#' hist(x,              main = "Raw")
#' hist(int_transform(x), main = "INT (Van der Waerden)")
#'
#' @importFrom stats qnorm
#' @export
int_transform <- function(x, method = c("vdw", "blom"), ties = "average") {
  method <- match.arg(method)
  n <- sum(!is.na(x))
  r <- rank(x, ties.method = ties, na.last = "keep")
  if (method == "vdw") {
    qnorm(r / (n + 1))
  } else {
    qnorm((r - 0.5) / n)
  }
}
