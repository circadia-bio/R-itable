#' Inverse-normal (rank-based) transformation
#'
#' Applies a rank-based inverse-normal transformation (INT) to a numeric
#' vector, placing empirical quantiles onto a standard normal scale.
#'
#' @param x Numeric vector. `NA`s are preserved.
#' @param ties Method passed to [base::rank()]. Default `"average"`.
#'
#' @return A numeric vector of the same length as `x`, with `NA`s in the same
#'   positions and non-missing values transformed to approximate normality.
#'
#' @details
#' The Blom-style formula used is:
#' \deqn{\Phi^{-1}\!\left(\frac{r_i - 0.5}{n}\right)}
#' where \eqn{r_i} is the rank of observation \eqn{i} and \eqn{n} is the
#' number of non-missing observations. This is the standard transformation used
#' in SOLAR Eclipse and most variance-components heritability software.
#'
#' @examples
#' set.seed(1)
#' x <- rexp(200, rate = 0.5)   # right-skewed
#' hist(x,              main = "Raw")
#' hist(int_transform(x), main = "INT")
#'
#' @importFrom stats qnorm
#' @export
int_transform <- function(x, ties = "average") {
  n    <- sum(!is.na(x))
  r    <- rank(x, ties.method = ties, na.last = "keep")
  qnorm((r - 0.5) / n)
}
