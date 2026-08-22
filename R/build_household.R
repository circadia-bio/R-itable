#' Build a household (shared-environment) indicator matrix
#'
#' Constructs the block-diagonal indicator matrix **C** used as the shared
#' /common-environment ("household") random-effect design in [herit_ace()]:
#' `C[i, j] = 1` if individuals `i` and `j` belong to the same household
#' group, `0` otherwise (including the diagonal convention `C[i, i] = 1`).
#'
#' @param hhid Vector of household group labels, one per individual, in the
#'   same order as `ids`. Individuals sharing an `hhid` value are treated as
#'   co-resident (e.g. twin pairs reared together).
#' @param ids Vector of individual IDs, the same length as `hhid`. Used to
#'   name the rows/columns of the returned matrix so it aligns with a GRM
#'   from [build_grm()].
#'
#' @return A symmetric numeric matrix of dimension `length(ids) x
#'   length(ids)`, with row/column names `as.character(ids)`.
#'
#' @examples
#' hh  <- c("A", "A", "B", "B", "C")
#' ids <- 1:5
#' build_household(hh, ids)
#'
#' @export
build_household <- function(hhid, ids) {
  if (length(hhid) != length(ids)) {
    rlang::abort("`hhid` and `ids` must be the same length.")
  }
  ids <- as.character(ids)
  n   <- length(ids)
  C   <- matrix(0, n, n, dimnames = list(ids, ids))
  for (h in unique(hhid)) {
    idx <- which(hhid == h)
    C[idx, idx] <- 1
  }
  C
}
