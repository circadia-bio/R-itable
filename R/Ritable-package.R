#' @keywords internal
"_PACKAGE"

## Suppress R CMD CHECK notes for rlang .data pronoun
NULL

# ---- Package colour palette -------------------------------------------------

#' R-itable colour palette
#'
#' A named character vector of the four brand colours used throughout
#' R-itable figures and documentation.
#'
#' \describe{
#'   \item{`pink`}{`#FE9EC7` -- primary accent; significant traits, sex effect}
#'   \item{`blue`}{`#44ACFF` -- primary brand colour; estimates, CI lines}
#'   \item{`sky`}{`#89D4FF` -- secondary blue; muted elements, CI shading}
#'   \item{`cream`}{`#F9F6C4` -- soft highlight; code backgrounds}
#' }
#'
#' @examples
#' Ritable_colours
#'
#' @export
Ritable_colours <- c(
  pink  = "#FE9EC7",
  blue  = "#44ACFF",
  sky   = "#89D4FF",
  cream = "#F9F6C4"
)
