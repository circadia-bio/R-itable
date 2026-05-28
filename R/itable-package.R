#' @keywords internal
"_PACKAGE"

## Suppress R CMD CHECK notes for rlang .data pronoun
## usethis::use_import_from("rlang", ".data") handles this at build time
NULL

# ---- Package colour palette -------------------------------------------------
# Four-colour palette derived from the R-itable brand.
# Used as default discrete colours in plot_forest() and available to users
# for consistent figure styling in downstream analyses.

#' R-itable colour palette
#'
#' A named character vector of the four brand colours used throughout
#' R-itable figures and documentation.
#'
#' \describe{
#'   \item{`pink`}{`#FE9EC7` — primary accent; significant traits, sex effect}
#'   \item{`blue`}{`#44ACFF` — primary brand colour; estimates, CI lines}
#'   \item{`sky`}{`#89D4FF` — secondary blue; muted elements, CI shading}
#'   \item{`cream`}{`#F9F6C4` — soft highlight; code backgrounds}
#' }
#'
#' @examples
#' itable_colours
#' scales::show_col(itable_colours)
#'
#' @export
itable_colours <- c(
  pink  = "#FE9EC7",
  blue  = "#44ACFF",
  sky   = "#89D4FF",
  cream = "#F9F6C4"
)
