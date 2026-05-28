#' Forest plot of heritability estimates
#'
#' Produces a ggplot2 forest plot from the output of [herit_batch()] or a
#' list of [herit_vc()] results coerced to a data frame.
#'
#' @param results Data frame as returned by [herit_batch()], containing at
#'   least columns `label`, `trait`, `h2`, `ci_lo`, `ci_hi`, `pval`.
#' @param model_filter Optional character vector of model name substrings to
#'   keep (matched against `label`). E.g. `"cov2"` to show only the
#'   age + sex + age-squared model.
#' @param colour_by Column to colour points by. Default `"trait"`. Set to
#'   `NULL` for a monochrome plot.
#' @param sig_threshold Numeric. Traits with `pval` below this threshold are
#'   shown with a filled point; others with an open point. Default `0.05`.
#' @param title Optional plot title string.
#' @param x_limits Numeric vector of length 2 for the x-axis range.
#'   Default `c(0, 1)`.
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @details
#' Requires **ggplot2** (listed in `Suggests`). An informative error is thrown
#' if it is not installed.
#'
#' @examples
#' \dontrun{
#' res <- herit_batch(c("bmi", "hdl", "systolic_bp"),
#'                    grm = A, data = my_data,
#'                    covs_list = list(unadj = NULL,
#'                                     cov2  = c("age", "sex", "age2")))
#'
#' plot_forest(res, model_filter = "cov2", title = "Adjusted heritability")
#' }
#'
#' @seealso [herit_batch()]
#' @importFrom rlang .data abort
#' @export
plot_forest <- function(results,
                        model_filter  = NULL,
                        colour_by     = "trait",
                        sig_threshold = 0.05,
                        title         = NULL,
                        x_limits      = c(0, 1)) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    rlang::abort(c(
      "`plot_forest()` requires the ggplot2 package.",
      "i" = 'Install it with: install.packages("ggplot2")'
    ))
  }

  df <- results

  # Optional model filter
  if (!is.null(model_filter)) {
    keep <- vapply(model_filter, function(f) grepl(f, df$label), logical(nrow(df)))
    df   <- df[rowSums(keep) > 0, , drop = FALSE]
    if (nrow(df) == 0L)
      rlang::abort("No rows remaining after applying `model_filter`.")
  }

  df$sig   <- df$pval < sig_threshold
  df$label <- factor(df$label, levels = rev(unique(df$label)))

  # Base aesthetics
  aes_call <- if (!is.null(colour_by) && colour_by %in% names(df)) {
    ggplot2::aes(
      x      = .data[["h2"]],
      y      = .data[["label"]],
      colour = .data[[colour_by]],
      shape  = .data[["sig"]]
    )
  } else {
    ggplot2::aes(
      x     = .data[["h2"]],
      y     = .data[["label"]],
      shape = .data[["sig"]]
    )
  }

  ggplot2::ggplot(df, aes_call) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.4) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data[["ci_lo"]], xmax = .data[["ci_hi"]]),
      orientation = "y", width = 0.3, linewidth = 0.7
    ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_shape_manual(
      values = c(`FALSE` = 1L, `TRUE` = 16L),
      labels = c(`FALSE` = paste0("p >= ", sig_threshold),
                 `TRUE`  = paste0("p < ",  sig_threshold)),
      name   = NULL
    ) +
    ggplot2::scale_x_continuous(
      limits = x_limits,
      breaks = seq(x_limits[1], x_limits[2], by = 0.2)
    ) +
    ggplot2::labs(
      x     = "Narrow-sense heritability (h\u00b2)",
      y     = NULL,
      title = title
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor  = ggplot2::element_blank(),
      panel.grid.major  = ggplot2::element_line(colour = "grey92", linewidth = 0.3),
      legend.background = ggplot2::element_blank(),
      legend.position   = "bottom"
    )
}
