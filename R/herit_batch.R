#' Batch heritability estimation over multiple traits
#'
#' A convenience wrapper around [herit_vc()] that iterates over a vector of
#' trait names, optionally across multiple covariate models, and returns a
#' tidy data frame.
#'
#' @param traits Character vector of trait column names in `data`.
#' @param grm Numeric matrix: additive genetic relationship matrix as returned
#'   by [build_grm()].
#' @param data Data frame containing ID, trait, and covariate columns.
#' @param covs_list A named list of covariate vectors, where each element
#'   defines one covariate model. If `NULL`, a single unadjusted model is run.
#'   Example: `list(unadj = NULL, cov1 = c("age", "sex"),
#'                  cov2 = c("age", "sex", "age2"))`.
#' @param id_col Name of the individual ID column. Default `"IID"`.
#' @param min_n Minimum sample size to attempt estimation. Default `80`.
#' @param ci_level Profile-likelihood CI level. Default `0.95`.
#' @param .progress Logical. Show a cli progress bar. Default `TRUE`.
#'
#' @return A data frame (tibble-compatible) with one row per successfully
#'   fitted model and columns: `label`, `trait`, `covariates`, `n`, `h2`,
#'   `se`, `ci_lo`, `ci_hi`, `pval`, `sigma2_a`, `sigma2_e`.
#'   Failed / skipped models are silently omitted.
#'
#' @examples
#' \dontrun{
#' A    <- build_grm(my_pedigree, study_ids = my_data$IID)
#'
#' res  <- herit_batch(
#'   traits     = c("bmi", "systolic_bp", "hdl"),
#'   grm        = A,
#'   data       = my_data,
#'   covs_list  = list(
#'     unadj = NULL,
#'     cov1  = c("age", "sex"),
#'     cov2  = c("age", "sex", "age2")
#'   )
#' )
#'
#' # Significant adjusted models
#' subset(res, grepl("cov2", label) & pval < 0.05)
#' }
#'
#' @seealso [herit_vc()], [build_grm()], [plot_forest()]
#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done
#' @importFrom rlang warn
#' @export
herit_batch <- function(traits,
                        grm,
                        data,
                        covs_list  = NULL,
                        id_col     = "IID",
                        min_n      = 80L,
                        ci_level   = 0.95,
                        .progress  = TRUE) {

  if (is.null(covs_list)) covs_list <- list(unadj = NULL)

  # Build a grid of (trait x model) combinations
  grid <- expand.grid(
    trait = traits,
    model = names(covs_list),
    stringsAsFactors = FALSE
  )

  n_jobs <- nrow(grid)
  if (.progress) {
    pb <- cli::cli_progress_bar("Estimating heritability",
                                total = n_jobs, clear = FALSE)
  }

  results <- vector("list", n_jobs)
  for (i in seq_len(n_jobs)) {
    tr    <- grid$trait[i]
    mod   <- grid$model[i]
    covs  <- covs_list[[mod]]
    lbl   <- paste0(tr, "_", mod)

    results[[i]] <- herit_vc(
      trait    = tr,
      grm      = grm,
      data     = data,
      covs     = covs,
      id_col   = id_col,
      label    = lbl,
      min_n    = min_n,
      ci_level = ci_level,
      verbose  = FALSE
    )

    if (.progress) cli::cli_progress_update(id = pb)
  }

  if (.progress) cli::cli_progress_done(id = pb)

  # Bind non-NULL results into a data frame
  ok <- Filter(Negate(is.null), results)
  if (length(ok) == 0L) {
    rlang::warn("All models were skipped (n < min_n or degenerate GRM).")
    return(data.frame())
  }

  out <- do.call(rbind, lapply(ok, as.data.frame))
  rownames(out) <- NULL
  out
}
