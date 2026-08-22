#' Bivariate genetic and environmental correlation between two traits
#'
#' Fits a bivariate AE variance-components model for two quantitative traits
#' and estimates the additive-genetic correlation (rhoG), the non-shared
#' environmental correlation (rhoE), and the derived phenotypic correlation
#' (rhoP = rhoG*sqrt(h2_1*h2_2) + rhoE*sqrt((1-h2_1)*(1-h2_2))), together
#' with likelihood-ratio tests of rhoG = 0, rhoE = 0, and the joint test of
#' rhoP = 0 (rhoG = rhoE = 0). This mirrors SOLAR Eclipse's
#' `polygenic -testrhog` / `-testrhoe` bivariate output (see
#' `solar_full_analysis.tcl`, Part 5).
#'
#' @param trait1,trait2 Character strings: names of the two trait columns in
#'   `data`.
#' @param grm Numeric matrix: additive genetic relationship matrix, as
#'   returned by [build_grm()].
#' @param data Data frame containing `id_col`, `trait1`, and `trait2`.
#' @param id_col Name of the individual ID column in `data`. Default `"IID"`.
#' @param transform Logical. Apply [int_transform()] to both traits before
#'   fitting. Default `TRUE`.
#' @param min_n Minimum number of complete observations (non-missing on both
#'   traits) required. Default `80L`.
#' @param verbose Logical. Print progress via `cli`. Default `TRUE`.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{`trait1`, `trait2`}{Trait names.}
#'   \item{`n`}{Sample size after dropping missing values.}
#'   \item{`h2_1`, `h2_2`}{Heritability of each trait under the joint
#'     bivariate model.}
#'   \item{`rhoG`, `rhoE`, `rhoP`}{Genetic, environmental, and derived
#'     phenotypic correlations.}
#'   \item{`p_rhoG`, `p_rhoE`}{LRT p-values for rhoG = 0 and rhoE = 0
#'     (standard two-sided chi-squared(1); these correlations are not
#'     boundary parameters).}
#'   \item{`p_rhoP`}{LRT p-value for the joint null rhoG = rhoE = 0
#'     (chi-squared(2)).}
#' }
#' Returns `NULL` if `n < min_n`.
#'
#' @details
#' **Model:** traits are stacked into a `2n`-vector with block covariance
#' `Sigma = [[sigma2_g1*A + sigma2_e1*I, rhoG*sqrt(sigma2_g1*sigma2_g2)*A +
#' rhoE*sqrt(sigma2_e1*sigma2_e2)*I], [..., sigma2_g2*A + sigma2_e2*I]]`,
#' where `sigma2_g_i = h2_i * sigma2_p_i` and `sigma2_e_i = (1-h2_i) *
#' sigma2_p_i`. All six free parameters (`h2_1, h2_2, rhoG, rhoE, sigma2_p1,
#' sigma2_p2`) are estimated jointly by multi-start Nelder-Mead ML on a
#' reparametrised scale (logit for `h2`, `atanh` for the correlations, `log`
#' for the variances) so that all constraints hold automatically.
#'
#' As with [herit_ace()], this uses direct numerical ML over the full `2n x
#' 2n` covariance matrix rather than `herit_vc()`'s eigendecomposition
#' shortcut, since the cross-trait blocks do not diagonalise jointly with the
#' within-trait blocks in general. In practice this is optimised by
#' diagonalising `A` once and reducing every likelihood evaluation to `n`
#' independent 2x2 blocks (both covariance blocks are linear combinations of
#' `A` and `I`, so they share `A`'s eigenvectors) -- this is what makes
#' [herit_bivar_batch()] over many pairs fast.
#'
#' Results are validated against, but not bit-identical to, SOLAR Eclipse
#' output (differences arise from ML vs SOLAR's exact optimiser and variance
#' parametrisation). Validated on the Gomes et al. twin dataset: mean
#' |delta rhoG| ~ 0.05, mean |delta rhoE| ~ 0.02 vs SOLAR, with matching
#' nominal-significance calls on 29/31 and 31/31 of 31 trait pairs
#' respectively (see NEWS.md).
#'
#' @seealso [build_grm()], [herit_vc()], [herit_ace()], [herit_bivar_batch()]
#' @importFrom stats complete.cases pchisq
#' @importFrom rlang abort
#' @importFrom cli cli_alert_success cli_alert_warning
#' @export
herit_bivar <- function(trait1,
                        trait2,
                        grm,
                        data,
                        id_col    = "IID",
                        transform = TRUE,
                        min_n     = 80L,
                        verbose   = TRUE) {

  needed <- unique(c(id_col, trait1, trait2))
  absent <- setdiff(needed, names(data))
  if (length(absent)) {
    rlang::abort(c("Column(s) not found in `data`:", paste0("  ", paste(absent, collapse = ", "))))
  }

  dat <- data[stats::complete.cases(data[, needed, drop = FALSE]), needed, drop = FALSE]
  n   <- nrow(dat)
  if (n < min_n) {
    if (verbose) cli::cli_alert_warning("Skipping {trait1} x {trait2}: n = {n} < {min_n}.")
    return(NULL)
  }

  ids <- as.character(dat[[id_col]])
  if (!all(ids %in% rownames(grm))) {
    rlang::abort("Some `data` IDs are absent from `grm`.")
  }
  A <- grm[ids, ids]

  y1 <- if (transform) int_transform(dat[[trait1]]) else dat[[trait1]]
  y2 <- if (transform) int_transform(dat[[trait2]]) else dat[[trait2]]

  full <- .fit_bivar_full(y1, y2, A)

  ll_rhoG0 <- .fit_bivar_constrained(y1, y2, A, fix_rhoG = TRUE,  start = full$par, prep = full$prep)
  ll_rhoE0 <- .fit_bivar_constrained(y1, y2, A, fix_rhoE = TRUE,  start = full$par, prep = full$prep)
  ll_00    <- .fit_bivar_constrained(y1, y2, A, fix_rhoG = TRUE, fix_rhoE = TRUE,
                                    start = full$par, prep = full$prep)

  p_rhoG <- stats::pchisq(max(0, 2 * (full$loglik - ll_rhoG0)), df = 1, lower.tail = FALSE)
  p_rhoE <- stats::pchisq(max(0, 2 * (full$loglik - ll_rhoE0)), df = 1, lower.tail = FALSE)
  p_rhoP <- stats::pchisq(max(0, 2 * (full$loglik - ll_00)),    df = 2, lower.tail = FALSE)

  rhoP <- full$rhoG * sqrt(full$h2_1 * full$h2_2) +
    full$rhoE * sqrt((1 - full$h2_1) * (1 - full$h2_2))

  if (verbose) {
    cli::cli_alert_success(
      "{trait1} x {trait2}  n={n}  rhoG={round(full$rhoG,3)} (p={signif(p_rhoG,3)})  rhoE={round(full$rhoE,3)} (p={signif(p_rhoE,3)})"
    )
  }

  list(
    trait1 = trait1, trait2 = trait2, n = n,
    h2_1   = round(full$h2_1, 4), h2_2 = round(full$h2_2, 4),
    rhoG   = round(full$rhoG, 4), rhoE = round(full$rhoE, 4), rhoP = round(rhoP, 4),
    p_rhoG = signif(p_rhoG, 5), p_rhoE = signif(p_rhoE, 5), p_rhoP = signif(p_rhoP, 5)
  )
}

#' Batch bivariate genetic/environmental correlations with FDR correction
#'
#' Iterates [herit_bivar()] over a set of trait pairs and applies
#' Benjamini-Hochberg FDR correction separately within each of the three
#' p-value families (rhoP, rhoG, rhoE), matching `parse_and_fdr.py`'s
#' post-processing of SOLAR output.
#'
#' @param pairs A two-column character matrix or data frame of trait name
#'   pairs (column 1 = trait1, column 2 = trait2), or a list of length-2
#'   character vectors.
#' @inheritParams herit_bivar
#' @param .progress Logical. Show a cli progress bar. Default `TRUE`.
#'
#' @return A data frame with one row per successfully fitted pair and columns
#'   `trait1`, `trait2`, `n`, `h2_1`, `h2_2`, `rhoG`, `rhoE`, `rhoP`,
#'   `p_rhoG`, `p_rhoE`, `p_rhoP`, `q_rhoG`, `q_rhoE`, `q_rhoP` (BH-adjusted).
#'   Failed / skipped pairs are silently omitted.
#'
#' @seealso [herit_bivar()]
#' @importFrom stats p.adjust
#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done
#' @importFrom rlang warn
#' @export
herit_bivar_batch <- function(pairs,
                              grm,
                              data,
                              id_col    = "IID",
                              transform = TRUE,
                              min_n     = 80L,
                              .progress = TRUE) {

  pairs <- as.data.frame(pairs, stringsAsFactors = FALSE)
  names(pairs)[1:2] <- c("trait1", "trait2")
  n_jobs <- nrow(pairs)

  if (.progress) {
    pb <- cli::cli_progress_bar("Estimating bivariate correlations", total = n_jobs, clear = FALSE)
  }

  results <- vector("list", n_jobs)
  for (i in seq_len(n_jobs)) {
    results[[i]] <- herit_bivar(
      trait1    = pairs$trait1[i],
      trait2    = pairs$trait2[i],
      grm       = grm,
      data      = data,
      id_col    = id_col,
      transform = transform,
      min_n     = min_n,
      verbose   = FALSE
    )
    if (.progress) cli::cli_progress_update(id = pb)
  }
  if (.progress) cli::cli_progress_done(id = pb)

  ok <- Filter(Negate(is.null), results)
  if (length(ok) == 0L) {
    rlang::warn("All pairs were skipped (n < min_n).")
    return(data.frame())
  }

  out <- do.call(rbind, lapply(ok, as.data.frame))
  rownames(out) <- NULL

  out$q_rhoG <- stats::p.adjust(out$p_rhoG, method = "BH")
  out$q_rhoE <- stats::p.adjust(out$p_rhoE, method = "BH")
  out$q_rhoP <- stats::p.adjust(out$p_rhoP, method = "BH")

  out
}
