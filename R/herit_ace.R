#' A vs C identifiability check (ACE model comparison)
#'
#' Fits four nested variance-component models for a single quantitative
#' trait -- sporadic (no familial effect), polygenic/AE (additive genetic
#' only), household/CE (common environment only), and full ACE (A and C
#' estimated jointly) -- and returns their log-likelihoods together with the
#' two likelihood-ratio tests needed to judge whether A and C are separately
#' identifiable in the sample: CE vs sporadic (is there a familial signal at
#' all?) and ACE vs AE (does C add anything once A is in the model?). This
#' mirrors SOLAR Eclipse's `house` / `polygenic -screen` model-comparison
#' workflow (see `solar_full_analysis.tcl`, Part 3).
#'
#' Unlike [herit_vc()], which diagonalises a single variance component via
#' `eigen()`, this function estimates one or two variance components jointly
#' by direct numerical maximum likelihood over the full `n x n` covariance
#' matrix (Cholesky-based GLS), since `A` and household indicator `C` do not
#' in general share eigenvectors.
#'
#' @param trait Character string: name of the trait column in `data`.
#' @param grm Numeric matrix: additive genetic relationship matrix, as
#'   returned by [build_grm()].
#' @param household Numeric matrix: household indicator matrix, as returned
#'   by [build_household()]. Must have the same row/column names as `grm`.
#' @param data Data frame containing `id_col` and `trait`.
#' @param id_col Name of the individual ID column in `data`. Default `"IID"`.
#' @param transform Logical. Apply [int_transform()] to `trait` before
#'   fitting. Default `TRUE`, matching SOLAR Eclipse convention for skewed
#'   phenotypes.
#' @param min_n Minimum number of complete observations required. Default
#'   `80L`.
#' @param verbose Logical. Print progress via `cli`. Default `TRUE`.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{`trait`}{Trait name.}
#'   \item{`n`}{Sample size after dropping missing values.}
#'   \item{`loglik_sporadic`, `loglik_ae`, `loglik_ce`, `loglik_ace`}{
#'     Log-likelihoods of the four nested models.}
#'   \item{`h2_ae`}{Heritability under the AE model.}
#'   \item{`c2_ce`}{Common-environment proportion under the CE model.}
#'   \item{`h2_ace`, `c2_ace`}{Joint estimates under the full ACE model.}
#'   \item{`chisq_c_vs_sporadic`, `p_c_vs_sporadic`}{One-sided boundary
#'     LRT of CE vs sporadic (evidence of a familial effect at all).}
#'   \item{`chisq_ace_vs_ae`, `p_ace_vs_ae`}{One-sided boundary LRT of ACE
#'     vs AE (evidence C adds anything once A is estimated). A near-zero
#'     statistic with `c2_ace` at or near the 0 boundary indicates A and C
#'     are not jointly identifiable in this sample -- report the simpler AE
#'     model (see package vignette / Leocadio-Miguel et al. companion
#'     analyses).}
#' }
#' Returns `NULL` if `n < min_n`.
#'
#' @details
#' **Model:** `Omega = sigma2_p * [h2*A + c2*C + (1 - h2 - c2)*I]`, jointly
#' estimated for the ACE model via multi-start Nelder-Mead on a
#' logit-reparametrised `(h2, c2)` (guaranteeing `h2, c2 >= 0` and
#' `h2 + c2 < 1`); AE and CE are 1-D special cases.
#'
#' **LRTs:** one-sided chi-squared(1) boundary correction, as in
#' [herit_vc()], since both null hypotheses (`c2 = 0`, or `h2 = 0`) sit on
#' the boundary of the parameter space.
#'
#' @seealso [build_grm()], [build_household()], [herit_vc()], [herit_bivar()]
#' @importFrom stats complete.cases
#' @importFrom rlang abort
#' @importFrom cli cli_alert_success cli_alert_warning
#' @export
herit_ace <- function(trait,
                      grm,
                      household,
                      data,
                      id_col    = "IID",
                      transform = TRUE,
                      min_n     = 80L,
                      verbose   = TRUE) {

  needed <- unique(c(id_col, trait))
  absent <- setdiff(needed, names(data))
  if (length(absent)) {
    rlang::abort(c("Column(s) not found in `data`:", paste0("  ", paste(absent, collapse = ", "))))
  }
  if (!identical(rownames(grm), rownames(household))) {
    rlang::abort("`grm` and `household` must have identical row/column names, in the same order.")
  }

  dat <- data[stats::complete.cases(data[, needed, drop = FALSE]), needed, drop = FALSE]
  n   <- nrow(dat)
  if (n < min_n) {
    if (verbose) cli::cli_alert_warning("Skipping {trait}: n = {n} < {min_n}.")
    return(NULL)
  }

  ids <- as.character(dat[[id_col]])
  if (!all(ids %in% rownames(grm))) {
    rlang::abort("Some `data` IDs are absent from `grm`/`household`.")
  }
  A <- grm[ids, ids]
  C <- household[ids, ids]

  y <- if (transform) int_transform(dat[[trait]]) else dat[[trait]]
  X <- matrix(1, n, 1)

  m_sp  <- .fit_ace_uni(y, X, A, C, "none")
  m_ae  <- .fit_ace_uni(y, X, A, C, "A")
  m_ce  <- .fit_ace_uni(y, X, A, C, "C")
  m_ace <- .fit_ace_uni(y, X, A, C, "AC")

  p_c_vs_sp  <- .lrt_boundary(m_ce$loglik,  m_sp$loglik, df = 1)
  p_ace_vs_ae <- .lrt_boundary(m_ace$loglik, m_ae$loglik, df = 1)

  if (verbose) {
    cli::cli_alert_success(
      "{trait}  n={n}  h2(AE)={round(m_ae$h2,3)}  c2(CE)={round(m_ce$c2,3)}  c2(ACE)={round(m_ace$c2,4)}  p(ACE vs AE)={signif(p_ace_vs_ae,3)}"
    )
  }

  list(
    trait               = trait,
    n                   = n,
    loglik_sporadic     = round(m_sp$loglik, 5),
    loglik_ae           = round(m_ae$loglik, 5),
    loglik_ce           = round(m_ce$loglik, 5),
    loglik_ace          = round(m_ace$loglik, 5),
    h2_ae               = round(m_ae$h2, 4),
    c2_ce               = round(m_ce$c2, 4),
    h2_ace              = round(m_ace$h2, 4),
    c2_ace              = round(m_ace$c2, 4),
    chisq_c_vs_sporadic = round(max(0, 2 * (m_ce$loglik - m_sp$loglik)), 4),
    p_c_vs_sporadic     = signif(p_c_vs_sp, 5),
    chisq_ace_vs_ae     = round(max(0, 2 * (m_ace$loglik - m_ae$loglik)), 4),
    p_ace_vs_ae         = signif(p_ace_vs_ae, 5)
  )
}
