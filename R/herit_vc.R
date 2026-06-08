# Internal: profile log-likelihood for a single h2 value.
# Fits the mixed model Omega = sigma2_p * [h2*A + (1-h2)*I_n] via
# the eigendecomposition Vt already computed outside the optimiser loop.
.profile_ll <- function(h2v, y_int, Vt, ev, X) {
  n     <- length(y_int)
  denom <- h2v * ev + (1 - h2v)
  yr    <- Vt %*% y_int
  Xr    <- Vt %*% X
  w     <- 1 / denom
  beta  <- tryCatch(
    solve(t(Xr) %*% (w * Xr), t(Xr) %*% (w * yr)),
    error = function(e) NULL
  )
  if (is.null(beta)) return(-1e10)
  r  <- yr - Xr %*% beta
  sp <- sum(w * r^2) / n
  -0.5 * (sum(log(denom)) + n * log(sp) + n)
}


#' Profile-likelihood variance-components heritability estimation
#'
#' Estimates narrow-sense heritability (h\ifelse{html}{\out{&sup2;}}{\eqn{^2}})
#' for a single quantitative trait using a profile-likelihood variance-components
#' approach equivalent to SOLAR Eclipse. The phenotype is inverse-normal
#' transformed internally.
#'
#' @param trait Character string: name of the trait column in `data`.
#' @param grm Numeric matrix: the additive genetic relationship matrix for
#'   all individuals in `data`, as returned by [build_grm()]. Row and column
#'   names must match the ID column of `data`.
#' @param data Data frame containing `id_col`, `trait`, and any covariate
#'   columns.
#' @param covs Character vector of covariate column names, or `NULL` for an
#'   intercept-only model. Covariates are mean-centred and scaled to unit
#'   variance before fitting.
#' @param id_col Name of the individual ID column in `data`. Default `"IID"`.
#' @param label Optional string label for this model (used in batch output).
#'   If `NULL`, defaults to `"<trait>_adj"` or `"<trait>_unadj"`.
#' @param min_n Minimum number of complete observations required to attempt
#'   estimation. Models with fewer observations return `NULL` silently.
#'   Default `80`.
#' @param ci_level Confidence level for the profile-likelihood interval.
#'   Default `0.95`.
#' @param verbose Logical. Print progress to the console. Default `TRUE`.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{`label`}{Model label.}
#'   \item{`trait`}{Trait name.}
#'   \item{`covariates`}{Covariate names joined by `"+"`, or `""`.}
#'   \item{`n`}{Sample size after dropping missing values.}
#'   \item{`h2`}{MLE of narrow-sense heritability.}
#'   \item{`se`}{Standard error from profile-likelihood curvature (Wald).}
#'   \item{`ci_lo`, `ci_hi`}{Profile-likelihood confidence interval bounds.}
#'   \item{`pval`}{One-sided LRT p-value with chi-squared(1) boundary correction.}
#'   \item{`var_covariates`}{Proportion of phenotypic variance explained by
#'     fixed-effect covariates (R² on INT-transformed phenotype). `NA` for
#'     unadjusted models. Corresponds to the "variance explained" column in
#'     Leocadio-Miguel et al. (2025).}
#'   \item{`sigma2_a`}{Additive genetic variance (sigma²_g in SOLAR notation).}
#'   \item{`sigma2_e`}{Residual environmental variance (sigma²_e).}
#' }
#' Returns `NULL` if `n < min_n` or if the GRM subset is degenerate.
#'
#' @details
#' **Model:** Omega = sigma2_p \[h2 A + (1 - h2) I_n\]
#'
#' **Optimisation:** eigendecomposition of A followed by 1-D
#' profile-likelihood optimisation over h2 in (0, 1) using a coarse grid
#' to seed [stats::optimize()].
#'
#' **LRT:** one-sided chi-squared(1) boundary correction -- the null is on the
#' boundary of the parameter space (h2 = 0), so the p-value is halved relative
#' to a standard chi-squared test. This matches the SOLAR Eclipse convention.
#'
#' **CIs:** derived by `uniroot()` on the profile log-likelihood, falling
#' back to a Wald interval if the root-finding fails.
#'
#' @examples
#' \dontrun{
#' # Build GRM then estimate heritability
#' A   <- build_grm(my_pedigree, study_ids = my_data$IID)
#' res <- herit_vc("bmi", grm = A, data = my_data,
#'                 covs = c("age", "sex", "age2"))
#' str(res)
#' }
#'
#' @seealso [build_grm()], [herit_batch()], [int_transform()]
#' @importFrom stats optimize uniroot pchisq qchisq sd complete.cases lm residuals
#' @importFrom rlang abort warn .data
#' @importFrom cli cli_alert_success cli_alert_warning
#' @export
herit_vc <- function(trait,
                     grm,
                     data,
                     covs      = NULL,
                     id_col    = "IID",
                     label     = NULL,
                     min_n     = 80L,
                     ci_level  = 0.95,
                     verbose   = TRUE) {

  # -- Column checks ----------------------------------------------------------
  needed <- unique(c(id_col, trait, covs))
  absent <- setdiff(needed, names(data))
  if (length(absent)) {
    rlang::abort(c(
      "Column(s) not found in `data`:",
      paste0("  ", paste(absent, collapse = ", "))
    ))
  }

  # Drop zero-variance covariates (e.g. sex in sex-stratified runs)
  if (!is.null(covs)) {
    tmp <- data[complete.cases(data[, needed, drop = FALSE]), , drop = FALSE]
    bad <- covs[vapply(covs, function(cv) sd(tmp[[cv]], na.rm = TRUE) < 1e-10, logical(1))]
    if (length(bad)) {
      rlang::warn(paste0("Dropping zero-variance covariate(s): ",
                         paste(bad, collapse = ", ")))
      covs  <- setdiff(covs, bad)
      needed <- unique(c(id_col, trait, covs))
    }
  }

  dat <- data[complete.cases(data[, needed, drop = FALSE]), needed, drop = FALSE]
  n   <- nrow(dat)
  tag <- label %||% paste0(trait, ifelse(is.null(covs), "_unadj", "_adj"))

  if (n < min_n) {
    if (verbose) cli::cli_alert_warning("Skipping {tag}: n = {n} < {min_n}.")
    return(NULL)
  }

  # -- GRM subset -------------------------------------------------------------
  ids <- as.character(dat[[id_col]])
  if (!all(ids %in% rownames(grm))) {
    missing_ids <- setdiff(ids, rownames(grm))
    rlang::abort(c(
      paste0(length(missing_ids), " ID(s) in `data` are absent from `grm`."),
      "i" = "Ensure `grm` was built from a pedigree containing all study subjects."
    ))
  }
  A_s <- grm[ids, ids]

  # -- Eigendecomposition -----------------------------------------------------
  eig <- tryCatch(eigen(A_s, symmetric = TRUE),
                  error = function(e) NULL)
  if (is.null(eig)) {
    rlang::warn(paste0("Eigendecomposition failed for ", tag, ". Skipping."))
    return(NULL)
  }
  Vt <- t(eig$vectors)
  ev <- eig$values

  # -- INT-transform phenotype ------------------------------------------------
  y_int <- int_transform(dat[[trait]])

  # -- Design matrix ----------------------------------------------------------
  X <- if (is.null(covs)) {
    matrix(1, n)
  } else {
    cv <- as.matrix(dat[, covs, drop = FALSE])
    for (j in seq_len(ncol(cv))) {
      s <- sd(cv[, j])
      if (s > 0) cv[, j] <- (cv[, j] - mean(cv[, j])) / s
    }
    cbind(1, cv)
  }
  X <- matrix(as.numeric(X), nrow = n)

  # -- Profile-likelihood optimisation ----------------------------------------
  prof    <- function(h2) .profile_ll(h2, y_int, Vt, ev, X)
  h2_grid <- seq(0.001, 0.999, by = 0.005)
  hi      <- h2_grid[which.max(vapply(h2_grid, prof, numeric(1)))]
  fit     <- optimize(prof,
                      interval = c(max(0.001, hi - 0.15), min(0.999, hi + 0.15)),
                      maximum  = TRUE,
                      tol      = 1e-8)
  h2_hat  <- fit$maximum

  # -- Variance components ----------------------------------------------------
  denom <- h2_hat * ev + (1 - h2_hat)
  w     <- 1 / denom
  Xr    <- Vt %*% X
  yr    <- Vt %*% y_int
  bh    <- solve(t(Xr) %*% (w * Xr), t(Xr) %*% (w * yr))
  r_vec <- yr - Xr %*% bh
  sp    <- sum(w * r_vec^2) / n
  sa    <- h2_hat * sp
  se2   <- (1 - h2_hat) * sp

  # -- Variance explained by covariates (R2 of fixed effects on INT y) --------
  var_covariates <- if (is.null(covs)) {
    NA_real_
  } else {
    # Fit intercept-only vs full fixed-effects model on INT y
    ss_tot <- sum((y_int - mean(y_int))^2)
    fit_ols <- lm(y_int ~ X[, -1])  # X already has intercept in col 1
    ss_res  <- sum(residuals(fit_ols)^2)
    round(1 - ss_res / ss_tot, 4)
  }

  # -- LRT (one-sided boundary correction) ------------------------------------
  lrt  <- max(0, 2 * (fit$objective - prof(1e-6)))
  pval <- 0.5 * pchisq(lrt, df = 1, lower.tail = FALSE)

  # -- SE from profile-LL curvature -------------------------------------------
  h2_se <- tryCatch({
    eps <- 1e-4
    d2  <- (prof(h2_hat + eps) - 2 * prof(h2_hat) + prof(h2_hat - eps)) / eps^2
    if (d2 < 0) 1 / sqrt(-d2) else NA_real_
  }, error = function(e) NA_real_)

  # -- Profile-likelihood CI --------------------------------------------------
  ll_thr <- fit$objective - 0.5 * qchisq(ci_level, df = 1)
  ci_lo  <- tryCatch(
    uniroot(function(h) prof(h) - ll_thr, c(0.001, h2_hat), tol = 1e-6)$root,
    error = function(e) max(0, h2_hat - 1.96 * h2_se)
  )
  ci_hi  <- tryCatch(
    uniroot(function(h) prof(h) - ll_thr, c(h2_hat, 0.999), tol = 1e-6)$root,
    error = function(e) min(1, h2_hat + 1.96 * h2_se)
  )

  if (verbose) {
    cli::cli_alert_success(
      "{tag}  n={n}  h2={round(h2_hat,3)} [{round(ci_lo,3)},{round(ci_hi,3)}]  p={signif(pval,3)}"
    )
  }

  list(
    label          = tag,
    trait          = trait,
    covariates     = paste(covs, collapse = "+"),
    n              = n,
    h2             = round(h2_hat, 4),
    se             = round(h2_se,  4),
    ci_lo          = round(ci_lo,  4),
    ci_hi          = round(ci_hi,  4),
    pval           = round(pval,   5),
    var_covariates = var_covariates,
    sigma2_a       = round(sa,     5),
    sigma2_e       = round(se2,    5)
  )
}


# Null-coalescing operator -- available internally only
`%||%` <- function(x, y) if (is.null(x)) y else x
