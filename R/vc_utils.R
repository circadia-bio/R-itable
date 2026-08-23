# Internal helpers for general (non-eigendecomposition) ML variance-component
# models: household (C) component and bivariate genetic/environmental
# correlation. Unlike herit_vc()'s 1-D profile likelihood (which relies on A
# alone being diagonalisable), these models involve two or more variance
# components / covariance blocks that do not jointly diagonalise, so they use
# direct numerical ML over the full n x n (or 2n x 2n) covariance matrix via
# Cholesky decomposition. Not exported.

# -- Univariate GLS log-likelihood for Omega = h2*A + c2*C + (1-h2-c2)*I ------
.loglik_uni_ace <- function(h2, c2, y, X, A, C) {
  n  <- length(y)
  Om <- h2 * A + c2 * C + (1 - h2 - c2) * diag(n)
  ch <- tryCatch(chol(Om), error = function(e) NULL)
  if (is.null(ch)) return(-1e10)
  Oinv   <- chol2inv(ch)
  logdet <- 2 * sum(log(diag(ch)))
  XtOiX  <- t(X) %*% Oinv %*% X
  beta   <- tryCatch(solve(XtOiX, t(X) %*% Oinv %*% y), error = function(e) NULL)
  if (is.null(beta)) return(-1e10)
  resid <- y - X %*% beta
  sp    <- as.numeric(t(resid) %*% Oinv %*% resid) / n
  if (!is.finite(sp) || sp <= 0) return(-1e10)
  -0.5 * (n * log(sp) + logdet + n)
}

# -- Fit one of: none (E only), A (AE), C (CE), AC (full ACE) -----------------
.fit_ace_uni <- function(y, X, A, C, estimate = c("none", "A", "C", "AC")) {
  estimate <- match.arg(estimate)

  if (estimate == "none") {
    return(list(h2 = 0, c2 = 0, loglik = .loglik_uni_ace(0, 0, y, X, A, C)))
  }
  if (estimate == "A") {
    op <- stats::optimize(function(p) -.loglik_uni_ace(stats::plogis(p), 0, y, X, A, C),
                          interval = c(-12, 12))
    return(list(h2 = stats::plogis(op$minimum), c2 = 0, loglik = -op$objective))
  }
  if (estimate == "C") {
    op <- stats::optimize(function(p) -.loglik_uni_ace(0, stats::plogis(p), y, X, A, C),
                          interval = c(-12, 12))
    return(list(h2 = 0, c2 = stats::plogis(op$minimum), loglik = -op$objective))
  }

  # estimate == "AC": h2, c2 >= 0, h2 + c2 < 1, via multi-start Nelder-Mead
  obj <- function(par) {
    h2 <- par[1]; c2 <- par[2]
    if (!is.finite(h2) || !is.finite(c2) || h2 < 0 || c2 < 0 || h2 + c2 > 0.999) return(1e10)
    -.loglik_uni_ace(h2, c2, y, X, A, C)
  }
  starts <- list(c(0.01, 0.01), c(0.3, 0.3), c(0.5, 0.1), c(0.1, 0.5), c(0.4, 0), c(0, 0.3))
  best <- NULL
  for (st in starts) {
    op <- tryCatch(
      stats::optim(st, obj, method = "Nelder-Mead",
                  control = list(reltol = 1e-10, maxit = 2000)),
      error = function(e) NULL
    )
    if (!is.null(op) && (is.null(best) || op$value < best$value)) best <- op
  }
  list(h2 = best$par[1], c2 = best$par[2], loglik = -best$value)
}

# -- Boundary-corrected one-sided LRT (H0 on parameter-space boundary) --------
.lrt_boundary <- function(ll_full, ll_null, df = 1) {
  lrt <- max(0, 2 * (ll_full - ll_null))
  0.5 * stats::pchisq(lrt, df = df, lower.tail = FALSE)
}

# -- Bivariate AE model, eigenbasis reduction ---------------------------------
# Both covariance blocks (within-trait and cross-trait) are linear
# combinations of A and I, so they share A's eigenvectors. Diagonalising A
# ONCE up front reduces every likelihood evaluation from an O(n^3) 2n x 2n
# Cholesky to n independent, cheap 2x2 blocks -- the same trick herit_vc()
# uses for the univariate case, extended to two traits. This is what makes
# herit_bivar_batch() over ~30 pairs tractable.
#
# sigma2_p1/sigma2_p2 ARE free ML-searched parameters here (unlike an
# intermediate version of this file, which fixed them to the sample
# variance to mirror SOLAR's "sd initialized in c++" architecture). That
# turned out to be a regression: SOLAR's sd/mean appear to be given a good
# starting value from the data but still refined by the outer optimiser,
# not permanently fixed -- free search matches SOLAR's actual numeric
# output more closely on well-conditioned traits. Numerical robustness for
# extreme-raw-scale traits (e.g. a corrupted column with variance ~1e18) is
# instead handled by standardising both traits to unit variance BEFORE
# calling into this code (see herit_bivar()) -- h2/rho are scale-invariant,
# so this changes nothing about the true ML answer, but keeps every
# covariance-matrix entry in a numerically well-conditioned O(1) range
# regardless of the traits' native units.
#
# rhoG/rhoE are bounded to [-0.9, 0.9], matching SOLAR's actual default
# parameter bounds (confirmed from source, lib/solar.tcl proc polymod:
# `set rholow -0.9; set rhoup 0.9`) rather than the mathematical [-1, 1].

# Precompute the pieces of the bivariate problem that don't depend on the
# variance-component parameters: A's eigendecomposition and the traits
# rotated into that eigenbasis.
.bivar_eigen_prep <- function(y1, y2, A) {
  eig <- eigen(A, symmetric = TRUE)
  Vt  <- t(eig$vectors)
  list(ev = eig$values, y1r = as.numeric(Vt %*% y1), y2r = as.numeric(Vt %*% y2),
      x1r = as.numeric(Vt %*% rep(1, length(y1))),
      x2r = as.numeric(Vt %*% rep(1, length(y2))))
}

# par = (logit h2_1, logit h2_2, atanh-scaled rhoG in [-0.9,0.9],
#        atanh-scaled rhoE in [-0.9,0.9], log sigma2_p1, log sigma2_p2)
.nll_bivar_eigen <- function(par, prep) {
  h2_1 <- stats::plogis(par[1]); h2_2 <- stats::plogis(par[2])
  rG   <- 0.9 * tanh(par[3]);    rE   <- 0.9 * tanh(par[4])
  s1   <- exp(par[5]);           s2   <- exp(par[6])

  sg1 <- h2_1 * s1; sg2 <- h2_2 * s2
  se1 <- (1 - h2_1) * s1; se2 <- (1 - h2_2) * s2
  cg  <- rG * sqrt(sg1 * sg2); ce <- rE * sqrt(se1 * se2)

  ev <- prep$ev
  # Per-channel 2x2 covariance: [[sg1*ev+se1, cg*ev+ce], [cg*ev+ce, sg2*ev+se2]]
  a11 <- sg1 * ev + se1
  a22 <- sg2 * ev + se2
  a12 <- cg  * ev + ce
  det <- a11 * a22 - a12^2
  if (any(!is.finite(det)) || any(det <= 0)) return(1e10)
  i11 <- a22 / det; i22 <- a11 / det; i12 <- -a12 / det   # per-channel inverse entries

  # GLS for (mu1, mu2): accumulate the 2x2 information matrix and score
  # across channels (each channel contributes a rank-1-per-trait term).
  x1 <- prep$x1r; x2 <- prep$x2r; y1 <- prep$y1r; y2 <- prep$y2r
  I11 <- sum(i11 * x1^2);  I22 <- sum(i22 * x2^2);  I12 <- sum(i12 * x1 * x2)
  s1v <- sum(i11 * x1 * y1 + i12 * x1 * y2)
  s2v <- sum(i12 * x2 * y1 + i22 * x2 * y2)
  Imat <- matrix(c(I11, I12, I12, I22), 2, 2)
  mu <- tryCatch(solve(Imat, c(s1v, s2v)), error = function(e) NULL)
  if (is.null(mu)) return(1e10)

  r1 <- y1 - x1 * mu[1]; r2 <- y2 - x2 * mu[2]
  quad   <- sum(i11 * r1^2 + 2 * i12 * r1 * r2 + i22 * r2^2)
  logdet <- sum(log(det))
  n <- length(ev)
  ll <- -0.5 * (logdet + quad + 2 * n * log(2 * pi))
  -ll
}

.fit_bivar_full <- function(y1, y2, A) {
  prep <- .bivar_eigen_prep(y1, y2, A)
  # Data-driven scale starts: assuming unit variance (log(1)=0) as the ONLY
  # start fails badly for traits with a large or small raw variance; a
  # data-driven start plus SOLAR's own conventional start (h2r=0.1) both
  # guard against that, and against local optima on weakly-identified pairs.
  s1_0 <- log(max(var(y1), 1e-8))
  s2_0 <- log(max(var(y2), 1e-8))
  h2_0 <- qlogis(0.1)
  starts <- list(
    c(h2_0, h2_0, 0, 0, s1_0, s2_0),
    c(0, 0, 0, 0, s1_0, s2_0),
    c(0.5, 0.5, 0.5, 0.5, s1_0, s2_0),
    c(-0.5, -0.5, 0.2, 0.2, s1_0, s2_0),
    c(0, 0, 0, 0, 0, 0)  # unit-variance fallback
  )
  best <- NULL
  for (st in starts) {
    op <- tryCatch(
      stats::optim(st, .nll_bivar_eigen, prep = prep, method = "Nelder-Mead",
                  control = list(maxit = 3000, reltol = 1e-12)),
      error = function(e) NULL
    )
    if (!is.null(op) && (is.null(best) || op$value < best$value)) best <- op
  }
  # Polish the best candidate with a second Nelder-Mead pass from its own
  # optimum -- guards against the flat/weakly-identified likelihood surfaces
  # seen for low-signal trait pairs, where a single pass can stop early.
  op2 <- tryCatch(
    stats::optim(best$par, .nll_bivar_eigen, prep = prep, method = "Nelder-Mead",
                control = list(maxit = 3000, reltol = 1e-13)),
    error = function(e) NULL
  )
  if (!is.null(op2) && op2$value < best$value) best <- op2

  h2_1 <- stats::plogis(best$par[1]); h2_2 <- stats::plogis(best$par[2])
  rG   <- 0.9 * tanh(best$par[3]);    rE   <- 0.9 * tanh(best$par[4])
  list(h2_1 = h2_1, h2_2 = h2_2, rhoG = rG, rhoE = rE,
      loglik = -best$value, par = best$par, prep = prep)
}

# Fit constrained models where rhoG and/or rhoE are pinned to 0. Reuses the
# eigen-prep (and, as a warm start, the full model's fitted parameters) from
# .fit_bivar_full() so no repeated eigendecomposition is needed.
.fit_bivar_constrained <- function(y1, y2, A, fix_rhoG = FALSE, fix_rhoE = FALSE,
                                   start, prep = NULL) {
  if (is.null(prep)) prep <- .bivar_eigen_prep(y1, y2, A)
  free_idx <- setdiff(seq_len(6), c(if (fix_rhoG) 3, if (fix_rhoE) 4))
  full_par <- function(p) {
    par <- numeric(6)
    par[free_idx] <- p
    if (fix_rhoG) par[3] <- 0
    if (fix_rhoE) par[4] <- 0
    par
  }
  op <- stats::optim(start[free_idx], function(p) .nll_bivar_eigen(full_par(p), prep),
                     method = "Nelder-Mead", control = list(maxit = 2000, reltol = 1e-10))
  -op$value
}

# -- Unbalanced bivariate model (SOLAR's default UnbalancedTraits = 1) --------
# When one trait has missing values, SOLAR's default behaviour is to keep
# individuals missing exactly one of the two traits in the analysis sample,
# using their observed trait's marginal contribution to the likelihood,
# rather than dropping them via a complete-case restriction. This matters
# for any trait pair where one phenotype has missing data: dropping those
# individuals both wastes information and, empirically, shifts rhoG/rhoE
# estimates noticeably (verified against SOLAR ground truth). Unlike the
# balanced case, the joint covariance's missingness pattern is individual-
# specific, so this does not reduce to n independent 2x2 eigen-blocks; it
# uses a general Cholesky-based GLS over the variable-size joint covariance
# of just the observed trait-slots (fast enough at this sample size; only
# triggered when there IS missingness in one of the two traits).
.nll_bivar_unbalanced <- function(par, y1, y2, A, has_y1, has_y2) {
  h2_1 <- stats::plogis(par[1]); h2_2 <- stats::plogis(par[2])
  rG   <- 0.9 * tanh(par[3]);    rE   <- 0.9 * tanh(par[4])
  s1   <- exp(par[5]);           s2   <- exp(par[6])
  n <- nrow(A)
  sg1 <- h2_1 * s1; sg2 <- h2_2 * s2
  se1 <- (1 - h2_1) * s1; se2 <- (1 - h2_2) * s2
  I <- diag(n)
  S11 <- sg1 * A + se1 * I
  S22 <- sg2 * A + se2 * I
  S12 <- rG * sqrt(sg1 * sg2) * A + rE * sqrt(se1 * se2) * I

  idx1 <- which(has_y1); idx2 <- which(has_y2)
  Sigma <- rbind(cbind(S11[idx1, idx1, drop = FALSE], S12[idx1, idx2, drop = FALSE]),
                cbind(t(S12[idx1, idx2, drop = FALSE]), S22[idx2, idx2, drop = FALSE]))
  ch <- tryCatch(chol(Sigma), error = function(e) NULL)
  if (is.null(ch)) return(1e10)
  logdet <- 2 * sum(log(diag(ch)))
  Sinv   <- chol2inv(ch)
  yv <- c(y1[idx1], y2[idx2])
  X  <- cbind(c(rep(1, length(idx1)), rep(0, length(idx2))),
             c(rep(0, length(idx1)), rep(1, length(idx2))))
  XtSiX <- t(X) %*% Sinv %*% X
  beta  <- tryCatch(solve(XtSiX, t(X) %*% Sinv %*% yv), error = function(e) NULL)
  if (is.null(beta)) return(1e10)
  resid <- yv - X %*% beta
  ll <- -0.5 * (logdet + as.numeric(t(resid) %*% Sinv %*% resid) + length(yv) * log(2 * pi))
  -ll
}

.fit_bivar_unbalanced_full <- function(y1, y2, A, has_y1, has_y2) {
  y1z <- y1; y1z[!has_y1] <- 0
  y2z <- y2; y2z[!has_y2] <- 0
  s1_0 <- log(max(stats::var(y1[has_y1]), 1e-8))
  s2_0 <- log(max(stats::var(y2[has_y2]), 1e-8))
  h2_0 <- qlogis(0.1)
  starts <- list(
    c(h2_0, h2_0, 0, 0, s1_0, s2_0),
    c(-0.5, -0.5, 0.2, 0.2, s1_0, s2_0),
    c(0, 0, 0, 0, 0, 0)
  )
  best <- NULL
  for (st in starts) {
    op <- tryCatch(
      stats::optim(st, .nll_bivar_unbalanced, y1 = y1z, y2 = y2z, A = A,
                  has_y1 = has_y1, has_y2 = has_y2, method = "Nelder-Mead",
                  control = list(maxit = 3000, reltol = 1e-11)),
      error = function(e) NULL
    )
    if (!is.null(op) && (is.null(best) || op$value < best$value)) best <- op
  }

  h2_1 <- stats::plogis(best$par[1]); h2_2 <- stats::plogis(best$par[2])
  rG   <- 0.9 * tanh(best$par[3]);    rE   <- 0.9 * tanh(best$par[4])
  list(h2_1 = h2_1, h2_2 = h2_2, rhoG = rG, rhoE = rE,
      loglik = -best$value, par = best$par)
}

.fit_bivar_unbalanced_constrained <- function(y1, y2, A, has_y1, has_y2,
                                              fix_rhoG = FALSE, fix_rhoE = FALSE, start) {
  y1z <- y1; y1z[!has_y1] <- 0
  y2z <- y2; y2z[!has_y2] <- 0
  free_idx <- setdiff(seq_len(6), c(if (fix_rhoG) 3, if (fix_rhoE) 4))
  full_par <- function(p) {
    par <- numeric(6)
    par[free_idx] <- p
    if (fix_rhoG) par[3] <- 0
    if (fix_rhoE) par[4] <- 0
    par
  }
  op <- stats::optim(start[free_idx],
                     function(p) .nll_bivar_unbalanced(full_par(p), y1z, y2z, A, has_y1, has_y2),
                     method = "Nelder-Mead", control = list(maxit = 2000, reltol = 1e-10))
  -op$value
}
