helper_ace_setup <- function() {
  d <- make_twin_data(n_mz = 20, n_dz = 20)
  A <- build_grm(d$full_ped, study_ids = d$study_ids,
                id_col = "id", pat_col = "pat", mom_col = "mom",
                sex_col = "sex", mz_col = "mztwin")
  C <- build_household(d$ped_df$hhid, d$ped_df$id)
  dat <- d$data
  dat$IID <- dat$id
  list(A = A, C = C, dat = dat)
}

test_that("herit_ace returns a list with expected elements", {
  s   <- helper_ace_setup()
  res <- herit_ace("trait1", s$A, s$C, s$dat, id_col = "IID",
                   min_n = 30, verbose = FALSE)
  expect_type(res, "list")
  expect_named(res, c("trait", "n", "loglik_sporadic", "loglik_ae", "loglik_ce",
                      "loglik_ace", "h2_ae", "c2_ce", "h2_ace", "c2_ace",
                      "chisq_c_vs_sporadic", "p_c_vs_sporadic",
                      "chisq_ace_vs_ae", "p_ace_vs_ae"))
})

test_that("herit_ace h2 and c2 estimates are in [0, 1]", {
  s   <- helper_ace_setup()
  res <- herit_ace("trait1", s$A, s$C, s$dat, id_col = "IID",
                   min_n = 30, verbose = FALSE)
  expect_gte(res$h2_ae, 0); expect_lte(res$h2_ae, 1)
  expect_gte(res$c2_ce, 0); expect_lte(res$c2_ce, 1)
  expect_gte(res$h2_ace, 0); expect_lte(res$h2_ace, 1)
  expect_gte(res$c2_ace, 0); expect_lte(res$c2_ace, 1)
  expect_lte(res$h2_ace + res$c2_ace, 1)
})

test_that("herit_ace p-values are in [0, 1]", {
  s   <- helper_ace_setup()
  res <- herit_ace("trait1", s$A, s$C, s$dat, id_col = "IID",
                   min_n = 30, verbose = FALSE)
  expect_gte(res$p_c_vs_sporadic, 0); expect_lte(res$p_c_vs_sporadic, 1)
  expect_gte(res$p_ace_vs_ae, 0); expect_lte(res$p_ace_vs_ae, 1)
})

test_that("herit_ace nested log-likelihoods are ordered as expected", {
  s   <- helper_ace_setup()
  res <- herit_ace("trait1", s$A, s$C, s$dat, id_col = "IID",
                   min_n = 30, verbose = FALSE)
  expect_gte(res$loglik_ace, res$loglik_ae - 1e-6)
  expect_gte(res$loglik_ace, res$loglik_ce - 1e-6)
  expect_gte(res$loglik_ae,  res$loglik_sporadic - 1e-6)
  expect_gte(res$loglik_ce,  res$loglik_sporadic - 1e-6)
})

test_that("herit_ace simulated household-only signal is detected", {
  s   <- helper_ace_setup()
  res <- herit_ace("trait1", s$A, s$C, s$dat, id_col = "IID",
                   min_n = 30, verbose = FALSE)
  expect_lt(res$p_c_vs_sporadic, 0.05)
})

test_that("herit_ace returns NULL when n < min_n", {
  s <- helper_ace_setup()
  expect_null(
    herit_ace("trait1", s$A, s$C, s$dat, id_col = "IID",
             min_n = 1000, verbose = FALSE)
  )
})

test_that("herit_ace errors on missing trait column", {
  s <- helper_ace_setup()
  expect_error(
    herit_ace("nonexistent", s$A, s$C, s$dat, id_col = "IID", verbose = FALSE),
    regexp = "not found"
  )
})

test_that("herit_ace errors when grm and household dimnames don't match", {
  s <- helper_ace_setup()
  C_bad <- s$C
  rownames(C_bad) <- colnames(C_bad) <- rev(rownames(C_bad))
  expect_error(
    herit_ace("trait1", s$A, C_bad, s$dat, id_col = "IID", verbose = FALSE),
    regexp = "identical row/column names"
  )
})

# -- Covariate screening (SOLAR -screen semantics) ---------------------------

test_that("herit_ace without covs matches prior behaviour exactly (backward compatible)", {
  s   <- helper_ace_setup()
  res <- herit_ace("trait1", s$A, s$C, s$dat, id_col = "IID", min_n = 30, verbose = FALSE)
  expect_named(res, c("trait", "n", "loglik_sporadic", "loglik_ae", "loglik_ce",
                      "loglik_ace", "h2_ae", "c2_ce", "h2_ace", "c2_ace",
                      "chisq_c_vs_sporadic", "p_c_vs_sporadic",
                      "chisq_ace_vs_ae", "p_ace_vs_ae"))
})

test_that("herit_ace with covs adds screening fields to the output", {
  s <- helper_ace_setup()
  set.seed(1)
  s$dat$cov1 <- rnorm(nrow(s$dat))
  res <- herit_ace("trait1", s$A, s$C, s$dat, covs = "cov1", id_col = "IID",
                   min_n = 30, verbose = FALSE)
  expect_true(all(c("covs_tested", "covs_kept", "covs_dropped", "covariate_lrt_p") %in% names(res)))
  expect_equal(res$covs_tested, "cov1")
})

test_that("herit_ace screening keeps a strongly predictive covariate over a p-value threshold", {
  # A covariate near-perfectly correlated with the trait should always have
  # an astronomically small LRT p-value, regardless of any other randomness
  # in the simulated data -- a robust (non-flaky) way to test the keep path.
  s <- helper_ace_setup()
  set.seed(2)
  s$dat$strong_cov <- s$dat$trait1 * 5 + rnorm(nrow(s$dat), sd = 0.01)
  res <- herit_ace("trait1", s$A, s$C, s$dat, covs = "strong_cov", id_col = "IID",
                   min_n = 30, verbose = FALSE)
  expect_lt(res$covariate_lrt_p[["strong_cov"]], 1e-10)
  expect_true("strong_cov" %in% res$covs_kept)
})

test_that("herit_ace screen = FALSE keeps all covariates unconditionally", {
  s <- helper_ace_setup()
  set.seed(3)
  s$dat$noise_cov <- rnorm(nrow(s$dat))
  res <- herit_ace("trait1", s$A, s$C, s$dat, covs = "noise_cov", screen = FALSE,
                   id_col = "IID", min_n = 30, verbose = FALSE)
  expect_false("covs_kept" %in% names(res))
})

test_that("herit_ace drops zero-variance covariates with a warning", {
  s <- helper_ace_setup()
  s$dat$const_cov <- 5
  expect_warning(
    herit_ace("trait1", s$A, s$C, s$dat, covs = "const_cov", id_col = "IID",
             min_n = 30, verbose = FALSE),
    regexp = "zero-variance"
  )
})

test_that("herit_ace covariate LRT p-values are in [0, 1]", {
  s <- helper_ace_setup()
  set.seed(4)
  s$dat$cov1 <- rnorm(nrow(s$dat))
  s$dat$cov2 <- s$dat$trait1 + rnorm(nrow(s$dat))
  res <- herit_ace("trait1", s$A, s$C, s$dat, covs = c("cov1", "cov2"), id_col = "IID",
                   min_n = 30, verbose = FALSE)
  expect_true(all(res$covariate_lrt_p >= 0 & res$covariate_lrt_p <= 1))
})
