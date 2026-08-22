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
  # ACE nests AE and CE, both of which nest sporadic, so loglik_ace should be
  # the largest and loglik_sporadic the smallest.
  s   <- helper_ace_setup()
  res <- herit_ace("trait1", s$A, s$C, s$dat, id_col = "IID",
                   min_n = 30, verbose = FALSE)
  expect_gte(res$loglik_ace, res$loglik_ae - 1e-6)
  expect_gte(res$loglik_ace, res$loglik_ce - 1e-6)
  expect_gte(res$loglik_ae,  res$loglik_sporadic - 1e-6)
  expect_gte(res$loglik_ce,  res$loglik_sporadic - 1e-6)
})

test_that("herit_ace simulated household-only signal is detected", {
  # trait1 is simulated with a real shared-household component, so CE should
  # fit meaningfully better than the sporadic (no familial effect) model.
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
