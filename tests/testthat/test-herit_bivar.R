helper_bivar_setup <- function() {
  d <- make_twin_data(n_mz = 20, n_dz = 20)
  A <- build_grm(d$full_ped, study_ids = d$study_ids,
                id_col = "id", pat_col = "pat", mom_col = "mom",
                sex_col = "sex", mz_col = "mztwin")
  dat <- d$data
  dat$IID <- dat$id
  list(A = A, dat = dat)
}

test_that("herit_bivar returns a list with expected elements", {
  s   <- helper_bivar_setup()
  res <- herit_bivar("trait1", "trait2", s$A, s$dat, id_col = "IID",
                     min_n = 30, verbose = FALSE)
  expect_type(res, "list")
  expect_named(res, c("trait1", "trait2", "n", "h2_1", "h2_2",
                      "rhoG", "rhoE", "rhoP", "p_rhoG", "p_rhoE", "p_rhoP"))
})

test_that("herit_bivar correlations are in [-1, 1]", {
  s   <- helper_bivar_setup()
  res <- herit_bivar("trait1", "trait2", s$A, s$dat, id_col = "IID",
                     min_n = 30, verbose = FALSE)
  for (r in c(res$rhoG, res$rhoE, res$rhoP)) {
    expect_gte(r, -1); expect_lte(r, 1)
  }
})

test_that("herit_bivar heritabilities are in [0, 1]", {
  s   <- helper_bivar_setup()
  res <- herit_bivar("trait1", "trait2", s$A, s$dat, id_col = "IID",
                     min_n = 30, verbose = FALSE)
  expect_gte(res$h2_1, 0); expect_lte(res$h2_1, 1)
  expect_gte(res$h2_2, 0); expect_lte(res$h2_2, 1)
})

test_that("herit_bivar p-values are in [0, 1]", {
  s   <- helper_bivar_setup()
  res <- herit_bivar("trait1", "trait2", s$A, s$dat, id_col = "IID",
                     min_n = 30, verbose = FALSE)
  for (p in c(res$p_rhoG, res$p_rhoE, res$p_rhoP)) {
    expect_gte(p, 0); expect_lte(p, 1)
  }
})

test_that("herit_bivar rhoP matches its derivation from rhoG, rhoE, h2_1, h2_2", {
  s   <- helper_bivar_setup()
  res <- herit_bivar("trait1", "trait2", s$A, s$dat, id_col = "IID",
                     min_n = 30, verbose = FALSE)
  expected_rhoP <- res$rhoG * sqrt(res$h2_1 * res$h2_2) +
    res$rhoE * sqrt((1 - res$h2_1) * (1 - res$h2_2))
  expect_equal(res$rhoP, round(expected_rhoP, 4), tolerance = 1e-3)
})

test_that("herit_bivar detects the strong simulated genetic/environmental correlation", {
  # trait2 is trait1 plus small noise, so both rhoG and rhoE should be large
  # and highly significant.
  s   <- helper_bivar_setup()
  res <- herit_bivar("trait1", "trait2", s$A, s$dat, id_col = "IID",
                     min_n = 30, verbose = FALSE)
  expect_gt(abs(res$rhoG), 0.5)
  expect_gt(abs(res$rhoE), 0.5)
  expect_lt(res$p_rhoE, 0.001)
})

test_that("herit_bivar returns NULL when n < min_n", {
  s <- helper_bivar_setup()
  expect_null(
    herit_bivar("trait1", "trait2", s$A, s$dat, id_col = "IID",
               min_n = 1000, verbose = FALSE)
  )
})

test_that("herit_bivar errors on missing trait column", {
  s <- helper_bivar_setup()
  expect_error(
    herit_bivar("nonexistent", "trait2", s$A, s$dat, id_col = "IID", verbose = FALSE),
    regexp = "not found"
  )
})

test_that("herit_bivar_batch returns a data frame with BH-adjusted q-values", {
  s     <- helper_bivar_setup()
  pairs <- data.frame(trait1 = "trait1", trait2 = "trait2")
  res   <- herit_bivar_batch(pairs, s$A, s$dat, id_col = "IID",
                             min_n = 30, .progress = FALSE)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1L)
  expect_true(all(c("q_rhoG", "q_rhoE", "q_rhoP") %in% names(res)))
})

test_that("herit_bivar_batch q-values are >= raw p-values (BH correction)", {
  s     <- helper_bivar_setup()
  pairs <- rbind(
    data.frame(trait1 = "trait1", trait2 = "trait2"),
    data.frame(trait1 = "trait2", trait2 = "trait1")
  )
  res <- herit_bivar_batch(pairs, s$A, s$dat, id_col = "IID", min_n = 30, .progress = FALSE)
  expect_true(all(res$q_rhoG >= res$p_rhoG - 1e-8))
  expect_true(all(res$q_rhoE >= res$p_rhoE - 1e-8))
})
