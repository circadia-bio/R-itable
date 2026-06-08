test_that("herit_vc returns a list with expected elements", {
  d    <- make_family_data()
  A    <- build_grm(d$ped, study_ids = d$study_ids)
  res  <- herit_vc("trait1", grm = A, data = d$data, verbose = FALSE)

  expect_type(res, "list")
  expect_named(res, c("label","trait","covariates","n","h2","se",
                      "ci_lo","ci_hi","pval","var_covariates",
                      "sigma2_a","sigma2_e"))
})

test_that("herit_vc h2 is in [0, 1]", {
  d   <- make_family_data()
  A   <- build_grm(d$ped, study_ids = d$study_ids)
  res <- herit_vc("trait1", grm = A, data = d$data, verbose = FALSE)
  expect_gte(res$h2, 0)
  expect_lte(res$h2, 1)
})

test_that("herit_vc ci_lo <= h2 <= ci_hi", {
  d   <- make_family_data()
  A   <- build_grm(d$ped, study_ids = d$study_ids)
  res <- herit_vc("trait1", grm = A, data = d$data, verbose = FALSE)
  expect_lte(res$ci_lo, res$h2)
  expect_gte(res$ci_hi, res$h2)
})

test_that("herit_vc pval is in [0, 1]", {
  d   <- make_family_data()
  A   <- build_grm(d$ped, study_ids = d$study_ids)
  res <- herit_vc("trait1", grm = A, data = d$data, verbose = FALSE)
  expect_gte(res$pval, 0)
  expect_lte(res$pval, 1)
})

test_that("herit_vc returns NULL when n < min_n", {
  d    <- make_family_data(n_families = 5)   # only 10 offspring
  A    <- build_grm(d$ped, study_ids = d$study_ids)
  expect_null(
    suppressWarnings(herit_vc("trait1", grm = A, data = d$data,
                              min_n = 50, verbose = FALSE))
  )
})

test_that("herit_vc errors on missing trait column", {
  d <- make_family_data()
  A <- build_grm(d$ped, study_ids = d$study_ids)
  expect_error(herit_vc("nonexistent", grm = A, data = d$data, verbose = FALSE),
               regexp = "not found in")
})

test_that("herit_vc with covariates still returns valid h2", {
  d   <- make_family_data()
  A   <- build_grm(d$ped, study_ids = d$study_ids)
  res <- herit_vc("trait1", grm = A, data = d$data,
                  covs = c("age", "sex_num"), verbose = FALSE)
  expect_gte(res$h2, 0)
  expect_lte(res$h2, 1)
})

test_that("herit_vc var_covariates is NA for unadjusted, numeric for adjusted", {
  d   <- make_family_data()
  A   <- build_grm(d$ped, study_ids = d$study_ids)
  res_u <- herit_vc("trait1", grm = A, data = d$data, verbose = FALSE)
  res_a <- herit_vc("trait1", grm = A, data = d$data,
                    covs = c("age", "sex_num"), verbose = FALSE)
  expect_true(is.na(res_u$var_covariates))
  expect_gte(res_a$var_covariates, 0)
  expect_lte(res_a$var_covariates, 1)
})

test_that("herit_vc sigma2_a + sigma2_e > 0", {
  d   <- make_family_data()
  A   <- build_grm(d$ped, study_ids = d$study_ids)
  res <- herit_vc("trait1", grm = A, data = d$data, verbose = FALSE)
  expect_gt(res$sigma2_a + res$sigma2_e, 0)
})
