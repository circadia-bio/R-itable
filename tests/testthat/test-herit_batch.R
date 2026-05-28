test_that("herit_batch returns a data frame", {
  d    <- make_family_data()
  A    <- build_grm(d$ped, study_ids = d$study_ids)
  res  <- herit_batch(c("trait1", "trait2"), grm = A, data = d$data,
                      .progress = FALSE)
  expect_s3_class(res, "data.frame")
})

test_that("herit_batch row count equals traits × models (no skips)", {
  d   <- make_family_data()
  A   <- build_grm(d$ped, study_ids = d$study_ids)
  res <- herit_batch(
    c("trait1", "trait2"), grm = A, data = d$data,
    covs_list  = list(unadj = NULL, adj = c("age", "sex_num")),
    .progress  = FALSE
  )
  expect_equal(nrow(res), 4L)  # 2 traits × 2 models
})

test_that("herit_batch label follows <trait>_<model> pattern", {
  d   <- make_family_data()
  A   <- build_grm(d$ped, study_ids = d$study_ids)
  res <- herit_batch("trait1", grm = A, data = d$data,
                     covs_list = list(unadj = NULL), .progress = FALSE)
  expect_equal(res$label, "trait1_unadj")
})

test_that("herit_batch all h2 values in [0, 1]", {
  d   <- make_family_data()
  A   <- build_grm(d$ped, study_ids = d$study_ids)
  res <- herit_batch(c("trait1","trait2"), grm = A, data = d$data,
                     .progress = FALSE)
  expect_true(all(res$h2 >= 0 & res$h2 <= 1))
})
