test_that("build_grm returns a square matrix with correct dimnames", {
  ped <- make_test_pedigree()
  A   <- build_grm(ped, study_ids = 13:16)

  expect_true(is.matrix(A))
  expect_equal(nrow(A), 4L)
  expect_equal(ncol(A), 4L)
  expect_equal(rownames(A), as.character(13:16))
  expect_equal(colnames(A), as.character(13:16))
})

test_that("build_grm diagonal is 1 for non-inbred individuals", {
  ped <- make_test_pedigree()
  A   <- build_grm(ped, study_ids = 13:16)
  expect_equal(unname(diag(A)), rep(1, 4), tolerance = 1e-6)
})

test_that("build_grm full sibling pairs have off-diagonal 0.5", {
  ped <- make_test_pedigree()
  A   <- build_grm(ped, study_ids = 13:16)
  # 13 and 14 are full siblings (parents 9 x 10)
  expect_equal(A["13", "14"], 0.5, tolerance = 1e-6)
  # 15 and 16 are full siblings (parents 11 x 12)
  expect_equal(A["15", "16"], 0.5, tolerance = 1e-6)
})

test_that("build_grm truly unrelated pairs have off-diagonal 0", {
  ped <- make_test_pedigree()
  A   <- build_grm(ped, study_ids = 13:16)
  # 13 (family A x B) and 15 (family C x D) share no founders
  expect_equal(A["13", "15"], 0, tolerance = 1e-6)
  expect_equal(A["13", "16"], 0, tolerance = 1e-6)
  expect_equal(A["14", "15"], 0, tolerance = 1e-6)
  expect_equal(A["14", "16"], 0, tolerance = 1e-6)
})

test_that("build_grm off-diagonal >= 0", {
  ped <- make_test_pedigree()
  A   <- build_grm(ped, study_ids = 13:16)
  expect_true(all(A[lower.tri(A)] >= 0))
})

test_that("build_grm is symmetric", {
  ped <- make_test_pedigree()
  A   <- build_grm(ped, study_ids = 13:16)
  expect_equal(A, t(A))
})

test_that("build_grm errors on missing pedigree columns", {
  ped <- make_test_pedigree()
  expect_error(build_grm(ped[, -1]),
               regexp = "Required pedigree columns not found")
})

test_that("build_grm errors on study_ids not in pedigree", {
  ped <- make_test_pedigree()
  expect_error(build_grm(ped, study_ids = c(13, 99)),
               regexp = "not found in the pedigree")
})

test_that("build_grm warns on unrecognised sex values", {
  ped        <- make_test_pedigree()
  ped$sex[1] <- 9L
  expect_warning(build_grm(ped, study_ids = 13:16),
                 regexp = "unrecognised sex")
})

test_that("build_grm works with default study_ids (all pedigree members)", {
  ped <- make_test_pedigree()
  A   <- build_grm(ped)
  expect_equal(nrow(A), nrow(ped))
})

test_that("build_grm accepts custom column names", {
  ped        <- make_test_pedigree()
  names(ped) <- c("ID", "DAD", "MOM", "SEX")
  A <- build_grm(ped, study_ids = 13:16,
                 id_col = "ID", pat_col = "DAD",
                 mom_col = "MOM", sex_col = "SEX")
  expect_equal(dim(A), c(4L, 4L))
})
