test_that("build_grm returns a square matrix with correct dimnames", {
  ped <- make_test_pedigree()
  A   <- build_grm(ped, study_ids = 11:15)

  expect_true(is.matrix(A))
  expect_equal(nrow(A), 5L)
  expect_equal(ncol(A), 5L)
  expect_equal(rownames(A), as.character(11:15))
  expect_equal(colnames(A), as.character(11:15))
})

test_that("build_grm diagonal is 1 for non-inbred individuals", {
  ped <- make_test_pedigree()
  A   <- build_grm(ped, study_ids = 11:15)
  # unname() because diag() on a named matrix returns a named vector
  expect_equal(unname(diag(A)), rep(1, 5), tolerance = 1e-6)
})

test_that("build_grm sibling pairs have off-diagonal ~ 0.5", {
  ped <- make_test_pedigree()
  A   <- build_grm(ped, study_ids = 11:15)
  # IDs 11 & 12 are full siblings (parents 7 x 10); expected kinship = 0.5
  expect_equal(A["11", "12"], 0.5, tolerance = 1e-6)
  # IDs 13 & 14 are full siblings (parents 8 x 9); expected kinship = 0.5
  expect_equal(A["13", "14"], 0.5, tolerance = 1e-6)
})

test_that("build_grm unrelated pairs have off-diagonal 0", {
  ped <- make_test_pedigree()
  A   <- build_grm(ped, study_ids = 11:15)
  # IDs 11 & 13 have no shared ancestry in this pedigree
  expect_equal(A["11", "13"], 0, tolerance = 1e-6)
})

test_that("build_grm off-diagonal >= 0", {
  ped <- make_test_pedigree()
  A   <- build_grm(ped, study_ids = 11:15)
  expect_true(all(A[lower.tri(A)] >= 0))
})

test_that("build_grm is symmetric", {
  ped <- make_test_pedigree()
  A   <- build_grm(ped, study_ids = 11:15)
  expect_equal(A, t(A))
})

test_that("build_grm errors on missing pedigree columns", {
  ped <- make_test_pedigree()
  expect_error(build_grm(ped[, -1]),  # drop 'id' column
               regexp = "Required pedigree columns not found")
})

test_that("build_grm errors on study_ids not in pedigree", {
  ped <- make_test_pedigree()
  expect_error(build_grm(ped, study_ids = c(11, 99)),
               regexp = "not found in the pedigree")
})

test_that("build_grm warns on unrecognised sex values", {
  ped        <- make_test_pedigree()
  ped$sex[1] <- 9L
  expect_warning(build_grm(ped, study_ids = 11:15),
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
  A <- build_grm(ped, study_ids = 11:15,
                 id_col = "ID", pat_col = "DAD",
                 mom_col = "MOM", sex_col = "SEX")
  expect_equal(dim(A), c(5L, 5L))
})
