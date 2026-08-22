test_that("build_household returns a matrix with correct dimnames", {
  hh  <- c("A", "A", "B", "B", "C")
  ids <- 1:5
  C <- build_household(hh, ids)
  expect_true(is.matrix(C))
  expect_equal(dim(C), c(5L, 5L))
  expect_equal(rownames(C), as.character(ids))
})

test_that("build_household diagonal is always 1", {
  hh  <- c("A", "A", "B", "B", "C")
  ids <- 1:5
  C <- build_household(hh, ids)
  expect_equal(unname(diag(C)), rep(1, 5))
})

test_that("build_household same-household pairs are 1", {
  hh  <- c("A", "A", "B", "B", "C")
  ids <- 1:5
  C <- build_household(hh, ids)
  expect_equal(C["1", "2"], 1)
  expect_equal(C["3", "4"], 1)
})

test_that("build_household different-household pairs are 0", {
  hh  <- c("A", "A", "B", "B", "C")
  ids <- 1:5
  C <- build_household(hh, ids)
  expect_equal(C["1", "3"], 0)
  expect_equal(C["1", "5"], 0)
})

test_that("build_household is symmetric", {
  hh  <- c("A", "A", "B", "B", "C")
  ids <- 1:5
  C <- build_household(hh, ids)
  expect_equal(C, t(C))
})

test_that("build_household errors on mismatched lengths", {
  expect_error(build_household(c("A", "B"), 1:3), regexp = "same length")
})
