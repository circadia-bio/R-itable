test_that("int_transform returns vector of same length", {
  x <- c(1, 5, 3, 2, 4)
  expect_length(int_transform(x), 5L)
})

test_that("int_transform output is approximately standard normal", {
  set.seed(1)
  x   <- rexp(500)
  out <- int_transform(x)
  expect_lt(abs(mean(out, na.rm = TRUE)), 0.05)
  expect_lt(abs(sd(out, na.rm = TRUE) - 1), 0.05)
})

test_that("int_transform preserves NA positions", {
  x <- c(1, NA, 3, NA, 5)
  out <- int_transform(x)
  expect_true(is.na(out[2]))
  expect_true(is.na(out[4]))
  expect_false(is.na(out[1]))
})

test_that("int_transform is monotone (rank-preserving)", {
  set.seed(2)
  x   <- rnorm(50)
  out <- int_transform(x)
  expect_equal(rank(x), rank(out))
})
