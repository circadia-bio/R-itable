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

test_that("int_transform matches Van der Waerden exactly when there are no ties", {
  set.seed(3)
  x   <- rnorm(40)
  n   <- length(x)
  out <- int_transform(x)
  expected <- qnorm(rank(x) / (n + 1))
  expect_equal(out, expected, tolerance = 1e-10)
})

test_that("int_transform averages z-scores (not ranks) within tied blocks", {
  # SOLAR's inormal averages qnorm(i/(n+1)) over each tied block's individual
  # rank positions, which differs from qnorm(mean(rank)/(n+1)) whenever the
  # tied block isn't symmetric around the middle of the sample -- e.g. a
  # large block of ties at one tail (zero-inflated data).
  x <- c(rep(0, 40), round(runif(60, 1, 100)))
  n <- length(x)
  out <- int_transform(x)

  tied_val <- out[x == 0][1]
  expect_true(all(out[x == 0] == tied_val))  # all tied individuals get the same value

  naive_rank_then_z <- qnorm(mean(rank(x)[x == 0]) / (n + 1))
  expect_false(isTRUE(all.equal(tied_val, naive_rank_then_z)))

  # The 40 zeros occupy sequential rank positions 1..40 (lowest values);
  # SOLAR-exact = mean of qnorm(i/(n+1)) over those individual positions.
  correct_z_then_average <- mean(qnorm((1:40) / (n + 1)))
  expect_equal(tied_val, correct_z_then_average, tolerance = 1e-10)
})

test_that("int_transform blom method matches the Blom formula with rank-average ties", {
  set.seed(4)
  x <- c(rep(2, 10), round(runif(30, 0, 20)))
  n <- length(x)
  out <- int_transform(x, method = "blom")
  expected <- qnorm((rank(x, ties.method = "average") - 0.5) / n)
  expect_equal(out, expected, tolerance = 1e-10)
})
