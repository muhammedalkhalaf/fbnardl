test_that("fbnardl FNARDL returns expected structure", {
  set.seed(1)
  n <- 80
  x <- cumsum(rnorm(n))
  y <- 0.6 * x + rnorm(n, sd = 0.5)
  res <- fbnardl(y, x, type = "fnardl", max_lag = 2, max_k = 2)
  expect_s3_class(res, "fbnardl")
  expect_equal(res$type, "fnardl")
  expect_true(is.numeric(res$coefficients))
  expect_true(res$nobs > 0)
})

test_that("fbnardl bootstrap runs with small reps", {
  set.seed(2)
  n <- 60
  x <- cumsum(rnorm(n))
  y <- 0.5 * x + rnorm(n)
  res <- fbnardl(y, x, type = "fbnardl", max_lag = 1, max_k = 1, reps = 20L)
  expect_s3_class(res, "fbnardl")
  expect_true(!is.null(res$bootstrap))
  expect_true(is.numeric(res$bootstrap$Fov_cv05))
})

test_that("fbnardl long-run multipliers are numeric", {
  set.seed(3)
  n <- 80
  x <- cumsum(rnorm(n))
  y <- 0.7 * x + rnorm(n, sd = 0.3)
  res <- fbnardl(y, x, type = "fnardl", max_lag = 2, max_k = 2)
  expect_true(is.numeric(res$lr_pos))
  expect_true(is.numeric(res$lr_neg))
})

test_that("fbnardl with control variable works", {
  set.seed(4)
  n <- 80
  x1 <- cumsum(rnorm(n))
  x2 <- rnorm(n)
  y  <- 0.5 * x1 + 0.2 * x2 + rnorm(n)
  res <- fbnardl(y, x1, x_ctrl = x2, type = "fnardl", max_lag = 1, max_k = 1)
  expect_s3_class(res, "fbnardl")
})

test_that("fbnardl no_fourier option works", {
  set.seed(5)
  n <- 70
  x <- cumsum(rnorm(n))
  y <- 0.5 * x + rnorm(n)
  res <- fbnardl(y, x, type = "fnardl", max_lag = 1, no_fourier = TRUE)
  expect_equal(res$best_kstar, 0)
})

test_that("print.fbnardl produces output", {
  set.seed(6)
  n <- 60
  x <- cumsum(rnorm(n))
  y <- 0.5 * x + rnorm(n)
  res <- fbnardl(y, x, type = "fnardl", max_lag = 1, max_k = 1)
  out <- capture.output(print(res))
  expect_true(length(out) > 0)
})
