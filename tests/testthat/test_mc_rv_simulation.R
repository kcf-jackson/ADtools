context("Test random variable simulation")

test_that("Inverse transform is consistent with the base package implementation", {
  n <- 1000
  sample_1 <- rgamma0(n, 1, 1, method = "base")
  sample_2 <- rgamma0(n, 1, 1, method = "inv_tf")
  testthat::expect_lt(abs(mean(sample_1) - mean(sample_2)), 1 / sqrt(n) * 10)
  testthat::expect_error(rgamma0(n, 1, 1, method = "cause_trouble"))

  sample_1 <- rchisq0(n, 10, method = "base")
  sample_2 <- rchisq0(n, 10, method = "inv_tf")
  testthat::expect_lt(abs(mean(sample_1) - mean(sample_2)), 1 / sqrt(n) * 10)

  mean_list <- function(x) purrr::reduce(x, `+`) / length(x)
  sample_1 <- purrr::map(1:n, ~rWishart0(10, diag(3), method = "base"))
  sample_2 <- purrr::map(1:n, ~rWishart0(10, diag(3), method = "inv_tf"))
  testthat::expect_lt(max(abs(mean_list(sample_1) - mean_list(sample_2))), 1 / sqrt(n) * 10)
})


test_that("Test Exponential simulation", {
  purrr::map(2:10, function(i) {
    f <- function(rate) {
      set.seed(123)
      rexp0(n = i, rate = rate)
    }
    fs <- list(f)
    inputs <- list(
      list(rate = runif(i)),
      list(rate = runif(1))
    )
    ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-6)
    test_fs(fs, inputs, ctrl)
  })
})


test_that("Test Chi-squared simulation", {
  purrr::map(2:10, function(i) {
    f <- function(df) {
      set.seed(123)
      rchisq0(n = i, df = df)
    }
    fs <- list(f)
    inputs <- list(
      list(df = sample(10, i)),
      list(df = sample(10, 1))
    )
    ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-6)
    test_fs(fs, inputs, ctrl)
  })
})


test_that("Test Wishart simulation", {
  # (The maximum error of) Matrix samples are more sensitive than the
  # Scalar samples as they have an order more of entries (n^2 vs n).
  purrr::map(2:10, function(i) {
    # dual-dual case
    f <- function(v, M) {
      set.seed(123)
      rWishart0(v = v, M = M)
    }
    fs <- list(f)
    inputs <- list(list(v = i, M = diag(i) + crossprod(randn(i, i)) / i^2))
    ctrl <- list(display = F, err_fun = abs_err, epsilon = 1e-5)
    test_fs(fs, inputs, ctrl)
  })

  purrr::map(2:10, function(i) {
    # numeric-dual case
    f2 <- function(M) {
      set.seed(123)
      rWishart0(v = i, M = M)
    }
    fs <- list(f2)
    inputs <- list(list(M = diag(i) + crossprod(randn(i, i)) / i^2))
    ctrl <- list(display = F, err_fun = abs_err, epsilon = 1e-5)
    test_fs(fs, inputs, ctrl)
  })

  purrr::map(2:10, function(i) {
    # dual-numeric case
    m0 <- crossprod(randn(i, i))
    f3 <- function(v) {
      set.seed(123)
      rWishart0(v = v, M = m0)
    }
    fs <- list(f3)
    inputs <- list(list(v = i))
    ctrl <- list(display = F, err_fun = abs_err, epsilon = 1e-5)
    test_fs(fs, inputs, ctrl)
  })
})


test_that("Test Gamma simulation", {
  purrr::map(2:10, function(i) {
    # dual-dual case
    f <- function(shape, scale) {
      set.seed(123)
      rgamma0(n = i, shape = shape, scale = scale)
    }
    inputs <- list(
      list(shape = 3 + runif(i), scale = 3 + runif(i)),
      list(shape = 3 + runif(1), scale = 3 + runif(i)),
      list(shape = 3 + runif(i), scale = 3 + runif(1)),
      list(shape = 3 + runif(1), scale = 3 + runif(1))
    )
    ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-6)
    test_fs(list(f), inputs, ctrl)

    # dual-numeric case
    f2 <- function(shape) {
      set.seed(123)
      rgamma0(n = i, shape = shape, scale = 1)
    }
    inputs <- list(
      list(shape = 3 + runif(i)),
      list(shape = 3 + runif(1))
    )
    ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-6)
    test_fs(list(f2), inputs, ctrl)

    # dual-numeric case
    f3 <- function(scale) {
      set.seed(123)
      rgamma0(n = i, shape = 1, scale = scale)
    }
    f4 <- function(scale) {
      set.seed(123)
      rgamma0(n = i, shape = 1:i, scale = scale)
    }
    inputs <- list(
      list(scale = runif(i)),
      list(scale = runif(1))
    )
    ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-6)
    test_fs(list(f3, f4), inputs, ctrl)
  })
})
