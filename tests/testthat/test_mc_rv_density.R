context("Test differentiation of density function")

testthat::test_that("Test multivariate t density", {
  set.seed(123)
  p <- 10
  df <- 4
  delta <- rnorm(p)
  sigma <- diag(runif(p, min = 0.5))
  x <- matrix(rnorm(2*p, mean = 2), nrow = 2, ncol = p)

  # Check function implementation
  testthat::expect_equal(
    mvtnorm::dmvt(x, delta = delta, sigma = sigma, df = df),
    dmvt0(x, delta = delta, sigma = sigma, df = df)
  )

  # Auto-differentiation
  f <- function(df, delta, sigma) {
    dmvt0(x = x, df = df, delta = delta, sigma = sigma)
  }
  inputs <- list(list(df = df, delta = delta, sigma = sigma))
  ctrl <- list(display = F, err_fun = abs_err, epsilon = 1e-6)
  test_fs(list(f), inputs, ctrl)
})


testthat::test_that("Test multivariate Normal density", {
  set.seed(123)
  p <- 5
  mean <- rnorm(p)
  sigma <- diag(runif(p, min = 0.5))
  n <- 20
  x <- matrix(rnorm(n*p, mean = 2), nrow = n, ncol = p)

  # Check function implementation
  testthat::expect_equal(
    mvtnorm::dmvnorm(x, mean = mean, sigma = sigma),
    dmvnorm0(x, mean = mean, sigma = sigma)
  )

  # Auto-differentiation
  f <- function(mean, sigma) {
    dmvnorm0(x = x, mean = mean, sigma = sigma)
  }
  inputs <- list(list(mean = mean, sigma = sigma))
  ctrl <- list(display = F, err_fun = abs_err, epsilon = 1e-6)
  test_fs(list(f), inputs, ctrl)
})


testthat::test_that("Test Normal density", {
  set.seed(123)
  mean <- 1
  sd <- 2
  n <- 20
  x <- rnorm(n, mean = mean, sd = sd)

  # Check function implementation
  testthat::expect_equal(
    dnorm(x, mean = mean, sd = sd),
    dnorm0(x, mean = mean, sd = sd)
  )

  # Auto-differentiation
  f <- function(mean, sd) {
    dnorm0(x = x, mean = mean, sd = sd)
  }
  inputs <- list(list(mean = mean, sd = sd))
  ctrl <- list(display = F, err_fun = abs_err, epsilon = 1e-6)
  test_fs(list(f), inputs, ctrl)
})


testthat::test_that("Test student-t density", {
  set.seed(123)
  df <- 5
  n <- 20
  x <- rt(n, df = df)

  # Check function implementation
  testthat::expect_equal(dt(x, df = df), dt0(x, df = df))

  # Auto-differentiation
  f <- function(df) dt0(x = x, df = df, log = F)
  inputs <- list(list(df = df))
  ctrl <- list(display = F, err_fun = abs_err, epsilon = 1e-6)
  test_fs(list(f), inputs, ctrl)
})


testthat::test_that("Test Gamma density", {
  set.seed(123)
  shape <- 4
  rate <- 9
  n <- 20
  x <- rgamma(n, shape = shape, rate = rate)

  # Check function implementation
  testthat::expect_equal(
    dgamma(x, shape = shape, rate = rate),
    dgamma0(x, shape = shape, rate = rate)
  )

  # Auto-differentiation
  f <- function(shape, rate) {
    dgamma0(x = x, shape = shape, rate = rate)
  }
  inputs <- list(list(shape = shape, rate = rate))
  ctrl <- list(display = F, err_fun = abs_err, epsilon = 1e-6)
  test_fs(list(f), inputs, ctrl)
})


testthat::test_that("Test Chi-squared density", {
  set.seed(123)
  df <- 5
  n <- 20
  x <- rchisq(n, df = df)

  # Check function implementation
  testthat::expect_equal(dchisq(x, df = df), dchisq0(x, df = df))

  # Auto-differentiation
  f <- function(df) dchisq0(x = x, df = df)
  inputs <- list(list(df = df))
  ctrl <- list(display = F, err_fun = abs_err, epsilon = 1e-6)
  test_fs(list(f), inputs, ctrl)
})


testthat::test_that("Test Exponential density", {
  set.seed(123)
  rate <- 5
  n <- 20
  x <- rexp(n, rate = rate)

  # Check function implementation
  testthat::expect_equal(
    dexp(x, rate = rate), dexp0(x, rate = rate)
  )

  # Auto-differentiation
  f <- function(rate) dexp0(x = x, rate = rate)
  inputs <- list(list(rate = rate))
  ctrl <- list(display = F, err_fun = abs_err, epsilon = 1e-6)
  test_fs(list(f), inputs, ctrl)
})


testthat::test_that("Test Wishart density", {
  set.seed(123)
  df <- 12
  Sigma <- crossprod(randn(4, 4))
  n <- 20
  x <- rWishart(n, df = df, Sigma = Sigma)
  y <- purrr::map(seq(dim(x)[3]), ~x[,,.x])
  # Check function implementation
  testthat::expect_equal(
    purrr::map_dbl(y, ~MCMCpack::dwish(.x, v = df, S = Sigma)),
    purrr::map_dbl(y, ~dWishart0(.x, v = df, M = Sigma))
  )

  # Auto-differentiation
  f <- function(df, Sigma) {
    dWishart0(y[[1]], v = df, M = Sigma)
  }
  inputs <- list(list(df = df, Sigma = Sigma))
  ctrl <- list(display = F, err_fun = abs_err, epsilon = 1e-6)
  test_fs(list(f), inputs, ctrl)
})


testthat::test_that("Test inverse Wishart density", {
  set.seed(123)
  df <- 12
  Sigma <- crossprod(randn(4, 4))
  n <- 20
  x <- purrr::map(1:n, ~MCMCpack::riwish(v = df, S = Sigma))

  # Check function implementation
  testthat::expect_equal(
    purrr::map_dbl(x, ~MCMCpack::diwish(.x, v = df, S = Sigma)),
    purrr::map_dbl(x, ~dinvWishart0(.x, v = df, M = Sigma))
  )

  # Auto-differentiation
  f <- function(df, Sigma) {
    dinvWishart0(x[[1]], v = df, M = Sigma)
  }
  inputs <- list(list(df = df, Sigma = Sigma))
  ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-5)
  test_fs(list(f), inputs, ctrl)
})
