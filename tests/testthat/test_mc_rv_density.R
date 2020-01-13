context("Test differentiation of density function")

testthat::test_that("Test multivariate t distribution", {
  set.seed(123)
  p <- 10
  df <- 4
  delta <- rnorm(p)
  sigma <- diag(runif(p, min = 0.5))
  x <- matrix(rnorm(2*p, mean = 2), nrow = 2, ncol = p)

  # Check function implementation
  testthat::expect_equal(
    mvtnorm::dmvt(x, delta = delta, sigma = sigma, df = df),
    as.numeric(dmvt0(x, delta = delta, sigma = sigma, df = df))
  )

  # Auto-differentiation
  f <- function(df, delta, sigma) {
    dmvt0(x = x, df = df, delta = delta, sigma = sigma)
  }
  inputs <- list(list(df = df, delta = delta, sigma = sigma))
  ctrl <- list(display = F, err_fun = abs_err, epsilon = 1e-6)
  test_fs(list(f), inputs, ctrl)
})
