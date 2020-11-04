context("Test finite differencing")

testthat::test_that("Test finite difference main function", {
  f <- function(X, y) { X %*% y }
  set.seed(123)
  for (i in 1:10) {
    X <- randn(2, 2)
    y <- rnorm(2)
    compare(
      finite_diff(f, at = list(X = X, y = y)),
      cbind(t(y %x% diag(2)), X),
      display = F,
      err_fun = rel_err,
      epsilon = 1e-6
    )

    compare(
      finite_diff(f, at = list(X = X, y = y), method = "central"),
      cbind(t(y %x% diag(2)), X),
      display = F,
      err_fun = rel_err,
      epsilon = 1e-6
    )

    testthat::expect_error(
      finite_diff(f, at = list(X = X, y = y), method = "ERROR")
    )
  }
})

testthat::test_that("Test finite difference with randomness", {
  f <- function(k) { rnorm(5) * k }
  set.seed(123)
  expect_res <- rnorm(5)
  compare(
    finite_diff(f, at = list(k = 5), seed = 123),
    expect_res,
    display = F,
    err_fun = rel_err,
    epsilon = 1e-6
  )
})
