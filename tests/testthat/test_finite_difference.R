context("Test finite differencing")

testthat::test_that("Test finite difference main function", {
  f <- function(X, y) { X %*% y }
  set.seed(123)
  for (i in 1:10) {
    X <- randn(2, 2)
    y <- rnorm(2)
    expect_lt(
      max(relative_diff(
        finite_diff(f, list(X = X, y = y)),
        cbind(t(y %x% diag(2)), X)
      )),
      1e-6
    )
  }
})
