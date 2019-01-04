testthat::test_that("Test AD source-code overloading", {
  f <- function(y, X, beta) { y - X %*% beta }
  vary <- list(beta = c(5, 6))
  fix <- list(X = matrix(1:4, 2, 2), y = c(2,3))
  res_AD <- auto_diff(f, vary = vary, fix = fix)
  res_FD <- finite_diff(f, vary = vary, fix = fix)
  testthat::expect_lt(
    max(relative_diff(get_deriv(res_AD), res_FD)), 1e-7
  )

  g <- function(X, Y) { X %*% Y }
  X <- randn(2, 2)
  Y <- randn(2, 2)
  res_AD <- auto_diff(g, vary = list(X = X, Y = Y))
  res_FD <- finite_diff(g, vary = list(X = X, Y = Y))
  testthat::expect_lt(
    max(relative_diff(get_deriv(res_AD), res_FD)), 1e-7
  )
})


testthat::test_that("Single-layer destructuring", {
  list(x, y) %<-% list(1, 2)
  testthat::expect_equal(x, 1)
  testthat::expect_equal(y, 2)

  list(x, y) %<-% list(4, 5, 6)
  testthat::expect_equal(x, 4)
  testthat::expect_equal(y, 5)

  testthat::expect_error(c(x, y) %<-% list(3))
  testthat::expect_error(c(x, y) %<-% c(3, 3))
  testthat::expect_error(list(x, y) %<-% c(3, 6))
  testthat::expect_error(list(x, y, z) %<-% list(3, 6))
})
