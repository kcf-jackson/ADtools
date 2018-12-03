testthat::context("Test finite differencing")


approx_eq <- function(expr, output) {
  res <- eval(quote(expr))
  testthat::expect_lt(sum(abs(res - output)), 1e-6)
}


testthat::test_that("Test finite difference", {
  f1 <- function(a) { 2 * a + c(1,2,3) }  # R^3 -> R^3
  approx_eq(finite_diff(f1, rnorm(3)), 2 * diag(3))

  f2 <- function(a) { a * c(5,4,3) }  # R^3 -> R^3
  approx_eq(finite_diff(f2, rnorm(3)), diag(c(5,4,3)))

  f3 <- function(a) { c(2 * a, a[1]) }  # R^3 -> R^4
  approx_eq(finite_diff(f3, rnorm(3)), rbind(2 * diag(3), c(1,0,0)))
})
