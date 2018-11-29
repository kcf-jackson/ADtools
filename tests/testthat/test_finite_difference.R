testthat::context("Test finite differencing")


unit_test <- function(expr, output) {
  eval(quote(expr)) %T>%
    print() %>%
    testthat::expect_equal(., output)
}


testthat::test_that("Test finite difference", {
  f1 <- function(a) { 2 * a + c(1,2,3) }  # R^3 -> R^3
  unit_test(finite_diff(f1, rnorm(3)), 2 * diag(3))

  f2 <- function(a) { a * c(5,4,3) }  # R^3 -> R^3
  unit_test(finite_diff(f2, rnorm(3)), diag(c(5,4,3)))

  f3 <- function(a) { c(2 * a, a[1]) }  # R^3 -> R^4
  unit_test(finite_diff(f3, rnorm(3)), rbind(2 * diag(3), c(1,0,0)))
})
