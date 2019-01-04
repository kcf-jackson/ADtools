# Finite differencing for testing
# @param f A function of which the derivative is seeked.
# @param x The point at which the derivative is required.
# @param h The finite differencing parameter; the size of perturbation.
finite_diff_test <- function(f, x, h = 1e-8) {
  perturbate <- function(v, h) {
    purrr::map(seq_along(v), function(i) { v[i] <- v[i] + h; v })
  }
  diff_fun <- . %>% { as.numeric((f(.) - f(x)) / h) }

  perturbate(x, h) %>%
    purrr::map(diff_fun) %>%
    do.call(cbind, .)
}


# testthat::context("Test finite differencing")
#
# approx_eq <- function(expr, output) {
#   res <- eval(quote(expr))
#   testthat::expect_lt(sum(abs(res - output)), 1e-6)
# }
#
# testthat::test_that("Test finite difference", {
#   f1 <- function(a) { 2 * a + c(1,2,3) }  # R^3 -> R^3
#   approx_eq(finite_diff_test(f1, rnorm(3)), 2 * diag(3))
#
#   f2 <- function(a) { a * c(5,4,3) }  # R^3 -> R^3
#   approx_eq(finite_diff_test(f2, rnorm(3)), diag(c(5,4,3)))
#
#   f3 <- function(a) { c(2 * a, a[1]) }  # R^3 -> R^4
#   approx_eq(finite_diff_test(f3, rnorm(3)), rbind(2 * diag(3), c(1,0,0)))
# })
