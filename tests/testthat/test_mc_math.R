context("Test element-wise derivative")

eq_transform <- function(x, y, f, show = F) {
  f <- purrr::compose(as.numeric, f)
  if (show) print(cbind(f(x), f(y)))
  testthat::expect_equal(f(x), f(y))
}

testthat::test_that("Test element-wise derivative", {
  sin0 <- element_wise_derivative(sin, cos)
  A <- dual(1, list(A = 1, B = 1), 1)
  eq_transform(sin0(A)@dx[1], finite_diff_test(sin, 1), identity)
})

testthat::test_that("Test sqrt of a dual object", {
  A <- dual(runif(10), list(A = 10, B = 1), 1)
  eq_transform(sqrt(A * A), A, parent_of)
  eq_transform(sqrt(A * A), A, deriv_of)
})
