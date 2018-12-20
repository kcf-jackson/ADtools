testthat::context("Test generic functions")


check_len <- function(x, n) { testthat::expect_equal(length(x), n) }
A <- randn(4, 5)
A_dual <- dual(A, param_dim = c(length(A), 20), 1)


testthat::test_that("Test length", {
	check_len(A_dual, length(A))
})


testthat::test_that("Test as.vector", {
  vec_A_dual <- as.vector(A_dual)
  expect_true(is.vector(parent_of(vec_A_dual)))
  expect_identical(deriv_of(A_dual), deriv_of(vec_A_dual))
})


testthat::test_that("Test as.matrix", {
  mat_A_dual <- as.matrix(A_dual, 2, 10)
  expect_equal(length(parent_of(A_dual)), length(parent_of(mat_A_dual)))
  expect_identical(deriv_of(A_dual), deriv_of(mat_A_dual))
})
