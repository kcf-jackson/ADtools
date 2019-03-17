testthat::context("Test generic functions")

testthat::test_that("Test length, dim, nrow, ncol", {
  A <- randn(4, 5)
  A_dual <- dual(A, param_dim = c(length(A), 20), 1)
  testthat::expect_equal(length(A_dual), length(A))
  testthat::expect_equal(dim(A_dual), dim(A))
  testthat::expect_equal(nrow(A_dual), nrow(A))
  testthat::expect_equal(ncol(A_dual), ncol(A))
})
