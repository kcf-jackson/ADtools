testthat::context("Test generic functions")


check_len <- function(x, n) { testthat::expect_equal(length(x), n) }
A <- randn(5, 5)
A_dual <- dual(A, param_dim = c(length(A), 20), 1)


testthat::test_that("Test length", {
	check_len(A_dual, length(A))
})
