testthat::context("Test binding two dual numbers")


check_len <- function(x, n) { testthat::expect_equal(length(x), n) }
check_dim <- function(x, n) { testthat::expect_equal(dim(x), n) }
A <- randn(5, 5)
B <- randn(5, 10)
C <- randn(10, 5)
A_dual <- dual(A, param_dim = c(length(A), length(B), length(C)), 1)
B_dual <- dual(B, param_dim = c(length(A), length(B), length(C)), 2)
C_dual <- dual(C, param_dim = c(length(A), length(B), length(C)), 3)


testthat::test_that("Test cbinding two dual numbers", {
	x <- cbind2(A_dual, B_dual)
	check_dim(parent_of(x), c(nrow(A), ncol(A) + ncol(B)))
	check_dim(
	  deriv_of(x),
		c(length(A) + length(B), length(A) + length(B) + length(C))
	)
})


testthat::test_that("Test rbinding two dual numbers", {
	x <- rbind2(A_dual, C_dual)
	check_dim(parent_of(x), c(nrow(A) + nrow(C), ncol(A)))
	check_dim(
	  deriv_of(x),
		c(length(A) + length(C), length(A) + length(B) + length(C))
	)
})
