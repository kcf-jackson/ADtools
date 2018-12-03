testthat::context("Test matrix sum")
library(Matrix)


A <- randn(5, 5)
A_dual <- dual(A, param_dim = list(dA = length(A), dB = 20), 1)


check_len <- function(x, n) { testthat::expect_equal(length(x), n) }
check_dim <- function(x, n) { testthat::expect_equal(dim(x), n) }
check_mat_eq <- function(A, B) { testthat::expect_true(all(A == B)) }
check_eq <- function(x, y) { testthat::expect_true(all(x == y)) }


testthat::test_that("Test sum", {
	res <- sum(A_dual)
	check_len(parent_of(res), 1)
	check_dim(deriv_of(res), c(1, ncol(deriv_of(A_dual))))
})


testthat::test_that("Test rowSums", {
	res <- rowSums(A_dual)
	check_len(parent_of(res), ncol(A))
	check_dim(deriv_of(res), c(ncol(A), ncol(deriv_of(A_dual))))
})


testthat::test_that("Test colSums", {
	res <- colSums(A_dual)
	check_len(parent_of(res), nrow(A))
	check_dim(deriv_of(res), c(nrow(A), ncol(deriv_of(A_dual))))
})


testthat::test_that("Test matrix trace", {
	res <- tr(A_dual)
	check_len(parent_of(res), 1)
	check_dim(deriv_of(res), c(1, ncol(deriv_of(A_dual))))
})
