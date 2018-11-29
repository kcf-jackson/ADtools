testthat::context("Test subsetting a dual number")


check_len <- function(x, n) { testthat::expect_equal(length(x), n) }
check_dim <- function(x, n) { testthat::expect_equal(dim(x), n) }
A <- randn(5, 5)
A_dual <- dual(A, param_dim = c(length(A), 20), 1)


testthat::test_that("head", {
  for (s in 1:5) {
    ha <- head(A_dual, s)
    check_dim(parent_of(ha), c(s, ncol(A)))
    check_dim(deriv_of(ha), c(s * ncol(A), length(A) + 20))
    testthat::expect_equal(head(A, s), parent_of(ha))
  }
})


testthat::test_that("tail", {
  for (s in 1:5) {
    ta <- tail(A_dual, s)
    check_dim(parent_of(ta), c(s, ncol(A)))
    check_dim(deriv_of(ta), c(s * ncol(A), length(A) + 20))
    # tail gives different row indexes, so cannot use `expect_equal` directly.
    testthat::expect_true(all(tail(A, s) == parent_of(ta)))
  }
})
