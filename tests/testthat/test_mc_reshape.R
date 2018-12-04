testthat::context("Test reshaping dual numbers")


A <- randn(5, 5)
A_dual <- dual(A, param_dim = list(dA = length(A), dB = 20), 1)


check_len <- function(x, n) { testthat::expect_equal(length(x), n) }
check_dim <- function(x, n) { testthat::expect_equal(dim(x), n) }
check_mat_eq <- function(A, B) { testthat::expect_true(all(A == B)) }
check_eq <- function(x, y) { testthat::expect_true(all(x == y)) }


testthat::test_that("Test vectorisation and half-vectorisation", {
  # vec operator
  check_len(vec(A), length(A))  # basic check

  res <- vec(A_dual)
  check_dim(parent_of(res), c(length(A), 1))
  check_eq(parent_of(res)[2], A[2,1])

  # vech operator
  check_len(vech(A), sum(upper.tri(A, diag = T)))  # basic check
  testthat::expect_error(vech(randn(3, 5)))        # basic check

  res <- vech(A_dual)
  check_dim(parent_of(res), c(0.5 * nrow(A) * (nrow(A) + 1), 1))
  check_eq(parent_of(res)[2], A[2,1])
})


testthat::test_that("Test diag", {
  res <- diag(A_dual)
  check_eq(parent_of(res), diag(A))
  check_dim(deriv_of(res), c(nrow(A), ncol(deriv_of(A_dual))))

  res <- diag(res)
  check_eq(parent_of(res), diag(diag(A)))
  check_dim(deriv_of(res), dim(deriv_of(A_dual)))
})
