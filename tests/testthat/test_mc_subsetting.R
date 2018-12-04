testthat::context("Test subsetting a dual number")


check_len <- function(x, n) { testthat::expect_equal(length(x), n) }
check_dim <- function(x, n) { testthat::expect_equal(dim(x), n) }
check_mat_eq <- function(A, B) { testthat::expect_true(all(A == B)) }


testthat::test_that("Test subsetting X[i, ]", {
  A <- randn(5, 5)
  A_dual <- dual(A, param_dim = c(length(A), 20), 1)

  for (s in 1:5) {
    ha <- A_dual[s, , drop = F]
    check_dim(parent_of(ha), c(1, ncol(A)))
    check_dim(deriv_of(ha), c(ncol(A), length(A) + 20))
    check_mat_eq(A[s, , drop = F], parent_of(ha))
  }
  s <- c(1, 3)
  ha <- A_dual[s, , drop = F]
  check_dim(parent_of(ha), c(length(s), ncol(A)))
  check_dim(deriv_of(ha), c(length(s) * ncol(A), length(A) + 20))
  check_mat_eq(A[s, , drop = F], parent_of(ha))
})


testthat::test_that("Test subsetting X[, j]", {
  A <- randn(5, 5)
  A_dual <- dual(A, param_dim = c(length(A), 20), 1)

  for (s in 1:5) {
    ha <- A_dual[, s, drop = F]
    check_dim(parent_of(ha), c(nrow(A), 1))
    check_dim(deriv_of(ha), c(nrow(A), length(A) + 20))
    check_mat_eq(A[, s, drop = F], parent_of(ha))
  }
  s <- c(1, 3)
  ha <- A_dual[, s, drop = F]
  check_dim(parent_of(ha), c(nrow(A), length(s)))
  check_dim(deriv_of(ha), c(nrow(A) * length(s), length(A) + 20))
  check_mat_eq(A[, s, drop = F], parent_of(ha))
})


testthat::test_that("Test subsetting X[i,j]", {
  A <- randn(5, 5)
  A_dual <- dual(A, param_dim = c(length(A), 20), 1)

  for (m in 1:5) {
    for (n in 1:5) {
      ha <- A_dual[m, n, drop = F]
      check_dim(parent_of(ha), c(1,1))
      check_dim(deriv_of(ha), c(1, length(A) + 20))
      check_mat_eq(A[m, n, drop = F], parent_of(ha))
    }
  }
})


testthat::test_that("Test head", {
  A <- randn(5, 5)
  A_dual <- dual(A, param_dim = c(length(A), 20), 1)

  L <- nrow(A)
  testthat::expect_error(head(A_dual, 0))
  for (s in 1:(2 * L)) {
    ha <- head(A_dual, s)
    s_bdd <- min(s, nrow(A))
    check_dim(parent_of(ha), c(s_bdd, ncol(A)))
    check_dim(deriv_of(ha), c(s_bdd * ncol(A), length(A) + 20))
    testthat::expect_equal(head(A, s), parent_of(ha))
  }
  for (s in -1:-(L - 1)) {
    ha <- head(A_dual, s)
    check_dim(parent_of(ha), c(L + s, ncol(A)))
    check_dim(deriv_of(ha), c((L + s) * ncol(A), length(A) + 20))
    testthat::expect_equal(head(A, s), parent_of(ha))
  }
  for (s in -L:(-2 * L)) {
    testthat::expect_error(head(A_dual, s))
  }
})


testthat::test_that("Test tail", {
  A <- randn(5, 5)
  A_dual <- dual(A, param_dim = c(length(A), 20), 1)

  L <- nrow(A)
  testthat::expect_error(tail(A_dual, 0))
  for (s in 1:(2 * L)) {
    ta <- tail(A_dual, s)
    s_bdd <- min(s, L)
    check_dim(parent_of(ta), c(s_bdd, ncol(A)))
    check_dim(deriv_of(ta), c(s_bdd * ncol(A), length(A) + 20))
    # tail gives different row indexes, so cannot use `expect_equal` directly.
    testthat::expect_true(all(tail(A, s) - parent_of(ta) < 1e-8))
  }
  for (s in -1:-(L - 1)) {
    ha <- tail(A_dual, s)
    check_dim(parent_of(ha), c(L + s, ncol(A)))
    check_dim(deriv_of(ha), c((L + s) * ncol(A), length(A) + 20))
    testthat::expect_true(all(tail(A, s) - parent_of(ha) < 1e-8))
  }
  for (s in -L:(-2 * L)) {
    testthat::expect_error(tail(A_dual, s))
  }
})
