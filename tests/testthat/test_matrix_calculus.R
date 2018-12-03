testthat::context("Test matrix calculus")
library(Matrix)

# Helper functions
relative_diff <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  abs(x - y) / max(abs(x), abs(y))
}

compare_FD_and_AD <- function(FD, AD, show = F) {
  rel_err <- relative_diff(FD, AD)
  if (show) {
    print(cbind(as.numeric(FD), as.numeric(AD)))
    print(glue::glue("Maximum relative error over all entries: {max(rel_err)} (entry: {which.max(rel_err)})"))
  }
  max(rel_err)
}


test_unary_operation <- function(f, A, h = 1e-8, epsilon = 1e-6) {
  FD_res <- finite_diff(f, A, h)

  X <- dual(A, length(A), 1)
  AD_res <- deriv_of(f(X))

  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), epsilon)
}

testthat::test_that("Test d_solve against finite difference", {
  set.seed(123)
  test_unary_operation(f = solve, A = randn(3, 3) * 10)
})

testthat::test_that("Test d_transpose against finite difference", {
  set.seed(123)
  test_unary_operation(f = t, A = randn(3, 3) * 10)
})

testthat::test_that("Test d_tcrossprod against finite difference", {
  set.seed(123)
  test_unary_operation(f = tcrossprod, A = randn(3, 3) * 10)
})

testthat::test_that("Test crossprod against finite difference", {
  set.seed(123)
  test_unary_operation(f = crossprod, A = randn(3, 3) * 10)
})

testthat::test_that("Test d_chol against finite difference", {
  f <- chol0
  epsilon <- 1e-6
  A <- tcrossprod(randn(3,3))  # ensure positive-definite-ness
  h <- 1e-8

  # Note: The cholesky decomposition in R (`chol`)
  # 1. returns an upper-triangle matrix
  # 2. uses only the upper half of the matrix when the matrix is real
  # So suppose `chol0` returns a lower-triangle matrix, then
  # d vec(chol0(A)) / d vec(A) in AD is actually d vec(t(chol(A))) / d vec(t(A)) in FD.
  # => AD = [d vec(t(chol(A))) / d vec(A)] * [d vec(A) / d vec(t(A))]
  #       = K_nn %*% d_chol(A) %*% K_nn
  #       = K_nn %*% FD %*% K_nn
  # The other way to look at this is that given the symmetric input A,
  # the implmentation in AD discards entries 4, 7, 9, while the implmentation in FD
  # discards entries 2, 3, 6. Hence, there is an "extra" commutation matrix at the end
  # even though chol0 is merely t o chol.
  K_nn <- commutation_matrix0(nrow(A), ncol(A))
  FD_res <- K_nn %*% finite_diff(chol, A, h) %*% K_nn

  X <- dual(A, length(A), 1)
  AD_res <- deriv_of(f(X))

  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), epsilon)
})

testthat::test_that("Test determinant against finite difference", {
  set.seed(123)
  test_unary_operation(f = det0, A = randn(3, 3) * 10, h = 1e-7)
})
