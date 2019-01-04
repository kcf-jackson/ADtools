testthat::context("Test matrix calculus")

# Helper functions
compare_FD_and_AD <- function(FD, AD, show = F) {
  rel_err <- relative_diff(FD, AD)
  if (show) {
    print(cbind(FD = as.numeric(FD), AD = as.numeric(AD)))
    print(glue::glue("Maximum relative error over all entries: {max(rel_err)} (entry: {which.max(rel_err)})"))
  }
  max(rel_err)
}


test_binary_operation <- function(fun, h = 1e-8, epsilon = 1e-6) {
	set.seed(123)
	A <- randn(3, 3)
	B <- randn(3, 3)
	param <- list(A = length(A), B = length(B))
	A_dual <- dual(A, param, 1)
	B_dual <- dual(B, param, 2)

	f1 = function(x) { fun(x, B) }
	FD_res <- finite_diff_test(f1, A, h)
	AD_res <- get_deriv(fun(A_dual, B_dual), "A")
	testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), epsilon)

	f2 = function(x) { fun(A, x) }
	FD_res <- finite_diff_test(f2, B, h)
	AD_res <- get_deriv(fun(A_dual, B_dual), "B")
	testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), epsilon)
}

testthat::test_that("Test tcrossprod (two arguments) against finite difference", {
  test_binary_operation(tcrossprod)
})

testthat::test_that("Test crossprod (two arguments) against finite difference", {
  test_binary_operation(crossprod)
})



test_unary_operation <- function(f, A, h = 1e-8, epsilon = 1e-6) {
  FD_res <- finite_diff_test(f, A, h)

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
  # The other way to look at this is that given the symmetric 3 x 3 input A,
  # the implmentation in AD discards entries 4, 7, 9, while the implmentation in FD
  # discards entries 2, 3, 6. Hence, there is an "extra" commutation matrix at the end
  # even though chol0 is merely t o chol.
  K_nn <- commutation_matrix0(nrow(A), ncol(A))
  FD_res <- K_nn %*% finite_diff_test(chol, A, h) %*% K_nn

  X <- dual(A, length(A), 1)
  AD_res <- deriv_of(f(X))

  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), epsilon)
})

testthat::test_that("Test determinant against finite difference", {
  set.seed(123)
  test_unary_operation(f = det0, A = randn(3, 3) * 10, h = 1e-7)
})
