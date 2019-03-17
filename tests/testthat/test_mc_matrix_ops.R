testthat::context("Test solve, t, tcrossprod, crossprod, chol0, det")


# One-argument case (dual)
solve_f <- function(x) solve(x)  # change the argument name to match the other functions
fs <- list(solve_f, t, tcrossprod, crossprod, det)
inputs <- generate_inputs(2:15, lambda(list(x = randn(i, i))))
ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-6)
test_fs(fs, inputs, ctrl)


# Two-argument case (dual, dual) for tcrossprod, crossprod, solve
fs <- list(tcrossprod, crossprod, lambda(solve(x, y)))
inputs <- generate_inputs(
  2:15, lambda(list(x = 10 + randn(i, i), y = randn(i, i)))
)
ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-6)
test_fs(fs, inputs, ctrl)


# Two-argument case (dual, matrix) for tcrossprod, crossprod, solve
purrr::map(2:15, function(i) {
  m0 <- randn(i, i)
  fs <- list(
    lambda(x, solve(a = x, b = m0)),
    lambda(x, solve(a = 10 + m0, b = x)),
    lambda(x, tcrossprod(x, y = m0)),
    lambda(x, crossprod(x, y = m0))
  )
  inputs <- list( list(x = randn(i, i)) )
  ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-6)
  test_fs(fs, inputs, ctrl)
})


# The testing of Cholesky decomposition requires a special treatment.
# This is because the Cholesky decomposition in R (`chol`)
#   1. returns an upper-triangle matrix
#   2. uses only the upper half of the matrix when the matrix is real
# This does not match with the usual math notation of A = LL^T, where L is 
# a lower triangular matrix. (Note that the Cholesky-Banachiewicz algorithm
# for example only uses the lower triangular part of A to compute L, so one 
# should expect d L / d A to look like something you get by differentiating
# a lower triangular matrix w.r.t. a lower triangular matrix, i.e. the resulting
# Jacobian matrix is also lower triangular.)
# To convert the R version of chol to the usual math version of chol, we define
#   math_chol <- function(x) { t(chol(t(x))) }
# The first transpose ensures the output is lower-triangular, and the second
# ensures the lower-triangular part of the input is used. Now, math_chol
#   1. returns a lower-trianglar matrix
#   2. uses only the lower half of the matrix when the matrix is real
# and this is what we want.
# (Additional note: finite-differencing with Cholesky is not too accurate 
# due to the many floating point operations involved.)
set.seed(123)
inputs <- generate_inputs(2:15, lambda(list(x = 10 + crossprod(randn(i, i)))))  # avoid near zero.
math_chol <- function(x) t(chol(t(x)))
purrr::map(inputs, function(input) {
  x <- do.call(math_chol, input)
  FD <- finite_diff(math_chol, vary = input)
  AD <- auto_diff(chol0, vary = input)
  compare(AD@x, x, F, rel_err, 1e-8)
  compare(AD@dx, FD, F, rel_err, 1e-6)
})
