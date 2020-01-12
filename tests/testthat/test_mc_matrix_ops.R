testthat::context("Test solve, t, tcrossprod, crossprod, chol0, det")

set.seed(1234)


# Functions without S3 declarations need to be wrapped
wrap <- function(fs) purrr::map(fs, ~lambda(.x(x)))
wrap2 <- function(fs) purrr::map(fs, ~lambda(.x(x, y)))


# One-argument case (dual) ================================================
testthat::test_that("solve", {
  solve_f <- function(x) solve(x)  # change the argument name to match the other functions
  fs <- list(solve_f, t)
  inputs <- generate_inputs(2:15, lambda(list(x = 10 * diag(i) + randn(i, i))))
  ctrl <- list(display = F, err_fun = abs_err, epsilon = 1e-7)
  test_fs(fs, inputs, ctrl)
})


testthat::test_that("tcrossprod, crossprod", {
  fs <- list(tcrossprod, crossprod)
  inputs <- generate_inputs(2:15, lambda(list(x = randn(i, i))))
  ctrl <- list(display = F, err_fun = abs_err, epsilon = 1e-6)
  test_fs(wrap(fs), inputs, ctrl)
})


testthat::test_that("det", {
  fs <- list(det)
  eigen_to_matrix <- function(eig0) {
    n <- length(eig0)
    Q <- qr.Q(qr(randu(n, n)))
    t(Q) %*% diag(eig0) %*% Q
  }
  inputs <- generate_inputs(2:25, lambda(list(
    x = eigen_to_matrix(runif(i, min = 0.1, max = 1))
  )))
  ctrl <- list(display = F, err_fun = abs_err, epsilon = 1e-7)
  test_fs(wrap(fs), inputs, ctrl)
})


testthat::test_that("chol0", {
  fs <- list(chol0)  # requires positive semi-definite input
  inputs <- generate_inputs(2:25, lambda(list(x = diag(i) + crossprod(randn(i, i)))))
  ctrl <- list(display = F, err_fun = abs_err, epsilon = 1e-6)
  test_fs(fs, inputs, ctrl)
})


# Two-argument case (dual, dual) for tcrossprod, crossprod, solve =========
testthat::test_that("tcrossprod, crossprod", {
  fs <- list(tcrossprod, crossprod)
  inputs <- generate_inputs(
    2:15, lambda(list(x = randn(i, i), y = randn(i, i)))
  )
  ctrl <- list(display = F, err_fun = abs_err, epsilon = 1e-6)
  test_fs(wrap2(fs), inputs, ctrl)
})


testthat::test_that("solve", {
  fs <- list(lambda(solve(x, y)))
  inputs <- generate_inputs(
    2:15, lambda(list(x = 10 * diag(i) + randn(i, i), y = randn(i, i)))
  )
  ctrl <- list(display = F, err_fun = abs_err, epsilon = 1e-7)
  test_fs(fs, inputs, ctrl)
})


# Two-argument case (dual, matrix) for tcrossprod, crossprod, solve =====
testthat::test_that("tcrossprod, crossprod", {
  purrr::map(2:15, function(i) {
    m0 <- randn(i, i)
    fs <- list(
      lambda(x, tcrossprod(x, y = m0)),
      lambda(x, crossprod(x, y = m0))
    )
    inputs <- list( list(x = randn(i, i)) )
    ctrl <- list(display = F, err_fun = abs_err, epsilon = 1e-6)
    test_fs(wrap(fs), inputs, ctrl)
  })
})


testthat::test_that("solve", {
  purrr::map(2:15, function(i) {
    m0 <- randn(i, i)
    fs <- list(
      lambda(x, solve(a = x, b = m0)),
      lambda(x, solve(a = 10 * diag(i) + m0, b = x))
    )
    inputs <- list( list(x = 10 * diag(i) + randn(i, i)) )
    ctrl <- list(display = F, err_fun = abs_err, epsilon = 1e-6)
    test_fs(fs, inputs, ctrl)
  })
})
