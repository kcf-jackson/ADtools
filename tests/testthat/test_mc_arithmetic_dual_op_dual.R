testthat::context("Test dual-dual arithmetic")

# Helper functions
compare_FD_and_AD <- function(FD, AD, show = F) {
  rel_err <- relative_diff(FD, AD)
  if (show) {
    print(cbind(FD = as.numeric(FD), AD = as.numeric(AD)))
    print(glue::glue("Maximum relative error over all entries: {max(rel_err)} (entry: {which.max(rel_err)})"))
  }
  max(rel_err)
}

# Actual tests
h <- 1e-8
A <- randn(3,3)
B <- randn(3,3)
k <- 5
z <- 8
param_dim <- list(A = length(A), B = length(B), k = length(k), z = length(z))
A_dual <- dual(x = A, param_dim, 1)
B_dual <- dual(x = B, param_dim, 2)
k_dual <- dual(x = k, param_dim, 3)
z_dual <- dual(x = z, param_dim, 4)


testthat::test_that("Test dual + dual", {
  f1 <- function(x) { x + B }
  f2 <- function(x) { A + x }
  f3 <- function(x) { x + z }
  f4 <- function(x) { k + x }

  # Matrix + Matrix
  FD_res <- finite_diff(f1, A, h)
  AD_res <- get_deriv(A_dual + B_dual, "A")
  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

  FD_res <- finite_diff(f2, B, h)
  AD_res <- get_deriv(A_dual + B_dual, "B")
  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

  # Matrix + Scalar / Scalar + Matrix
  FD_res <- finite_diff(f1, k, h)
  AD_res <- get_deriv(k_dual + B_dual, "k")
  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

  FD_res <- finite_diff(f2, z, h)
  AD_res <- get_deriv(A_dual + z_dual, "z")
  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

  # Scalar + Scalar
  FD_res <- finite_diff(f3, k, h)
  AD_res <- get_deriv(k_dual + z_dual, "k")
  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

  FD_res <- finite_diff(f4, z, h)
  AD_res <- get_deriv(k_dual + z_dual, "z")
  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)
})


testthat::test_that("Test dual - dual", {
  f1 <- function(x) { x - B }
  f2 <- function(x) { A - x }
  f3 <- function(x) { x - z }
  f4 <- function(x) { k - x }

  # Matrix - Matrix
  FD_res <- finite_diff(f1, A, h)
  AD_res <- get_deriv(A_dual - B_dual, "A")
  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

  FD_res <- finite_diff(f2, B, h)
  AD_res <- get_deriv(A_dual - B_dual, "B")
  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

  # Matrix - Scalar / Scalar - Matrix
  FD_res <- finite_diff(f1, k, h)
  AD_res <- get_deriv(k_dual - B_dual, "k")
  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

  FD_res <- finite_diff(f2, z, h)
  AD_res <- get_deriv(A_dual - z_dual, "z")
  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

  # Scalar - Scalar
  FD_res <- finite_diff(f3, k, h)
  AD_res <- get_deriv(k_dual - z_dual, "k")
  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

  FD_res <- finite_diff(f4, z, h)
  AD_res <- get_deriv(k_dual - z_dual, "z")
  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)
})


testthat::test_that("Test dual * dual", {
  f1 <- function(x) { x * B }
  f2 <- function(x) { A * x }
  f3 <- function(x) { x * z }
  f4 <- function(x) { k * x }

  # Matrix * Matrix
  FD_res <- finite_diff(f1, A, h)
  AD_res <- get_deriv(A_dual * B_dual, "A")
  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

  FD_res <- finite_diff(f2, B, h)
  AD_res <- get_deriv(A_dual * B_dual, "B")
  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

  # Matrix * Scalar / Scalar * Matrix
  FD_res <- finite_diff(f1, k, h)
  AD_res <- get_deriv(k_dual * B_dual, "k")
  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

  FD_res <- finite_diff(f2, z, h)
  AD_res <- get_deriv(A_dual * z_dual, "z")
  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

  # Scalar * Scalar
  FD_res <- finite_diff(f3, k, h)
  AD_res <- get_deriv(k_dual * z_dual, "k")
  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

  FD_res <- finite_diff(f4, z, h)
  AD_res <- get_deriv(k_dual * z_dual, "z")
  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)
})


testthat::test_that("Test dual / dual", {
  f1 <- function(x) { x / B }
  f2 <- function(x) { A / x }
  f3 <- function(x) { x / z }
  f4 <- function(x) { k / x }

  # Matrix / Matrix
  FD_res <- finite_diff(f1, A, h)
  AD_res <- get_deriv(A_dual / B_dual, "A")
  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

  FD_res <- finite_diff(f2, B, h)
  AD_res <- get_deriv(A_dual / B_dual, "B")
  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

  # Matrix / Scalar / Scalar / Matrix
  FD_res <- finite_diff(f1, k, h)
  AD_res <- get_deriv(k_dual / B_dual, "k")
  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

  FD_res <- finite_diff(f2, z, h)
  AD_res <- get_deriv(A_dual / z_dual, "z")
  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

  # Scalar / Scalar
  FD_res <- finite_diff(f3, k, h)
  AD_res <- get_deriv(k_dual / z_dual, "k")
  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

  FD_res <- finite_diff(f4, z, h)
  AD_res <- get_deriv(k_dual / z_dual, "z")
  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)
})


testthat::test_that("Test dual %*% dual", {
  f <- function(x) { x %*% B }
  FD_res <- finite_diff(f, A, h)
  AD_res <- get_deriv(A_dual %*% B_dual, "A")
  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

  f2 <- function(x) { A %*% x }
  FD_res <- finite_diff(f2, B, h)
  AD_res <- get_deriv(A_dual %*% B_dual, "B")
  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)
})


testthat::test_that("Test dual %x% dual", {
  f <- function(x) { x %x% B }
  FD_res <- finite_diff(f, A, h)
  AD_res <- get_deriv(A_dual %x% B_dual, "A")
  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

  f2 <- function(x) { A %x% x }
  FD_res <- finite_diff(f2, B, h)
  AD_res <- get_deriv(A_dual %x% B_dual, "B")
  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)
})
