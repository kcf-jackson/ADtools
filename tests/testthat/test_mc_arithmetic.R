testthat::context("Test scalar-dual / dual-scalar arithmetic")
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

testthat::test_that("Test ANY + dual ; dual + ANY", {
  k <- 5
  K <- randn(3, 3)
  B <- randn(3, 3)
  param_dim <- list(B = length(B))
  B <- dual(x = B, param_dim, 1)

  # scalar + dual
  testthat::expect_equal(k + parent_of(B), parent_of(k + B))
  testthat::expect_equal(parent_of(B) + k, parent_of(B + k))
  testthat::expect_equal(deriv_of(B), deriv_of(B + k))
  testthat::expect_equal(deriv_of(B), deriv_of(k + B))

  # matrix + dual
  testthat::expect_equal(K + parent_of(B), parent_of(K + B))
  testthat::expect_equal(parent_of(B) + K, parent_of(B + K))
  testthat::expect_equal(deriv_of(B), deriv_of(B + K))
  testthat::expect_equal(deriv_of(B), deriv_of(K + B))
})

testthat::test_that("Test ANY - dual ; dual - ANY", {
  k <- 5
  K <- randn(3, 3)
  B <- randn(3, 3)
  param_dim <- list(B = length(B))
  B <- dual(x = B, param_dim, 1)

  # scalar - dual
  testthat::expect_equal(k - parent_of(B), parent_of(k - B))
  testthat::expect_equal(parent_of(B) - k, parent_of(B - k))
  testthat::expect_equal(deriv_of(B), deriv_of(B - k))
  testthat::expect_equal(-deriv_of(B), deriv_of(k - B))

  # matrix - dual
  testthat::expect_equal(K - parent_of(B), parent_of(K - B))
  testthat::expect_equal(parent_of(B) - K, parent_of(B - K))
  testthat::expect_equal(deriv_of(B), deriv_of(B - K))
  testthat::expect_equal(-deriv_of(B), deriv_of(K - B))
})

testthat::test_that("Test ANY * dual ; dual * ANY", {
  k <- 5
  K <- randn(3, 3)
  B <- randn(3, 3)
  param_dim <- list(B = length(B))
  B <- dual(x = B, param_dim, 1)

  # scalar * dual
  testthat::expect_equal(k * parent_of(B), parent_of(k * B))
  testthat::expect_equal(parent_of(B) * k, parent_of(B * k))
  testthat::expect_equal(deriv_of(B) * k, deriv_of(B * k))
  testthat::expect_equal(k * deriv_of(B), deriv_of(k * B))

  # matrix * dual
  testthat::expect_equal(K * parent_of(B), parent_of(K * B))
  testthat::expect_equal(parent_of(B) * K, parent_of(B * K))
  testthat::expect_equal(deriv_of(B) * as.numeric(K), deriv_of(B * K))
  testthat::expect_equal(as.numeric(K) * deriv_of(B), deriv_of(K * B))
})

testthat::test_that("Test ANY / dual ; dual / ANY", {
  k <- 5
  K <- randn(3, 3)
  B <- randn(3, 3)
  param_dim <- list(B = length(B))
  B <- dual(x = B, param_dim, 1)

  # scalar / dual
  testthat::expect_equal(k / parent_of(B), parent_of(k / B))
  testthat::expect_equal(parent_of(B) / k, parent_of(B / k))
  testthat::expect_equal(deriv_of(B) / k, deriv_of(B / k))

  # matrix / dual
  testthat::expect_equal(K / parent_of(B), parent_of(K / B))
  testthat::expect_equal(parent_of(B) / K, parent_of(B / K))
  testthat::expect_equal(deriv_of(B) / as.numeric(K), deriv_of(B / K))
})

testthat::test_that("Test ANY %*% dual ; dual %*% ANY", {
  K <- randn(3, 3)
  B <- randn(3, 3)
  param_dim <- list(B = length(B))
  B <- dual(x = B, param_dim, 1)

  # matrix %*% dual
  abs_diff <- function(x, y) { abs(x - y) }
  testthat::expect_equal(K %*% parent_of(B), parent_of(K %*% B))
  testthat::expect_equal(parent_of(B) %*% K, parent_of(B %*% K))
  testthat::expect_true(
    all(
      abs_diff(
        (t(K) %x% diag(3)) %*% deriv_of(B),
        deriv_of(B %*% K)
      ) < 1e-8
    )
  )
  testthat::expect_true(
    all(
      abs_diff(
        (diag(3) %x% K) %*% deriv_of(B),
        deriv_of(K %*% B)
      ) < 1e-8
    )
  )
})

testthat::test_that("Test ANY %x% dual ; dual %x% ANY", {
  K <- randn(4, 4)
  B <- randn(3, 3)
  param_dim <- list(B = length(B))
  B <- dual(x = B, param_dim, 1)

  # matrix %x% dual
  testthat::expect_equal(K %x% parent_of(B), parent_of(K %x% B))
  testthat::expect_equal(parent_of(B) %x% K, parent_of(B %x% K))
})
