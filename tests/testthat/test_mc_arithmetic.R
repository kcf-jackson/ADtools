testthat::context("Test scalar-dual / dual-scalar arithmetic")

equal_up_to_type <- function(x, y) {
  type_match <- function(x) {
    x@dx <- as.matrix(deriv_of(x));
    x
  }
  testthat::expect_equal(type_match(x), type_match(y))
}

testthat::test_that("Test ANY + dual ; dual + ANY", {
  k <- 5
  K <- randn(3, 3)
  b <- 3
  b <- dual(x = b, list(b = length(b), 1), ind = 1)
  B <- randn(3, 3)
  B <- dual(x = B, list(B = length(B), 1), ind = 1)

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

  # scalar + dual
  testthat::expect_equal(k + parent_of(b), parent_of(k + b))
  testthat::expect_equal(parent_of(b) + k, parent_of(b + k))
  testthat::expect_equal(deriv_of(b), deriv_of(b + k))
  testthat::expect_equal(deriv_of(b), deriv_of(k + b))

  # matrix + dual
  testthat::expect_equal(K + parent_of(b), parent_of(K + b))
  testthat::expect_equal(parent_of(b) + K, parent_of(b + K))
  testthat::expect_equal(deriv_of(b), deriv_of(b + K)[1, , drop = F])
  testthat::expect_equal(deriv_of(b), deriv_of(K + b)[1, , drop = F])

  # Additional tests - should be consistent with dual-dual arithmetic
  k_dual <- dual(k, list(b = length(b), 1), -1)
  testthat::expect_equal(k + b, k_dual + b)
  testthat::expect_equal(b + k, b + k_dual)

  k_dual <- dual(k, list(B = length(B), 1), -1)
  testthat::expect_equal(k + B, k_dual + B)
  testthat::expect_equal(B + k, B + k_dual)

  K_dual <- dual(K, list(b = length(b), 1), -1)
  equal_up_to_type(K + b, K_dual + b)
  equal_up_to_type(b + K, b + K_dual)

  K_dual <- dual(K, list(B = length(B), 1), -1)
  testthat::expect_equal(K + B, K_dual + B)
  testthat::expect_equal(B + K, B + K_dual)
})

testthat::test_that("Test ANY - dual ; dual - ANY", {
  k <- 5
  K <- randn(3, 3)
  b <- 3
  B <- randn(3, 3)
  b <- dual(x = b, list(b = length(b)), 1)
  B <- dual(x = B, list(B = length(B)), 1)

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

  # scalar - dual
  testthat::expect_equal(k - parent_of(b), parent_of(k - b))
  testthat::expect_equal(parent_of(b) - k, parent_of(b - k))
  testthat::expect_equal(deriv_of(b), deriv_of(b - k))
  testthat::expect_equal(-deriv_of(b), deriv_of(k - b))

  # matrix - dual
  testthat::expect_equal(K - parent_of(b), parent_of(K - b))
  testthat::expect_equal(parent_of(b) - K, parent_of(b - K))
  testthat::expect_equal(deriv_of(b), deriv_of(b - K)[1, , drop = F])
  testthat::expect_equal(-deriv_of(b), deriv_of(K - b)[1, , drop = F])

  # Additional tests - should be consistent with dual-dual arithmetic
  k_dual <- dual(k, list(b = length(b)), -1)
  testthat::expect_equal(k - b, k_dual - b)
  testthat::expect_equal(b - k, b - k_dual)

  k_dual <- dual(k, list(B = length(B)), -1)
  testthat::expect_equal(k - B, k_dual - B)
  testthat::expect_equal(B - k, B - k_dual)

  K_dual <- dual(K, list(b = length(b)), -1)
  equal_up_to_type(K - b, K_dual - b)
  equal_up_to_type(b - K, b - K_dual)

  K_dual <- dual(K, list(B = length(B)), -1)
  testthat::expect_equal(K - B, K_dual - B)
  testthat::expect_equal(B - K, B - K_dual)
})

testthat::test_that("Test negation of a dual number", {
  b <- 4
  b <- dual(x = b, list(b = length(b)), 1)
  testthat::expect_equal(parent_of(-b), -parent_of(b))
  testthat::expect_equal(-parent_of(b), parent_of(-b))
  testthat::expect_equal(deriv_of(-b), -deriv_of(b))
  testthat::expect_equal(-deriv_of(b), deriv_of(-b))

  B <- randn(3, 3)
  B <- dual(x = B, list(B = length(B)), 1)
  testthat::expect_equal(parent_of(-B), -parent_of(B))
  testthat::expect_equal(-parent_of(B), parent_of(-B))
  testthat::expect_equal(deriv_of(-B), -deriv_of(B))
  testthat::expect_equal(-deriv_of(B), deriv_of(-B))
})

testthat::test_that("Test ANY * dual ; dual * ANY", {
  k <- 5
  K <- randn(3, 3)
  b <- 3
  B <- randn(3, 3)
  b <- dual(x = b, list(b = length(b), other = 2), 1)
  B <- dual(x = B, list(B = length(B), other = 2), 1)

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

  # scalar * dual
  testthat::expect_equal(k * parent_of(b), parent_of(k * b))
  testthat::expect_equal(parent_of(b) * k, parent_of(b * k))
  testthat::expect_equal(deriv_of(b) * k, deriv_of(b * k))
  testthat::expect_equal(k * deriv_of(b), deriv_of(k * b))

  # matrix * dual
  testthat::expect_equal(K * parent_of(b), parent_of(K * b))
  testthat::expect_equal(parent_of(b) * K, parent_of(b * K))
  testthat::expect_equal(as.numeric(K) %*% deriv_of(b), deriv_of(b * K))
  testthat::expect_equal(as.numeric(K) %*% deriv_of(b), deriv_of(K * b))

  # Additional tests - should be consistent with dual-dual arithmetic
  k_dual <- dual(k, list(b = length(b), other = 2), -1)
  testthat::expect_equal(k * b, k_dual * b)
  testthat::expect_equal(b * k, b * k_dual)

  k_dual <- dual(k, list(B = length(B), other = 2), -1)
  testthat::expect_equal(k * B, k_dual * B)
  testthat::expect_equal(B * k, B * k_dual)

  K_dual <- dual(K, list(b = length(b), other = 2), -1)
  equal_up_to_type(K * b, K_dual * b)
  equal_up_to_type(b * K, b * K_dual)

  K_dual <- dual(K, list(B = length(B), other = 2), -1)
  testthat::expect_equal(K * B, K_dual * B)
  testthat::expect_equal(B * K, B * K_dual)
})

testthat::test_that("Test ANY / dual ; dual / ANY", {
  k <- 5
  K <- randn(3, 3)
  b <- 3
  B <- randn(3, 3)
  b <- dual(x = b, list(b = length(b)), 1)
  B <- dual(x = B, list(B = length(B)), 1)

  # scalar / dual
  testthat::expect_equal(k / parent_of(B), parent_of(k / B))
  testthat::expect_equal(parent_of(B) / k, parent_of(B / k))
  testthat::expect_equal(deriv_of(B) / k, deriv_of(B / k))

  # matrix / dual
  testthat::expect_equal(K / parent_of(B), parent_of(K / B))
  testthat::expect_equal(parent_of(B) / K, parent_of(B / K))
  testthat::expect_equal(deriv_of(B) / as.numeric(K), deriv_of(B / K))

  # scalar / dual
  testthat::expect_equal(k / parent_of(b), parent_of(k / b))
  testthat::expect_equal(parent_of(b) / k, parent_of(b / k))
  testthat::expect_equal(deriv_of(b) / k, deriv_of(b / k))

  # matrix / dual
  testthat::expect_equal(K / parent_of(b), parent_of(K / b))
  testthat::expect_equal(parent_of(b) / K, parent_of(b / K))
  testthat::expect_equal(as.numeric(1 / K) %*% deriv_of(b), deriv_of(b / K))

  # Additional tests - should be consistent with dual-dual arithmetic
  k_dual <- dual(k, list(b = length(b)), -1)
  equal_up_to_type(k / b, k_dual / b)
  testthat::expect_equal(b / k, b / k_dual)

  k_dual <- dual(k, list(B = length(B)), -1)
  testthat::expect_equal(k / B, k_dual / B)
  testthat::expect_equal(B / k, B / k_dual)

  K_dual <- dual(K, list(b = length(b)), -1)
  equal_up_to_type(K / b, K_dual / b)
  equal_up_to_type(b / K, b / K_dual)

  K_dual <- dual(K, list(B = length(B)), -1)
  equal_up_to_type(K / B, K_dual / B)
  testthat::expect_equal(B / K, B / K_dual)
})


testthat::test_that("Test dual ^ ANY", {
  k <- 5
  K <- randn(3, 3)
  b <- 3
  B <- randn(3, 3)
  b_dual <- dual(x = b, list(b = length(b)), 1)
  B_dual <- dual(x = B, list(B = length(B)), 1)

  testthat::expect_equal(parent_of(b_dual ^ k), parent_of(b_dual) ^ k)
  testthat::expect_equal(parent_of(B_dual ^ k), parent_of(B_dual) ^ k)

  # Compare with AD
  abs_rel_err <- function(x, y) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    sum(abs(x - y) / max(x, y))
  }
  .f <- function(x) { x^k }

  AD_res <- deriv_of(b_dual ^ k)
  FD_res <- finite_diff(.f, b)
  testthat::expect_lt(abs_rel_err(AD_res, FD_res), 1e-6)

  AD_res <- deriv_of(B_dual ^ k)
  FD_res <- finite_diff(.f, B)
  testthat::expect_lt(abs_rel_err(AD_res, FD_res), 1e-6)

  testthat::expect_error(b_dual ^ K)
  testthat::expect_error(B_dual ^ K)

  # Test taking to the power 0
  testthat::expect_true(all(parent_of(b_dual ^ 0) == 1))
  testthat::expect_true(all(parent_of(B_dual ^ 0) == 1))
  testthat::expect_true(all(deriv_of(b_dual ^ 0) == 0))
  testthat::expect_true(all(deriv_of(B_dual ^ 0) == 0))

  # Check consistence with sqrt
  eq <- function(x, y, transform) {
    testthat::expect_equal(transform(x), transform(y))
  }
  eq_up_to_eps <- function(x, y, transform, epsilon) {
    testthat::expect_lt(abs_rel_err(transform(x), transform(y)), epsilon)
  }

  B_dual@x <- abs(parent_of(B_dual))
  eq(sqrt(B_dual), B_dual^0.5, parent_of)
  eq_up_to_eps(sqrt(B_dual), B_dual^0.5, deriv_of, 1e-10)
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

  # Additional tests - should be consistent with dual-dual arithmetic
  K_dual <- dual(x = K, param_dim, -1)
  equal_up_to_type(B %*% K, B %*% K_dual)
  equal_up_to_type(K %*% B, K_dual %*% B)
})

testthat::test_that("Test ANY %x% dual ; dual %x% ANY", {
  K <- randn(4, 4)
  B <- randn(3, 3)
  param_dim <- list(B = length(B))
  B <- dual(x = B, param_dim, 1)

  # matrix %x% dual
  testthat::expect_equal(K %x% parent_of(B), parent_of(K %x% B))
  testthat::expect_equal(parent_of(B) %x% K, parent_of(B %x% K))

  # Additional tests - should be consistent with dual-dual arithmetic
  K_dual <- dual(x = K, param_dim, -1)
  equal_up_to_type(B %x% K, B %x% K_dual)
  equal_up_to_type(K %x% B, K_dual %x% B)
})
