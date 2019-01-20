context("Test random variable simulation")

# Helper functions
compare_FD_and_AD <- function(FD, AD, show = F) {
  AD <- AD[seq(nrow(FD)), seq(ncol(FD))]
  rel_err <- relative_diff(FD, AD)
  if (show) {
    print(cbind(FD = as.numeric(FD), AD = as.numeric(AD)))
    print(glue::glue("Maximum relative error over all entries: {max(rel_err)} (entry: {which.max(rel_err)})"))
  }
  max(rel_err)
}

eq_transform <- function(x, y, f, show = F) {
  f <- purrr::compose(as.numeric, f)
  if (show) print(cbind(f(x), f(y)))
  testthat::expect_equal(f(x), f(y))
}

aeq_transform <- function(x, y, f, show = F) {
  f <- purrr::compose(as.numeric, f)
  diff0 <- sum(relative_diff(f(x), f(y)))
  if (show) print(cbind(f(x), f(y)))
  testthat::expect_lt(diff0, 1e-8)
}

test_that("Test univariate normal simulation", {
  mu <- 20
  sigma <- 3
  for (n in 1:10) {
    # Check Derivative
    # AD
    mean0 <- dual(mu, list(A = 1, B = 1, C = 3), 1)
    sd0 <- dual(sigma, list(A = 1, B = 1, C = 3), 2)
    set.seed(123)
    AD_res <- rnorm0(n, mean = mean0, sd = sd0)@dx
    # FD
    f <- function(param) {
      set.seed(123)
      rnorm0(n, mean = param[1], sd = param[2])
    }
    FD_res <- finite_diff_test(f, c(mu, sigma))
    testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)
    # Check sample
    set.seed(123)
    AD_sample <- rnorm0(n, mean = mean0, sd = sd0)@x
    FD_sample <- f(c(mu, sigma))
    testthat::expect_lt(sum(abs(AD_sample - FD_sample)), 1e-8)
  }
  # Test cases where mean0 and sd0 have length > 1
  set.seed(123)
  param_dim <- list(A = 1, B = 2, C = 1)
  mu <- rnorm(10, 0, 10)
  sigma <- runif(10, 1, 10)

  mean0 <- dual(mu, param_dim, -1)
  mean0@dx[,1] <- 1
  sd0 <- dual(sigma, param_dim, -1)
  sd0@dx[,2] <- 1

  set.seed(123)
  AD_res <- rnorm0(10, mean0, sd0)
  set.seed(123)
  AD_res_2 <- mapreduce(1:10, ~rnorm0(1, mean0[.x], sd0[.x]), rbind2)
  eq_transform(AD_res, AD_res_2, parent_of)
  aeq_transform(AD_res, AD_res_2, deriv_of)

  # Test when mean0 and n have unequal length
  mean0 <- dual(c(-10, 5, 10), param_dim, -1)
  sd0 <- dual(c(5, 1, 3), param_dim, -1)
  testthat::expect_error(rnorm0(2, mean0, sd0))

  # Test when sd0 and n have unequal length
  mean0 <- dual(c(-10, 10), param_dim, -1)
  sd0 <- dual(c(5, 1, 3), param_dim, -1)
  testthat::expect_error(rnorm0(2, mean0, sd0))
})

test_that("Test multivariate normal simulation", {
  mu <- matrix(c(1, 2))
  sigma <- diag(2)
  for (n in 1:10) {
    # Check Derivative
    # AD
    param <- list(A = length(mu), B = length(sigma), C = 2)
    mean0 <- dual(mu, param, 1)
    sd0 <- dual(sigma, param, 2)
    set.seed(123)
    AD_res <- rmvnorm0(n, mean0, sd0)@dx
    # FD
    f <- function(param) {
      set.seed(123)
      Z <- matrix(rnorm(n * nrow(param)), nrow = nrow(param))
      param[,1] + chol0(param[,-1]) %*% Z
    }
    # see test_matrix_calculus.R: `d_chol` for more detail why an
    # extra commutation matrix is needed.
    K_nn <- commutation_matrix0(nrow(sigma), ncol(sigma))
    FD_res <- finite_diff_test(f, cbind(mu, sigma))
    FD_res[, -seq_along(mu)] <- as.matrix(FD_res[, -seq_along(mu)] %*% K_nn)
    testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

    # Check sample
    set.seed(123)
    AD_sample <- rmvnorm0(n, mean0, sd0)@x
    FD_sample <- f(cbind(mu, sigma))
    testthat::expect_lt(sum(abs(AD_sample - FD_sample)), 1e-8)
  }
  # Check sample 2
  set.seed(123)
  s1 <- rmvnorm0(1, numeric(2), diag(2))
  set.seed(123)
  s2 <- rnorm(2)
  eq_transform(s1, s2, identity)
})
