context("Test random variable simulation")

# Helper functions
relative_diff <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  abs(x - y) / max(abs(x), abs(y))
}

compare_FD_and_AD <- function(FD, AD, show = F) {
  AD <- AD[seq(nrow(FD)), seq(ncol(FD))]
  rel_err <- relative_diff(FD, AD)
  if (show) {
    print(cbind(FD = as.numeric(FD), AD = as.numeric(AD)))
    print(glue::glue("Maximum relative error over all entries: {max(rel_err)} (entry: {which.max(rel_err)})"))
  }
  max(rel_err)
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
      rnorm(n, mean = param[1], sd = param[2])
    }
    FD_res <- finite_diff(f, c(mu, sigma))
    testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)
    # Check sample
    set.seed(123)
    AD_sample <- rnorm0(n, mean = mean0, sd = sd0)@x
    FD_sample <- f(c(mu, sigma))
    testthat::expect_lt(sum(abs(AD_sample - FD_sample)), 1e-8)
  }
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
    FD_res <- finite_diff(f, cbind(mu, sigma))
    FD_res[, -seq_along(mu)] <- as.matrix(FD_res[, -seq_along(mu)] %*% K_nn)
    testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

    # Check sample
    set.seed(123)
    AD_sample <- rmvnorm0(n, mean0, sd0)@x
    FD_sample <- f(cbind(mu, sigma))
    testthat::expect_lt(sum(abs(AD_sample - FD_sample)), 1e-8)
  }
})

test_that("Test exponential simulation", {
  lambda <- 3
  # Check derivative
  for (n in 1:10) {
    # AD
    rate0 <- dual(lambda, list(A = 1, B = 4), 1)
    set.seed(123)
    AD_res <- rexp0(n, rate0)@dx
    # FD
    f <- function(param) {
      set.seed(123)
      rexp0(n, rate = param)
    }
    FD_res <- finite_diff(f, lambda)
    testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

    # Check sample
    set.seed(123)
    AD_sample <- rexp0(n, rate0)@x
    FD_sample <- f(lambda)
    testthat::expect_lt(sum(abs(AD_sample - FD_sample)), 1e-8)
  }
})

test_that("Test inverse transform gamma", {
  set.seed(123)
  s1 <- rgamma0(10000, shape = 3, scale = 5, method = "base")
  s2 <- rgamma0(10000, shape = 3, scale = 5, method = "inv_tf")
  expect_lt(relative_diff(mean(s1), mean(s2)), 0.1)
  expect_lt(relative_diff(sd(s1), sd(s2)), 0.1)
})

test_that("Test gamma simulation", {
  for (n in 1:10) {
    for (i in 1:10) {
      for (j in 1:10) {
        set.seed(10 * i + j)
        shape <- as.numeric(sample(20, 1))
        scale <- as.numeric(sample(20, 1))
        # ==== Check derivative dual / dual ====
        # AD
        param <- list(shape = 1, scale = 1, C = 2)
        shape0 <- dual(shape, param, 1)
        scale0 <- dual(scale, param, 2)
        set.seed(123)
        AD_res <- rgamma0(n, shape0, scale0, "inv_tf")@dx
        # FD
        f <- function(param) {
          set.seed(123)
          rgamma0(n, param[1], param[2], "inv_tf")
        }
        FD_res <- finite_diff(f, c(shape, scale), h = 1e-7)
        testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

        # ==== Check derivative dual / num ====
        # AD
        param <- list(shape = 1, scale = 1, C = 2)
        shape0 <- dual(shape, param, 1)
        set.seed(123)
        AD_res <- rgamma0(n, shape0, scale, "inv_tf")@dx
        # FD
        f <- function(shape) {
          set.seed(123)
          rgamma0(n, shape, scale, "inv_tf")
        }
        FD_res <- finite_diff(f, shape, h = 1e-6)
        testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

        # ==== Check derivative num / dual ====
        # AD
        param <- list(shape = 1, scale = 1, C = 2)
        scale0 <- dual(scale, param, 2)
        set.seed(123)
        AD_res <- rgamma0(n, shape, scale0, "inv_tf")@dx[, 2, drop = F]
        # FD
        f <- function(scale) {
          set.seed(123)
          rgamma0(n, shape, scale, "inv_tf")
        }
        FD_res <- finite_diff(f, scale, h = 1e-7)
        testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

        # ==== Check sample ====
        f <- function(param) {
          set.seed(123)
          rgamma0(n, param[1], param[2])
        }
        set.seed(123)
        AD_sample <- rgamma0(n, shape0, scale0)@x
        FD_sample <- f(c(shape, scale))
        testthat::expect_lt(sum(abs(AD_sample - FD_sample)), 1e-6)
      }
    }
  }
})

testthat::test_that("Test Chi-squared simulation", {
  df <- 5
  for (n in 1:10) {
    # AD
    df0 <- dual(df, list(A = 1, B = 4), 1)
    set.seed(123)
    AD_res <- rchisq0(n, df0, method = "inv_tf")@dx

    # FD
    f <- function(param) {
      set.seed(123)
      rchisq0(n, param, method = "inv_tf")
    }
    FD_res <- finite_diff(f, df)
    testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

    # Check sample
    set.seed(123)
    AD_sample <- rchisq0(n, df0, method = "inv_tf")@x
    FD_sample <- f(df)
    testthat::expect_lt(sum(abs(AD_sample - FD_sample)), 1e-8)
  }
})

testthat::test_that("Wishart simulation", {
  v <- 5
  M <- diag(3)

  v_dual <- dual(5, list(v = 1, other = 9), 1)
  M_dual <- dual(diag(3), list(v = 1, other = 9), 2)

  f_v <- function(param) {
    set.seed(123)
    rWishart0(param, M, method = "inv_tf")
  }
  set.seed(123)
  AD_res <- rWishart0(v_dual, M, method = "inv_tf")@dx
  FD_res <- finite_diff(f_v, v)
  expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

  f_M <- function(param) {
    set.seed(123)
    rWishart0(v, param, method = "inv_tf")
  }
  set.seed(123)
  AD_res <- get_deriv(rWishart0(v, M_dual, method = "inv_tf"), "other")
  FD_res <- finite_diff(f_M, M)
  # An extra commutation matrix is needed because rWishart0 uses Cholesky
  # decomposition; see test_matrix_calculus.R: `d_chol` for more detail.
  K_nn <- commutation_matrix0(nrow(M), ncol(M))
  FD_res <- FD_res %*% K_nn
  expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

  f <- function(param) {
    set.seed(123)
    v <- param[1]
    M <- matrix(param[-1], nrow = sqrt(length(param[-1])))
    rWishart0(v, M, method = "inv_tf")
  }
  set.seed(123)
  AD_res <- rWishart0(v_dual, M_dual, method = "inv_tf")@dx
  FD_res <- finite_diff(f, c(v, M))
  # An extra commutation matrix is needed because rWishart0 uses Cholesky
  # decomposition; see test_matrix_calculus.R: `d_chol` for more detail.
  K_nn <- commutation_matrix0(nrow(M), ncol(M))
  FD_res[,-1] <- as.matrix(FD_res[,-1] %*% K_nn)
  expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)
})
