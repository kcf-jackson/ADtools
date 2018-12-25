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

# Tests
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
    FD_res <- finite_diff(f, c(mu, sigma))
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
    FD_res <- finite_diff(f, cbind(mu, sigma))
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
  # Test cases where rate has length > 1
  param_dim <- list(A = 1, B = 2, C = 1)
  s <- runif(10, 1, 10)
  rate0 <- dual(s, param_dim, -1)
  rate0@dx[,1] <- 1

  set.seed(123)
  AD_res <- rexp0(10, rate0)
  set.seed(123)
  AD_res_2 <- mapreduce(1:10, ~rexp0(1, rate0[.x]), rbind2)
  eq_transform(AD_res, AD_res_2, parent_of, F)
  aeq_transform(AD_res, AD_res_2, deriv_of, F)

  testthat::expect_error(rexp0(9, rate0))
})

test_that("Test gamma simulation with inverse-transform method", {
  testthat::expect_error(rgamma0(10, shape = 3, scale = 5, method = "a"))

  set.seed(123)
  s1 <- rgamma0(10000, shape = 3, scale = 5, method = "base")
  s2 <- rgamma0(10000, shape = 3, scale = 5, method = "inv_tf")
  expect_lt(relative_diff(mean(s1), mean(s2)), 0.1)
  expect_lt(relative_diff(sd(s1), sd(s2)), 0.1)

  shape <- c(3, 1)
  s1 <- rgamma0(10000, shape = shape, scale = 5, method = "base")
  s2 <- rgamma0(10000, shape = shape, scale = 5, method = "inv_tf")
  odd <- seq(1, 10000, 2)
  expect_lt(relative_diff(mean(s1[odd]), mean(s2[odd])), 0.1)
  expect_lt(relative_diff(sd(s1[odd]), sd(s2[odd])), 0.1)
  even <- seq(2, 10000, 2)
  expect_lt(relative_diff(mean(s1[even]), mean(s2[even])), 0.1)
  expect_lt(relative_diff(sd(s1[even]), sd(s2[even])), 0.1)

  scale <- c(2, 7)
  s1 <- rgamma0(10000, shape = 3, scale = scale, method = "base")
  s2 <- rgamma0(10000, shape = 3, scale = scale, method = "inv_tf")
  odd <- seq(1, 10000, 2)
  expect_lt(relative_diff(mean(s1[odd]), mean(s2[odd])), 0.1)
  expect_lt(relative_diff(sd(s1[odd]), sd(s2[odd])), 0.1)
  even <- seq(2, 10000, 2)
  expect_lt(relative_diff(mean(s1[even]), mean(s2[even])), 0.1)
  expect_lt(relative_diff(sd(s1[even]), sd(s2[even])), 0.1)

  shape <- c(3, 8)
  scale <- c(5, 9)
  s1 <- rgamma0(10000, shape = shape, scale = scale, method = "base")
  s2 <- rgamma0(10000, shape = shape, scale = scale, method = "inv_tf")
  odd <- seq(1, 10000, 2)
  expect_lt(relative_diff(mean(s1[odd]), mean(s2[odd])), 0.1)
  expect_lt(relative_diff(sd(s1[odd]), sd(s2[odd])), 0.1)
  even <- seq(2, 10000, 2)
  expect_lt(relative_diff(mean(s1[even]), mean(s2[even])), 0.1)
  expect_lt(relative_diff(sd(s1[even]), sd(s2[even])), 0.1)

  shape <- c(3, 8, 7)
  scale <- 5
  s1 <- rgamma0(10000, shape = shape, scale = scale, method = "base")
  s2 <- rgamma0(10000, shape = shape, scale = scale, method = "inv_tf")
  for (i in 1:3) {
    s_seq <- seq(i, 10000, 3)
    expect_lt(relative_diff(mean(s1[s_seq]), mean(s2[s_seq])), 0.1)
    expect_lt(relative_diff(sd(s1[s_seq]), sd(s2[s_seq])), 0.1)
  }
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
  # Check unequal length gamma simulation
  shape <- runif(10, 1, 10)
  scale <- runif(10, 1, 10)
  param <- list(shape = 1, scale = 1, C = 2)
  shape_dual <- dual(shape, param, -1)
  shape_dual@dx[,1] <- 1
  scale_dual <- dual(scale, param, -1)
  shape_dual@dx[,2] <- 1
  unit_shape_dual <- dual(5, param, 1)
  unit_scale_dual <- dual(3, param, 2)

  # dual / dual
  set.seed(123)
  AD_res <- rgamma0(10, shape_dual, scale_dual)
  set.seed(123)
  AD_res_2 <- mapreduce(1:10, ~rgamma0(1, shape_dual[.x], scale_dual[.x]), rbind2)
  eq_transform(AD_res, AD_res_2, parent_of)
  aeq_transform(AD_res, AD_res_2, deriv_of)

  set.seed(123)
  AD_res <- rgamma0(10, unit_shape_dual, scale_dual)
  set.seed(123)
  AD_res_2 <- mapreduce(1:10, ~rgamma0(1, unit_shape_dual, scale_dual[.x]), rbind2)
  eq_transform(AD_res, AD_res_2, parent_of)
  aeq_transform(AD_res, AD_res_2, deriv_of)

  set.seed(123)
  AD_res <- rgamma0(10, shape_dual, unit_scale_dual)
  set.seed(123)
  AD_res_2 <- mapreduce(1:10, ~rgamma0(1, shape_dual[.x], unit_scale_dual), rbind2)
  eq_transform(AD_res, AD_res_2, parent_of)
  aeq_transform(AD_res, AD_res_2, deriv_of)

  set.seed(123)
  AD_res <- rgamma0(10, unit_shape_dual, unit_scale_dual)
  set.seed(123)
  AD_res_2 <- mapreduce(1:10, ~rgamma0(1, unit_shape_dual, unit_scale_dual), rbind2)
  eq_transform(AD_res, AD_res_2, parent_of)
  aeq_transform(AD_res, AD_res_2, deriv_of)

  # num / dual
  set.seed(123)
  AD_res <- rgamma0(10, 5, scale_dual)
  set.seed(123)
  AD_res_2 <- mapreduce(1:10, ~rgamma0(1, 5, scale_dual[.x]), rbind2)
  identical(AD_res, AD_res_2)

  shape_vec <- runif(10, 1, 10)
  set.seed(123)
  AD_res <- rgamma0(10, shape_vec, scale_dual)
  set.seed(123)
  AD_res_2 <- mapreduce(1:10, ~rgamma0(1, shape_vec[.x], scale_dual[.x]), rbind2)
  identical(AD_res, AD_res_2)

  set.seed(123)
  AD_res <- rgamma0(10, shape_vec, unit_scale_dual)
  set.seed(123)
  AD_res_2 <- mapreduce(1:10, ~rgamma0(1, shape_vec[.x], unit_scale_dual), rbind2)
  identical(AD_res, AD_res_2)

  # dual / num
  set.seed(123)
  AD_res <- rgamma0(10, shape_dual, 4)
  set.seed(123)
  AD_res_2 <- mapreduce(1:10, ~rgamma0(1, shape_dual[.x], 4), rbind2)
  eq_transform(AD_res, AD_res_2, parent_of)
  aeq_transform(AD_res, AD_res_2, deriv_of)

  # num / num
  set.seed(123)
  AD_res <- rgamma0(10, 4, 5)
  set.seed(123)
  AD_res_2 <- mapreduce(1:10, ~rgamma0(1, 4, 5), c)
  eq_transform(AD_res, AD_res_2, identity)
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
  # Check consistency between base R and inverse transform
  s1 <- rchisq0(1e4, df, method = "inv_tf")
  s2 <- rchisq0(1e4, df, method = "base")
  testthat::expect_lt(abs(mean(s1) - mean(s2)), 0.1)
  testthat::expect_lt(abs(sd(s1) - sd(s2)), 0.1)
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
  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

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
  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

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
  testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)

  # Check consistency with base R implementation
  testthat::expect_identical(
    dim(rWishart0(5, diag(2), method = "inv_tf")),
    dim(rWishart0(5, diag(2)))
  )
  set.seed(123)
  inv_tf_wishart <- 1:2000 %>%
    purrr::map(~as.numeric(rWishart0(5, diag(2), method = "inv_tf"))) %>%
    do.call(rbind, .) %>%
    apply(2, mean)
  base_wishart <- 1:2000 %>%
    purrr::map(~as.numeric(rWishart0(5, diag(2)))) %>%
    do.call(rbind, .) %>%
    apply(2, mean)
  testthat::expect_lt(max(abs(inv_tf_wishart - base_wishart)), 0.2)
})
