context("Test random variable simulation")


test_that("Test univariate normal simulation", {
  purrr::map(2:15, function(i) {
    # dual-dual case
    f <- function(mean, sd) {
      set.seed(123)
      rnorm0(n = i, mean = mean, sd = sd)
    }
    fs <- list(f)
    inputs <- list(
      list(mean = runif(i), sd = runif(i)),
      list(mean = runif(i), sd = runif(1)),  # different input length
      list(mean = runif(1), sd = runif(i)),  # different input length
      list(mean = runif(1), sd = runif(1))   # different input length
    )
    ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-6)
    test_fs(fs, inputs, ctrl)

    # dual-numeric case
    f2 <- function(mean) {
      set.seed(123)
      rnorm0(n = i, mean = mean, sd = i)
    }
    f2b <- function(mean) {
      set.seed(123)
      rnorm0(n = i, mean = mean, sd = runif(i))
    }
    fs <- list(f2, f2b)
    inputs <- list(
      list(mean = runif(i)),
      list(mean = runif(1))
    )
    ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-6)
    test_fs(fs, inputs, ctrl)

    # numeric-dual case
    f3 <- function(sd) {
      set.seed(123)
      rnorm0(n = i, mean = i, sd = sd)
    }
    f3b <- function(sd) {
      set.seed(123)
      rnorm0(n = i, mean = runif(i), sd = sd)
    }
    fs <- list(f3, f3b)
    inputs <- list(
      list(sd = runif(i)),
      list(sd = runif(1))
    )
    ctrl <- list(display = F, err_fun = abs_err, epsilon = 1e-6)
    test_fs(fs, inputs, ctrl)

    # Test output dimension is correct
    m <- 10
    k <- 2
    scalar <- 1
    vector <- 1:10
    scalar_dual <- dual(m, c(1, k), 1)       # dual-dim: k + 1
    vector_dual <- dual(seq(m), c(m, k), 1)  # dual-dim: k + m
    check_output_dim <- function(x, d) {
      testthat::expect_equal(dim(x@x), c(m, 1))
      testthat::expect_equal(dim(x@dx), c(m, d))
    }
    check_output_dim(rnorm0(m, mean = scalar_dual, sd = scalar), 1 + k)
    check_output_dim(rnorm0(m, mean = scalar_dual, sd = vector), 1 + k)
    check_output_dim(rnorm0(m, mean = vector_dual, sd = scalar), m + k)
    check_output_dim(rnorm0(m, mean = vector_dual, sd = vector), m + k)
    check_output_dim(rnorm0(m, sd = scalar_dual, mean = scalar), 1 + k)
    check_output_dim(rnorm0(m, sd = scalar_dual, mean = vector), 1 + k)
    check_output_dim(rnorm0(m, sd = vector_dual, mean = scalar), m + k)
    check_output_dim(rnorm0(m, sd = vector_dual, mean = vector), m + k)
  })
})


test_that("Test multivariate normal simulation", {
  # Test rmvnorm0 implementation matches with rmvnorm
  purrr::map(1:15, function(i) {
    set.seed(123)
    ref <- mvtnorm::rmvnorm(1, numeric(i), diag(i))
    set.seed(123)
    imp <- rmvnorm0(1, numeric(i), diag(i))
    compare(ref, imp, display = F, err_fun = rel_err, epsilon = 1e-8)
  })
  purrr::map(1:15, function(i) {
    f <- function(mean, sigma) {
      set.seed(123)
      rmvnorm0(n = i, mean = mean, sigma = sigma)
    }
    fs <- list(f)
    inputs <- list(list(mean = runif(i), sigma = diag(i)))
    ctrl <- list(display = F, err_fun = abs_err, epsilon = 1e-6)
    test_fs(fs, inputs, ctrl)
  })
})
