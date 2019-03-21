context("Test random variable simulation")

test_that("Test Exponential simulation", {
  purrr::map(2:10, function(i) {
    f <- function(rate) {
      set.seed(123)
      rexp0(n = i, rate = rate)
    }
    fs <- list(f)
    inputs <- list(
      list(rate = runif(i)),
      list(rate = runif(1))
    )
    ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-6)
    test_fs(fs, inputs, ctrl)
  })
})


test_that("Test Chi-squared simulation", {
  purrr::map(2:10, function(i) {
    f <- function(df) {
      set.seed(123)
      rchisq0(n = i, df = df)
    }
    fs <- list(f)
    inputs <- list(
      list(df = 1 + sample(30, i)),
      list(df = 1 + sample(30, 1))
    )
    ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-5)
    test_fs(fs, inputs, ctrl)
  })
})


test_that("Test Wishart simulation", {
  purrr::map(2:10, function(i) {
    # dual-dual case
    f <- function(v, M) {
      set.seed(123)
      rWishart0(v = v, M = M)
    }
    fs <- list(f)
    inputs <- list(list(v = i, M = crossprod(randn(i, i))))
    ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-5)
    test_fs(fs, inputs, ctrl)


    # dual-numeric case
    f2 <- function(M) {
      set.seed(123)
      rWishart0(v = i + 5, M = M)
    }
    fs <- list(f2)
    inputs <- list(list(M = crossprod(randn(i, i))))
    ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-6)
    test_fs(fs, inputs, ctrl)


    # dual-numeric case
    m0 <- crossprod(randn(i, i))
    f3 <- function(v) {
      set.seed(123)
      rWishart0(v = v, M = m0)
    }
    fs <- list(f3)
    inputs <- list(list(v = i + 5))
    ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-6)
    test_fs(fs, inputs, ctrl)
  })
})


test_that("Test Gamma simulation", {
  purrr::map(2:10, function(i) {
    # dual-dual case
    f <- function(shape, scale) {
      set.seed(123)
      rgamma0(n = i, shape = shape, scale = scale)
    }
    inputs <- list(
      list(shape = 3 + runif(i), scale = 3 + runif(i)),
      list(shape = 3 + runif(1), scale = 3 + runif(i)),
      list(shape = 3 + runif(i), scale = 3 + runif(1)),
      list(shape = 3 + runif(1), scale = 3 + runif(1))
    )
    ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-6)
    test_fs(list(f), inputs, ctrl)

    # dual-numeric case
    f2 <- function(shape) {
      set.seed(123)
      rgamma0(n = i, shape = shape, scale = 1)
    }
    inputs <- list(
      list(shape = 3 + runif(i)),
      list(shape = 3 + runif(1))
    )
    ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-6)
    test_fs(list(f2), inputs, ctrl)

    # dual-numeric case
    f3 <- function(scale) {
      set.seed(123)
      rgamma0(n = i, shape = 1, scale = scale)
    }
    inputs <- list(
      list(scale = runif(i)),
      list(scale = runif(1))
    )
    ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-6)
    test_fs(list(f3), inputs, ctrl)
  })
})
