testthat::context("Test dual number arithmetic")

set.seed(123)
`%++%` <- append
f_plus <- function(e1, e2) { e1 + e2 }
f_minus <- function(e1, e2) { e1 - e2 }
f_times <- function(e1, e2) { e1 * e2 }
f_divide <- function(e1, e2) { e1 / e2 }
f_ktimes <- function(e1, e2) { e1 %x% e2 }
f_mtimes <- function(e1, e2) { e1 %*% e2 }


testthat::test_that("dual matrix AND dual matrix", {
  fs <- list(f_plus, f_minus, f_times, f_divide, f_ktimes)
  inputs <- generate_inputs(2:12,
    # add 10 to avoid dividing by near-zero
    lambda(list(e1 = randn(i, i), e2 = 10 + randn(i, i)))
  )
  ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-6)
  test_fs(fs, inputs, ctrl)

  # Matrix product is like taking squares; precision is often lost with
  # finite-difference method, so it is tested separately.
  fs <- list(f_mtimes)
  inputs <- generate_inputs(2:12,
    lambda(list(e1 = 1 + randn(i, i), e2 = 1 + randn(i, i)))
  )
  ctrl <- list(display = F, err_fun = abs_err, epsilon = 1e-6)
  test_fs(fs, inputs, ctrl)
})


testthat::test_that("2 x {dual scalar AND dual matrix} and {dual scalar AND dual scalar}", {
  fs <- list(f_plus, f_minus, f_times, f_divide)
  inputs <- generate_inputs(2:15, lambda(list(e1 = i, e2 = 10 + randn(i, i)))) %++%
    generate_inputs(2:15, lambda(list(e1 = 10 + randn(i, i), e2 = i))) %++%
    generate_inputs(2:15, lambda(list(e1 = i, e2 = i)))
  ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-6)
  test_fs(fs, inputs, ctrl)
})


testthat::test_that("numeric matrix AND dual matrix", {
  fs <- list(f_plus, f_minus, f_times, f_divide, f_ktimes)
  ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-6)
  purrr::map(2:15, function(i) {
    m0 <- randn(i, i)
    gs <- purrr::map(fs, function(f) lambda(e2, f(m0, e2)))
    inputs <- list( list(e2 = 10 + randn(i, i)) )
    test_fs(gs, inputs, ctrl)
  })
  # This needs to be done separately due to division.
  purrr::map(2:15, function(i) {
    m0 <- 10 + randn(i, i)
    gs <- purrr::map(fs, function(f) lambda(e1, f(e1, m0)))
    inputs <- list( list(e1 = randn(i, i)) )
    test_fs(gs, inputs, ctrl)
  })

  # Matrix multiplication
  fs <- list(f_mtimes)
  ctrl <- list(display = F, err_fun = abs_err, epsilon = 1e-6)
  purrr::map(2:15, function(i) {
    m0 <- randn(i, i)
    gs <- purrr::map(fs, function(f) lambda(e2, f(m0, e2)))
    inputs <- list( list(e2 = randn(i, i)) )
    test_fs(gs, inputs, ctrl)
  })
  purrr::map(2:15, function(i) {
    m0 <- randn(i, i)
    gs <- purrr::map(fs, function(f) lambda(e1, f(e1, m0)))
    inputs <- list( list(e1 = randn(i, i)) )
    test_fs(gs, inputs, ctrl)
  })
})


testthat::test_that("numeric scalar AND dual matrix", {
  fs <- list(f_plus, f_minus, f_times, f_divide)
  ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-6)
  purrr::map(2:15, function(i) {
    m0 <- i
    gs <- purrr::map(fs, function(f) lambda(e2, f(m0, e2)))
    inputs <- list( list(e2 = 10 + randn(i, i)) )
    test_fs(gs, inputs, ctrl)
  })
  purrr::map(2:15, function(i) {
    m0 <- 10 + i
    gs <- purrr::map(fs, function(f) lambda(e1, f(e1, m0)))
    inputs <- list( list(e1 = randn(i, i)) )
    test_fs(gs, inputs, ctrl)
  })
})


testthat::test_that("numeric matrix AND dual scalar", {
  fs <- list(f_plus, f_minus, f_times, f_divide)
  ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-6)
  purrr::map(2:15, function(i) {
    m0 <- randn(i, i)
    gs <- purrr::map(fs, function(f) lambda(e2, f(m0, e2)))
    inputs <- list( list(e2 = 10 + runif(1, max = i)) )
    test_fs(gs, inputs, ctrl)
  })
  purrr::map(2:15, function(i) {
    m0 <- 10 + randn(i, i)
    gs <- purrr::map(fs, function(f) lambda(e1, f(e1, m0)))
    inputs <- list( list(e1 = runif(1, max = i)) )
    test_fs(gs, inputs, ctrl)
  })
})


# Power (careful with low precision for high exponents)
testthat::test_that("Taking element-wise power", {
  ctrl <- list(display = F, err_fun = abs_err, epsilon = 1e-6)
  purrr::map(2:10, function(i) {
    # dual matrix - dual scalar
    # base must have only positive entries
    gs <- list( lambda(e1 ^ e2) )
    inputs <- list( list(e1 = randu(i, i, min = 0.1), e2 = i) )
    test_fs(gs, inputs, ctrl)

    # # numeric matrix - dual scalar
    m0 <- randu(i, i, min = 0.1)   # base must have only positive entries
    gs <- list( lambda(e2, m0 ^ e2) )
    inputs <- list( list(e2 = i) )
    test_fs(gs, inputs, ctrl)

    # dual matrix - numeric scalar
    gs <- list(lambda(e1, e1 ^ i), lambda(e1, e1 ^ 0))
    inputs <- list( list(e1 = randu(i, i, min = 0.1)) )
    test_fs(gs, inputs, ctrl)

    # Error case
    testthat::expect_error((-3)^dual(3, c(1, 2, 3), 1))                      # must have positive scalar base
    testthat::expect_error(3^dual(1:3, c(1, 2, 3), 1))                       # must have scalar exponent
    testthat::expect_error(dual(-3, c(1, 2, 3), 1)^(1:3))                    # must have scalar exponent
    testthat::expect_error(dual(-3, c(1, 2, 3), 1)^dual(3, c(1, 2, 3), 1))   # must have positive scalar base
    testthat::expect_error(dual(3, c(1, 2, 3), 1)^dual(1:3, c(1, 2, 3), 1))  # must have scalar exponent
  })
})


# Additional tests for cross-class arithmetic
testthat::test_that("Testing cross-class arithmetic", {
  ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-6)

  f <- function(x) x + 1:2
  test_fs(list(f), list(list(x = 5), list(x = 2:3)), ctrl)

  f2 <- function(x) x + matrix(1:4, 2, 2)
  f3 <- function(x) diag(x) + matrix(1:4, 2, 2)
  test_fs(list(f2, f3), list(list(x = matrix(1:4, 2, 2))), ctrl)
})
