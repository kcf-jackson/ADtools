testthat::context("Test matrix_calculus")
set.seed(123)

push <- function(res, el) {
  res[[length(res) + 1]] <- el
  res
}

perturb_vector <- function(v, h) {
  seq_along(v) %>%
    purrr::map(function(i) {
      u <- v
      u[i] <- u[i] + h
      u
    }) %>%
    list(origin = v, configs = ., h = h)
}

perturb_matrix <- function(A, h) {
  1:length(A) %>%
    purrr::map(function(i) {
      B <- as.numeric(A)
      B[i] <- B[i] + h
      matrix(B, nrow = nrow(A))
    }) %>%
    list(origin = A, configs = ., h = h)
}

finite_diff <- function(param, fun) {
  diff_fun <- . %>% { (fun(.) - fun(param$origin)) / param$h }
  purrr::map(param$configs, diff_fun) %>%
    purrr::map(~.x %>% as.numeric()) %>%
    do.call(cbind, .)
}


testthat::test_that("Test d_sum against finite difference", {
  A <- randn(3,3)
  B <- randn(3,3)
  h <- 1e-8
  f <- function(x) { x + B }
  f2 <- function(x) { B + x }
  finite_res <- finite_diff(param = perturb_matrix(A, h), fun = f)
  finite_res2 <- finite_diff(param = perturb_matrix(A, h), fun = f2)

  len_A <- length(A)
  len_B <- length(B)
  X <- new("dual", x = A, dx = init_dx(list(c(len_A, len_A)), 1))
  Y <- new("dual", x = B, dx = init_dx(list(c(len_B, len_B)), -1))
  auto_res <- element_of(deriv_of(X + Y))[[1]]
  auto_res2 <- element_of(deriv_of(Y + X))[[1]]

  testthat::expect_lt(
    max(abs(as.numeric(finite_res) - as.numeric(auto_res))),
    1e-6
  )

  testthat::expect_lt(
    max(abs(as.numeric(finite_res2) - as.numeric(auto_res2))),
    1e-6
  )
})

testthat::test_that("Test d_minus against finite difference", {
  A <- randn(3,3)
  B <- randn(3,3)
  h <- 1e-8
  f <- function(x) { x - B }
  f2 <- function(x) { B - x }
  finite_res <- finite_diff(param = perturb_matrix(A, h), fun = f)
  finite_res2 <- finite_diff(param = perturb_matrix(A, h), fun = f2)

  len_A <- length(A)
  len_B <- length(B)
  X <- new("dual", x = A, dx = init_dx(list(c(len_A, len_A)), 1))
  Y <- new("dual", x = B, dx = init_dx(list(c(len_B, len_B)), -1))
  auto_res <- element_of(deriv_of(X - Y))[[1]]
  auto_res2 <- element_of(deriv_of(Y - X))[[1]]

  testthat::expect_lt(
    max(abs(as.numeric(finite_res) - as.numeric(auto_res))),
    1e-6
  )

  testthat::expect_lt(
    max(abs(as.numeric(finite_res2) - as.numeric(auto_res2))),
    1e-6
  )
})

testthat::test_that("Test d_matrix_prod against finite difference", {
  A <- randn(3,3)
  B <- randn(3,3)
  h <- 1e-8
  f <- function(x) { x %*% B }
  f2 <- function(x) { B %*% x }
  finite_res <- finite_diff(param = perturb_matrix(A, h), fun = f)
  finite_res2 <- finite_diff(param = perturb_matrix(A, h), fun = f2)

  len_A <- length(A)
  len_B <- length(B)
  X <- new("dual", x = A, dx = init_dx(list(c(len_A, len_A)), 1))
  Y <- new("dual", x = B, dx = init_dx(list(c(len_B, len_B)), -1))
  auto_res <- element_of(deriv_of(X %*% Y))[[1]]
  auto_res2 <- element_of(deriv_of(Y %*% X))[[1]]

  testthat::expect_lt(
    max(abs(as.numeric(finite_res) - as.numeric(auto_res))),
    1e-6
  )

  testthat::expect_lt(
    max(abs(as.numeric(finite_res2) - as.numeric(auto_res2))),
    1e-6
  )
})

testthat::test_that("Test d_kronecker against finite difference", {
  A <- randn(3,3)
  B <- randn(3,3)
  h <- 1e-8
  f <- function(x) { x %x% B }
  f2 <- function(x) { B %x% x }
  finite_res <- finite_diff(param = perturb_matrix(A, h), fun = f)
  finite_res2 <- finite_diff(param = perturb_matrix(A, h), fun = f2)

  len_A <- length(A)
  len_B <- length(B)
  X <- new("dual", x = A, dx = init_dx(list(c(len_A, len_A)), 1))
  Y <- new("dual", x = B, dx = init_dx(list(c(len_B, len_B)), -1))
  auto_res <- element_of(deriv_of(X %x% Y))[[1]]
  auto_res2 <- element_of(deriv_of(Y %x% X))[[1]]

  testthat::expect_lt(
    max(abs(as.numeric(finite_res) - as.numeric(auto_res))),
    1e-6
  )

  testthat::expect_lt(
    max(abs(as.numeric(finite_res2) - as.numeric(auto_res2))),
    1e-6
  )
})

testthat::test_that("Test d_solve against finite difference", {
  A <- randn(3,3) * 10
  h <- 1e-8
  f <- solve
  finite_res <- finite_diff(param = perturb_matrix(A, h), fun = f)

  len_A <- length(A)
  X <- new("dual", x = A, dx = init_dx(list(c(len_A, len_A)), 1))
  auto_res <- element_of(deriv_of(f(X)))[[1]]

  testthat::expect_lt(
    max(abs(as.numeric(finite_res) - as.numeric(auto_res))),
    1e-6
  )
})

testthat::test_that("Test d_transpose against finite difference", {
  A <- randn(3,3)
  h <- 1e-8
  f <- t
  finite_res <- finite_diff(param = perturb_matrix(A, h), fun = f)

  len_A <- length(A)
  X <- new("dual", x = A, dx = init_dx(list(c(len_A, len_A)), 1))
  auto_res <- element_of(deriv_of(f(X)))[[1]]

  testthat::expect_lt(
    max(abs(as.numeric(finite_res) - as.numeric(auto_res))),
    1e-6
  )
})

testthat::test_that("Test d_XXT against finite difference", {
  A <- randn(3,3)
  h <- 1e-8
  f <- tcrossprod
  finite_res <- finite_diff(param = perturb_matrix(A, h), fun = f)

  len_A <- length(A)
  X <- new("dual", x = A, dx = init_dx(list(c(len_A, len_A)), 1))
  auto_res <- element_of(deriv_of(f(X)))[[1]]

  testthat::expect_lt(
    max(abs(as.numeric(finite_res) - as.numeric(auto_res))),
    1e-6
  )
})

testthat::test_that("Test crossprod against finite difference", {
  A <- randn(3,3)
  h <- 1e-8
  f <- crossprod
  finite_res <- finite_diff(param = perturb_matrix(A, h), fun = f)

  len_A <- length(A)
  X <- new("dual", x = A, dx = init_dx(list(c(len_A, len_A)), 1))
  auto_res <- element_of(deriv_of(f(X)))[[1]]

  testthat::expect_lt(
    max(abs(as.numeric(finite_res) - as.numeric(auto_res))),
    1e-6
  )
})


testthat::test_that("Test d_chol against finite difference", {
  A <- randn(3,3)
  A <- A %*% t(A)
  h <- 1e-8
  f <- chol0
  finite_res <- finite_diff(param = perturb_matrix(A, h), fun = f)
  finite_res <- finite_res %*% commutation_matrix(nrow(A), ncol(A))

  len_A <- length(A)
  X <- new("dual", x = A, dx = init_dx(list(c(len_A, len_A)), 1))
  auto_res <- element_of(deriv_of(chol0(X)))[[1]]

  testthat::expect_lt(
    max(abs(as.numeric(finite_res) - as.numeric(auto_res))),
    1e-6
  )
})
