testthat::context("Test reshaping dual numbers")

# Functions without S3 declarations need to be wrapped
wrap <- function(fs) purrr::map(fs, ~lambda(.x(x)))


testthat::test_that("diag, vec, vech, as.vector, as.matrix, matrix", {
  g <- function(x) diag(diag(x))
  g2 <- function(x) matrix(x)
  g3 <- function(x) matrix(x, nrow = 2)
  g4 <- function(x) matrix(x, nrow = 2, byrow = T)

  fs <- list(diag, vec, vech, g, as.vector, as.matrix, g, g2, g3, g4)
  inputs <- generate_inputs(seq(2, 20, 2), lambda( list(x = randn(i,i)) ))
  ctrl = list(display = F, err_fun = rel_err, epsilon = 1e-8)
  test_fs(wrap(fs), inputs, ctrl)

  fs <- list(diag)
  inputs <- generate_inputs(1:10, lambda( list(x = i) ))
  ctrl = list(display = F, err_fun = rel_err, epsilon = 1e-8)
  test_fs(wrap(fs), inputs, ctrl)

  # Expect error
  testthat::expect_error(diag(dual(array(1:10), c(10, 3), 1)))
  testthat::expect_error(vech(randn(5, 6)))
})


testthat::test_that("Check diag has correct behaviour for different classes of input", {
  k <- 2
  scalar0 <- 4
  vector0 <- 1:5
  matrix0 <- matrix(1:9, 3, 3)
  scalar_dual <- dual(scalar0, c(length(scalar0), k), 1)
  vector_dual <- dual(vector0, c(length(vector0), k), 1)
  matrix_dual <- dual(matrix0, c(length(matrix0), k), 1)
  equal_in_value <- function(x, y) {
    testthat::expect_equal(sum(abs(x - y)), 0)
  }
  equal_in_value(diag(scalar_dual)@x, diag(scalar0))
  equal_in_value(diag(vector_dual)@x, diag(vector0))
  equal_in_value(diag(matrix_dual)@x, diag(matrix0))

  testthat::expect_true(check_dual(diag(vector_dual)))
  testthat::expect_true(check_dual(diag(matrix_dual)))
})
