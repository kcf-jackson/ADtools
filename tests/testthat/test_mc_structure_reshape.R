testthat::context("Test reshaping dual numbers")

testthat::test_that("diag, vec, vech, as.vector, as.matrix, matrix", {
  g <- function(x) diag(diag(x))
  g2 <- function(x) matrix(x)
  g3 <- function(x) matrix(x, nrow = 2)
  g4 <- function(x) matrix(x, nrow = 2, byrow = T)

  fs <- list(diag, vec, vech, g, as.vector, as.matrix, g, g2, g3, g4)
  inputs <- generate_inputs(seq(2, 20, 2), lambda( list(x = randn(i,i)) ))
  ctrl = list(display = F, err_fun = rel_err, epsilon = 1e-8)
  test_fs(fs, inputs, ctrl)

  fs <- list(diag)
  inputs <- generate_inputs(1:10, lambda( list(x = i) ))
  ctrl = list(display = F, err_fun = rel_err, epsilon = 1e-8)
  test_fs(fs, inputs, ctrl)

  # Expect error
  testthat::expect_error(diag(dual(array(1:10), c(10, 3), 1)))
  testthat::expect_error(vech(randn(5, 6)))
})
