testthat::context("Test reshaping dual numbers")

testthat::test_that("diag, vec, vech, as.vector, as.matrix, matrix", {
  g <- function(x) diag(diag(x))
  g2 <- function(x) matrix(x)
  fs <- list(diag, vec, vech, g, as.vector, as.matrix, g2)
  inputs <- generate_inputs(2:10, lambda( list(x = randn(i,i))) )
  ctrl = list(display = F, err_fun = rel_err, epsilon = 1e-8)
  test_fs(fs, inputs, ctrl)
})
