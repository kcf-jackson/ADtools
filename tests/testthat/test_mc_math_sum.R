testthat::context("Test matrix sum")

testthat::test_that("test sum, colSums, rowSums, tr", {
  fs <- list(sum, colSums, rowSums, tr)
  inputs <- generate_inputs(2:10, lambda(list(x = randn(i, i))))
  test_fs(fs, inputs, list(display = F, err_fun = rel_err, epsilon = 1e-6))
})
