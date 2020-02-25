testthat::context("Test matrix sum")

# Functions without S3 declarations need to be wrapped
wrap <- function(fs) purrr::map(fs, ~lambda(.x(x)))

testthat::test_that("test sum, colSums, rowSums, tr, mean, colMeans, rowMeans", {
  fs <- list(sum, colSums, rowSums, tr, mean, colMeans, rowMeans)
  inputs <- generate_inputs(2:10, lambda(list(x = randn(i, i))))
  test_fs(wrap(fs), inputs, list(display = F, err_fun = rel_err, epsilon = 1e-6))

  testthat::expect_error(tr(randn(10, 11)))  # must be square
})
