testthat::context("Test element-wise derivative")

# Functions without S3 declarations need to be wrapped
wrap <- function(fs) purrr::map(fs, ~lambda(.x(x)))

testthat::test_that("test sqrt, sin, cos, tan, exp, log", {
  ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-6)
  fs <- list(sin, cos, exp)  # exp requires smallish numbers
  inputs <- generate_inputs(2:15, lambda(list(x = randu(i, i, min = 1))))
  test_fs(wrap(fs), inputs, ctrl)

  ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-6)
  fs <- list(sqrt, log)  # requires positive numbers
  inputs <- generate_inputs(2:15, lambda(list(x = 20 + randn(i, i))))
  test_fs(wrap(fs), inputs, ctrl)

  fs <- list(tan)  # explodes at pi / 2, so needs extra care
  inputs <- generate_inputs(2:15, lambda(list(x = matrix(runif(i^2), i, i))))
  test_fs(wrap(fs), inputs, ctrl)
})
