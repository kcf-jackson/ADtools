testthat::context("Test subsetting a dual number")

# Functions without S3 declarations need to be wrapped
wrap <- function(fs) purrr::map(fs, ~lambda(.x(x)))

testthat::test_that("test subset, head, tail", {
  i <- sample(10, 1)
  j <- sample(10, 1)
  fs <- list(
    lambda(x, x[i,]),
    lambda(x, x[, j]),
    lambda(x, x[i, j]),
    head, tail
  )
  inputs <- generate_inputs(
    g = lambda( list(x = randn(ind[[1]], ind[[2]])) ),
    config_ls = purrr::map2(10:20, 10:20, list)
  )
  ctrl = list(display = F, err_fun = rel_err, epsilon = 1e-8)
  test_fs(wrap(fs), inputs, ctrl)
})
