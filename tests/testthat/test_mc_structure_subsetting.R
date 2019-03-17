testthat::context("Test subsetting a dual number")


i <- sample(10, 1)
j <- sample(10, 1)
fs <- list(
  lambda(x, x[i,]), 
  lambda(x, x[, j]), 
  lambda(x, x[i, j]), 
  head, tail
)
inputs <- generate_input(
  g = lambda( list(x = randn(ind[[1]], ind[[2]])) ),
  config_ls = purrr::map2(10:20, 10:20, list)
)
ctrl = list(display = F, err_fun = rel_err, epsilon = 1e-7)
testthat::test_that("test subset, head, tail", {
  test_fs(fs, inputs, ctrl)
})
