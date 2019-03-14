testthat::context("Test matrix sum")

testthat::test_that("test sum, colSums, rowSums, tr", {
  fs <- list(sum, colSums, rowSums, tr)
  for (f in fs) {
    for (i in 2:10) {
      test_AD_with(
        f, input = list(x = randn(i, i)),
        ctrl = list(display = F, err_fun = rel_err, epsilon = 1e-6)
      )
    }
  }
})
