testthat::context("Test element-wise derivative")

testthat::test_that("test sqrt, sin, cos, tan, exp, log", {
  fs <- list(sqrt, sin, cos, tan, exp, log)
  for (f in fs) {
    for (i in 2:10) {
      test_AD_with(
        f, input = list(x = 100 + randn(i, i)),
        ctrl = list(display = F, err_fun = rel_err, epsilon = 1e-5)
      )
    }
  }
})
