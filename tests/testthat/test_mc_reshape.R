testthat::context("Test reshaping dual numbers")

testthat::test_that("diag, vec, vech", {
  g <- function(x) diag(diag(x))
  fs <- list(diag, vec, vech, g)
  for (f in fs) {
    for (i in 2:10) {
      test_AD_with(
        f, input = list(x = randn(i, i)),
        ctrl = list(display = F, err_fun = rel_err, epsilon = 1e-6)
      )
    }
  }
})
