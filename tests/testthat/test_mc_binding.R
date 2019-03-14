testthat::context("Test binding two dual numbers")

testthat::test_that("test cbind2, rbind2", {
  g <- function(x, y) rbind2(x * y, x + y)  # additional check
  fs <- list(cbind2, rbind2, g)
  for (f in fs) {
    for (i in 2:10) {
      test_AD_with(
        f, input = list(x = randn(i, i), y = randn(i, i)),
        ctrl = list(display = F, err_fun = rel_err, epsilon = 1e-7)
      )
    }
  }
})
