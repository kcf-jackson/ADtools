testthat::context("Test generic functions")

testthat::test_that("Test length", {
  A <- randn(4, 5)
  A_dual <- dual(A, param_dim = c(length(A), 20), 1)
  testthat::expect_equal(length(A_dual), length(A))
})

testthat::test_that("test as.vector, as.matrix", {
  fs <- list(as.vector, as.matrix)
  for (f in fs) {
    for (i in 2:10) {
      test_AD_with(
        f, input = list(x = randn(i, i)),
        ctrl = list(display = T, err_fun = rel_err, epsilon = 1e-7)
      )
    }
  }
})
