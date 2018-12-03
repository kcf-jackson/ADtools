testthat::context("Test class 'dual'")

testthat::test_that("Test Constructor and Extractor", {
  a <- dual(x = randn(2,2), list(x = 3, y = 4), 2)
  testthat::expect_s4_class(a, "dual")
  testthat::expect_true(is.numeric(a@x))
  testthat::expect_equal(dim(a@dx), c(4, 7))

  b <- dual(x = randn(2,2), c(x = 1, y = 4), 2)
  testthat::expect_s4_class(b, "dual")
  testthat::expect_true(is.numeric(b@x))
  testthat::expect_equal(dim(b@dx), c(4, 5))

  testthat::expect_error(new("dual", x = randn(2,2), dx = "Error!"))
})
