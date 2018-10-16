testthat::context("Test class 'dual'")

testthat::test_that("Test Constructor and Extractor", {
  a <- new(
    "dual", x = randn(2,2),
    dx = init_dx(list(c(4, 5), c(4, 2)), 1)
  )
  testthat::expect_s4_class(a, "dual")
  testthat::expect_true(is.numeric(a@x))
  testthat::expect_true((class(a@dx) == "matrix_list") || (a@dx == 0))

  b <- new("dual", x = randn(2,2), dx = 0)
  testthat::expect_true(is.numeric(b@x))
  testthat::expect_true((class(b@dx) == "matrix_list") || (b@dx == 0))

  testthat::expect_error(new("dual", x = randn(2,2), dx = "Error!"))
})

# Test Arithmetic -> test_matrix_calculus.R
