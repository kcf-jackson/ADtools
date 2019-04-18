testthat::context("Test class 'dual'")

testthat::test_that("Test Constructor and Extractor", {
  a <- dual(x = randn(2,2), list(x = 3, y = 4), 2)
  testthat::expect_s4_class(a, "dual")
  testthat::expect_true(is.numeric(parent_of(a)))
  testthat::expect_equal(dim(deriv_of(a)), c(4, 7))
  testthat::expect_true(ADtools:::check_dual(a))

  b <- dual(x = randn(2,2), c(x = 1, y = 4), 2)
  testthat::expect_s4_class(b, "dual")
  testthat::expect_true(is.numeric(parent_of(b)))
  testthat::expect_equal(dim(deriv_of(b)), c(4, 5))
  testthat::expect_true(ADtools:::check_dual(b))

  c <- dual(x = 1:10, c(x = 1, y = 4), 2)
  testthat::expect_s4_class(c, "dual")
  testthat::expect_true(is.numeric(parent_of(c)))
  testthat::expect_equal(dim(deriv_of(c)), c(10, 5))
  testthat::expect_true(ADtools:::check_dual(c))

  testthat::expect_equal(ncol(get_deriv(c)), 5)
  testthat::expect_error(new("dual", x = randn(2,2), dx = "Error!"))
})
