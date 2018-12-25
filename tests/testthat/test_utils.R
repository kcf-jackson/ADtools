testthat::test_that("Single-layer destructuring", {
  c(x, y) %<-% list(1, 2)
  testthat::expect_equal(x, 1)
  testthat::expect_equal(y, 2)

  c(x, y) %<-% list(4, 5, 6)
  testthat::expect_equal(x, 4)
  testthat::expect_equal(y, 5)

  testthat::expect_error(c(x, y) %<-% list(3))
  testthat::expect_error(list(x, y) %<-% list(3, 3))
  testthat::expect_error(c(x, y) %<-% c(3, 6))
})