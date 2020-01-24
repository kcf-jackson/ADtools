testthat::test_that("Test that the correct minimum cost for MCM is found", {
  expect_equal(mcm_optimal_order(c(2, 1, 2, 5))$cost          , 20)
  expect_equal(mcm_optimal_order(c(1, 10, 2, 10))$cost        , 40)
  expect_equal(mcm_optimal_order(c(30, 1, 40, 10, 25))$cost   , 1400)
  expect_equal(mcm_optimal_order(c(10, 20, 50, 1, 100))$cost  , 2200)
  expect_equal(mcm_optimal_order(c(2, 6, 1, 1))$cost          , 14)
  expect_equal(mcm_optimal_order(c(5, 2, 1, 1, 5))$cost       , 37)
  expect_equal(mcm_optimal_order(c(10, 100, 5, 50, 1))$cost   , 1750)
})
