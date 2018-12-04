testthat::context("Test special matrix construction")

testthat::test_that("Memoisation", {
  mdiag <- ADtools:::memoize(diag)
  first_runtime <- system.time(mdiag(1000))[3]
  second_runtime <-  system.time(mdiag(1000))[3]
  testthat::expect_true(first_runtime > second_runtime)
})

testthat::test_that("Commutation matrix", {
  testthat::expect_equal(commutation_matrix0(3), commutation_matrix0(3, 3))
})
