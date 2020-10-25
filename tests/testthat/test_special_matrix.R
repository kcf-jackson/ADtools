testthat::context("Test special matrix construction")

testthat::test_that("Memoisation", {
  mdiag <- ADtools:::memoize(diag)
  first_runtime <- system.time(m1 <- mdiag(4000))[3]
  second_runtime <-  system.time(m2 <- mdiag(4000))[3]
  testthat::expect_identical(m1, m2)
  testthat::expect_true(first_runtime > second_runtime)
})

expect_equal_matrix <- function(x, y) {
  # Convert a sparse matrix into a dense matrix and
  # remove the hidden list of dimnames.
  as_matrix <- function(x) {
    x <- as.matrix(x)
    dimnames(x) <- NULL
    x
  }
  testthat::expect_equal(as_matrix(x), as_matrix(y))
}

testthat::test_that("Commutation matrix", {
  A <- randn(3, 4)
  K <- commutation_matrix(3, 4)
  expect_equal_matrix(K %*% vec(A), vec(t(A)))
})

testthat::test_that("Elimination matrix", {
  A <- crossprod(randn(3, 3))
  E <- elimination_matrix(3)
  expect_equal_matrix(E %*% vec(A), vech(A))
})
