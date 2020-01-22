testthat::test_that("Fast kronecker product", {
  for (i in 1:10) {
    A <- randn(30, 40)
    B <- randn(60, 70)
    Z <- randn(2800, 10)
    testthat::expect_equal((B %x% A) %*% Z, BxAZ(B, A, Z))

    A <- randn(30, 40)
    B <- randn(60, 70)
    X <- randn(50, 1800)
    testthat::expect_equal(X %*% (B %x% A), XBxA(X, B, A))

    A <- randn(50, 1800)
    B <- randn(60, 70)
    I <- diag(30)
    C <- randn(60, 70)
    D <- randn(2100, 50)
    testthat::expect_equal((B %x% I) %*% D, BxID(B, D))
    testthat::expect_equal((I %x% C) %*% D, IxCD(C, D))
    testthat::expect_equal(A %*% (B %x% I), ABxI(A, B))
    testthat::expect_equal(A %*% (I %x% C), AIxC(A, C))
  }
})

# # Benchmark
# s <- function() sample(10:50, 1)
# A <- randn(s(), s())
# B <- randn(s(), s())
# X <- randn(s(), nrow(A) * nrow(B))
# Z <- randn(ncol(A) * ncol(B), s())
# microbenchmark::microbenchmark((B %x% A) %*% Z, BxAZ(B, A, Z))
# microbenchmark::microbenchmark(X %*% (B %x% A), XBxA(X, B, A))
#
# B <- C <- randn(s(), s())
# I <- Diagonal0(s())
# A <- randn(s(), nrow(B) * nrow(I))
# D <- randn(ncol(B) * ncol(I), s())
# microbenchmark::microbenchmark((B %x% I) %*% D, BxID(B, D))
# microbenchmark::microbenchmark((I %x% C) %*% D, IxCD(C, D))
# microbenchmark::microbenchmark(A %*% (B %x% I), ABxI(A, B))
# microbenchmark::microbenchmark(A %*% (I %x% C), AIxC(A, C))
