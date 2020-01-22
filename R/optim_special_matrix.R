commutation_matrix <- function(r, c) {
  if (missing(c)) c <- r
  entries <- expand.grid(1:r, 1:c)
  src <- r * (entries[,2] - 1) + entries[,1]
  tgt <- c * (entries[,1] - 1) + entries[,2]
  Matrix::sparseMatrix(tgt, src, x = 1)
}

elimination_matrix <- function(n) {
  k <- 1:(n*(n+1)*0.5)
  A <- n + 0.5
  b <- ceiling(A - sqrt(A^2 - 2*k))
  Matrix::sparseMatrix(i = k, j = k + 0.5 * b * (b - 1), x = 1)
}

zero_matrix <- function(nr, nc) {
  Matrix::Matrix(data = 0, nrow = nr, ncol = nc)
}

one_matrix <- function(nr, nc) {
  Matrix::Matrix(data = 1, nrow = nr, ncol = nc)
}

band_matrix <- function(nr, nc, x = 1) {
  dim0 <- c(nr, nc)
  seq_n <- seq(min(dim0))
  Matrix::sparseMatrix(seq_n, seq_n, x = x, dims = dim0)
}

memoize <- function(f) {
  table0 <- list()
  char_hash <- function(...) {
    paste0(as.character(list(...)), collapse = ",")
  }
  function(...) {
    char_x <- char_hash(...)
    res <- table0[[char_x]]
    if (is.null(res)) {
      res <- f(...)
      table0[[char_x]] <<- res
    }
    return(res)
  }
}

memo_zero_matrix <- memoize(zero_matrix)
memo_one_matrix <- memoize(one_matrix)
memo_band_matrix <- memoize(band_matrix)
memo_Diagonal <- memoize(Matrix::Diagonal)
memo_commutation_matrix <- memoize(commutation_matrix)
memo_elimination_matrix <- memoize(elimination_matrix)
