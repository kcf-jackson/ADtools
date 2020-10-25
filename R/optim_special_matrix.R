#' Commutation matrix
#'
#' @description Suppose A is a (r x c) matrix, then a (rc x rc) matrix K
#' is a commutation matrix if K %*% vec(A) = vec(t(A)).
#'
#' @param r integer; row dimension.
#' @param c integer; column dimension.
#'
#' @examples
#' A <- randn(3, 4)
#' K <- commutation_matrix(3, 4)
#' compare <- function(x, y) {
#'     cbind(x, y, err = abs(x - y))
#' }
#' compare(K %*% vec(A), vec(t(A)))
#'
#' @export
commutation_matrix <- function(r, c) {
  if (missing(c)) c <- r
  entries <- expand.grid(1:r, 1:c)
  src <- r * (entries[,2] - 1) + entries[,1]
  tgt <- c * (entries[,1] - 1) + entries[,2]
  Matrix::sparseMatrix(tgt, src, x = 1)
}


#' Elimination matrix
#'
#' @description Suppose A is a (n x n) symmetric matrix, then a (n(n-1)/2 x n^2) matrix E
#' is an elimination matrix if E %*% vec(A) = vech(A).
#'
#' @param n integer, row or column dimension.
#'
#' @examples
#' A <- crossprod(randn(3, 3))
#' E <- elimination_matrix(3)
#' compare <- function(x, y) {
#'     cbind(x, y, err = abs(x - y))
#' }
#' compare(E %*% vec(A), vech(A))
#'
#' @export
elimination_matrix <- function(n) {
  k <- 1:(n*(n+1)*0.5)
  A <- n + 0.5
  b <- ceiling(A - sqrt(A^2 - 2*k))
  Matrix::sparseMatrix(i = k, j = k + 0.5 * b * (b - 1), x = 1)
}


#' Matrix of zeroes
#'
#' @param nr integer; row dimension.
#' @param nc integer; column dimension.
#'
#' @examples
#' zero_matrix(1, 3)
#' zero_matrix(5, 1)
#' zero_matrix(2, 2)
#'
#' @export
zero_matrix <- function(nr, nc) {
  Matrix::Matrix(data = 0, nrow = nr, ncol = nc)
}


#' Matrix of ones
#'
#' @param nr integer; row dimension.
#' @param nc integer; column dimension.
#'
#' @examples
#' one_matrix(1, 3)
#' one_matrix(5, 1)
#' one_matrix(2, 2)
#'
#' @export
one_matrix <- function(nr, nc) {
  Matrix::Matrix(data = 1, nrow = nr, ncol = nc)
}


#' Band matrix
#'
#' @param nr integer; row dimension.
#' @param nc integer; column dimension.
#' @param x A scalar or a vector to be placed on the diagonal of the matrix.
#'
#' @examples
#' band_matrix(2, 3)
#' band_matrix(2, 3, c(999))
#' band_matrix(2, 3, c(999, 111))
#'
#' @export
band_matrix <- function(nr, nc, x = 1) {
  dim0 <- c(nr, nc)
  seq_n <- seq(min(dim0))
  Matrix::sparseMatrix(seq_n, seq_n, x = x, dims = dim0)
}


#' Diagonal matrix
#'
#' @param n integer; the dimension of the square matrix.
#' @param x A scalar or a vector to be placed on the diagonal of the matrix.
#'
#' @examples
#' Diagonal(3)
#' Diagonal(3, 999)
#' Diagonal(3, c(11, 22, 33))
#'
#' @export
Diagonal <- Matrix::Diagonal


memoize <- function(f) {
  table0 <- new.env()  # local cache
  char_hash <- function(...) {
    paste0("hash_", paste0(as.character(list(...)), collapse = ","))
  }
  function(...) {
    # Clear local cache if no argument is provided
    if (length(list(...)) == 0) {
      rm(list = names(table0), envir = table0)
      return(TRUE)
    }
    # Retrieve from cache or compute on the fly
    char_x <- char_hash(...)
    res <- table0[[char_x]]
    if (is.null(res)) {
      res <- f(...)
      table0[[char_x]] <- res
    }
    return(res)
  }
}

memo_zero_matrix <- memoize(zero_matrix)
memo_one_matrix <- memoize(one_matrix)
memo_band_matrix <- memoize(band_matrix)
memo_Diagonal <- memoize(Diagonal)
memo_commutation_matrix <- memoize(commutation_matrix)
memo_elimination_matrix <- memoize(elimination_matrix)
