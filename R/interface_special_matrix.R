#' @include optim_special_matrix.R
NULL

#' Matrix of zeroes (memoised)
#'
#' @param ...
#' `nr` integer; number of rows.
#'
#' `nc` integer; number of columns.
#'
#' Use empty argument to clear the cache.
#'
#' @export
zero_matrix0 <- memo_zero_matrix


#' Matrix of ones (memoised)
#'
#' @param ...
#' `nr` integer; number of rows.
#'
#' `nc` integer; number of columns.
#'
#' Use empty argument to clear the cache.
#'
#' @export
one_matrix0 <- memo_one_matrix


#' Band matrix (memoised)
#'
#' @param ...
#' `nr` integer; number of rows.
#'
#' `nc` integer; number of columns.
#'
#' `x` A scalar or a vector to be placed on the diagonal of the matrix.
#'
#' Use empty argument to clear the cache.
#'
#' @export
band_matrix0 <- memo_band_matrix


#' Diagonal matrix (memoised)
#'
#' @param ...
#' `n` integer; the dimension of the square matrix.
#'
#' `x` A scalar or a vector to be placed on the diagonal of the matrix.
#'
#' Use empty argument to clear the cache.
#'
#' @export
Diagonal0 <- memo_Diagonal


#' Commutation matrix (memoised)
#'
#' @param ...
#' `r` integer; row dimension.
#'
#' `c` integer; column dimension.
#'
#' Use empty argument to clear the cache.
#'
#' @export
commutation_matrix0 <- memo_commutation_matrix


#' Elimination matrix (memoised)
#'
#' @param ...
#' `n` integer, row or column dimension.
#'
#' Use empty argument to clear the cache.
#'
#' @export
elimination_matrix0 <- memo_elimination_matrix
