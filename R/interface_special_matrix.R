#' @include optim_special_matrix.R
NULL

#' Matrix of zeroes (Memoised)
#' @param ...
#' `nr` integer; number of rows.
#'
#' `nc` integer; number of columns.
#' @export
zero_matrix0 <- memo_zero_matrix


#' Matrix of ones (Memoised)
#' @param ...
#' `nr` integer; number of rows.
#'
#' `nc` integer; number of columns.
#' @export
one_matrix0 <- memo_one_matrix


#' Band matrix (Memoised)
#' @param ...
#' `nr` integer; number of rows.
#'
#' `nc` integer; number of columns.
#'
#' `x` A scalar or a vector to be placed on the diagonal of the matrix.
#' @export
band_matrix0 <- memo_band_matrix


#' Diagonal matrix (Memoised)
#' @param ...
#' `n` integer; the dimension of the square matrix.
#'
#' `x` A scalar or a vector to be placed on the diagonal of the matrix.
#' @export
Diagonal0 <- memo_Diagonal


#' Commutation matrix (Memoised)
#' @param ...
#' `r` integer; row dimension.
#'
#' `c` integer; column dimension.
#' @export
commutation_matrix0 <- memo_commutation_matrix


#' Elimination matrix (Memoised)
#' @param ...
#' `n` integer, row or column dimension.
#' @export
elimination_matrix0 <- memo_elimination_matrix
