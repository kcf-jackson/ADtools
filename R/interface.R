#' @include special_matrix.R
NULL

#' Matrix of zeroes
#' @param ...
#' `nr` integer; number of rows.
#'
#' `nc` integer; number of columns.
#' @export
zero_matrix0 <- memo_zero_matrix


#' Matrix of ones
#' @param ...
#' `nr` integer; number of rows.
#'
#' `nc` integer; number of columns.
#'
#' `x` A scalar or a vector to be placed on the diagonal of the matrix.
#' @export
one_matrix0 <- memo_one_matrix


#' Diagonal matrix
#' @param ...
#' `n` integer; the dimension of the square matrix.
#'
#' `x` A scalar or a vector to be placed on the diagonal of the matrix.
#' @export
Diagonal0 <- memo_Diagonal


#' Commutation matrix
#' @param ...
#' `r` integer; row dimension.
#'
#' `c` integer; column dimension.
#' @export
commutation_matrix0 <- memo_commutation_matrix


#' Elimination matrix
#' @param ...
#' `n` integer, row or column dimension.
#' @export
elimination_matrix0 <- memo_elimination_matrix
