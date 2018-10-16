#' Initialise the dual component, i.e. the storage matrices for derivatives
#'
#' @param dx_dims list of pairs; the list of matrix dimensions.
#' @param dx_ind An integer that indicates the location of the derivative
#' corresponding to x; input -1 if none.
#'
#' @examples
#' init_dx(list(c(4, 2), c(4, 5), c(4, 1)), -1)
#' init_dx(list(c(4, 2), c(4, 5), c(4, 1)), 1)
#' init_dx(list(c(4, 2), c(4, 5), c(4, 1)), 3)
#'
#' @export
init_dx <- function(dx_dims, dx_ind) {
  deriv <- purrr::map(dx_dims, zero_matrix)
  if (dx_ind != -1) {
    deriv[[dx_ind]] <- diag_matrix(dx_dims[[dx_ind]])
  }
  deriv <- new("matrix_list", matrices = deriv)
}

zero_matrix <- function(dim0) {
  Matrix::Matrix(0, nrow = dim0[1], ncol = dim0[2])
}

diag_matrix <- function(dim0, x = 1) {
  seq_n <- seq(min(dim0))
  Matrix::sparseMatrix(seq_n, seq_n, x = x, dims = dim0)
}

#' Generate a matrix randomly
#' @param ... Any number of dimensions.
#' @export
randn <- function(...) {
  structure(rnorm(prod(...)), dim = c(...))
}
