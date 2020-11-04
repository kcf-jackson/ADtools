#' Generate a matrix randomly from a normal distribution
#'
#' @param ... Any number of dimensions.
#' @param mean vector of means.
#' @param sd vector of standard deviations.
#'
#' @export
randn <- function(..., mean = 0, sd = 1) {
  structure(rnorm(prod(...), mean = mean, sd = sd), dim = c(...))
}


#' Generate a matrix randomly from a uniform distribution
#'
#' @param ... Any number of dimensions.
#' @param min lower limit of the distribution. Must be finite.
#' @param max upper limit of the distribution. Must be finite.
#'
#' @export
randu <- function(..., min = 0, max = 1) {
  structure(runif(prod(...), min = min, max = max), dim = c(...))
}


mapreduce <- function(x, f, g, ...) {
  x %>% purrr::map(f, ...) %>% purrr::reduce(g)
}


map2reduce <- function(x, y, f, g, ...) {
  purrr::reduce(purrr::map2(x, y, f, ...), g)
}


if_null_then <- function(x, y) {
  if (is.null(x)) y else x
}


`%||%` <- if_null_then


scale_columns_by_vector <- function(m0, v0) {
  t(v0 * t(m0))
}


scale_rows_by_vector <- function(m0, v0) {
  v0 * m0
}


diag_v0_times_m0 <- scale_rows_by_vector


#' Add a column vector to each column of a matrix
#'
#' @param v0 Column vector
#' @param m0 Matrix
#'
# #' @examples
# #' add_vector_to_matrix_column(1:4, matrix(1:12, nrow = 4))
add_vector_to_matrix_column <- function(v0, m0) {
  n <- ncol(m0)
  m0 + v0 %*% one_matrix0(1, n)
}


#' Add a column vector to each row of a matrix
#'
#' @param v0 Column vector
#' @param m0 Matrix
#'
# #' @examples
# #' add_vector_to_matrix_row(1:4, matrix(1:12, nrow = 3))
add_vector_to_matrix_row <- function(v0, m0) {
  n <- nrow(m0)
  m0 + one_matrix0(n, 1) %*% t(v0)
}
