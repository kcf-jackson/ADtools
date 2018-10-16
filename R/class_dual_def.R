#' @include class_matrix_list_def.R
NULL

#=========================================================================
# Class "dual"
#-------------------------------------------------------------------------

setClassUnion("list_matrices_or_zero", c("matrix_list", "numeric"))
setClassUnion("array_or_numeric", c("array", "numeric"))

check_dual <- function(object) {
  x <- object@x
  dx <- object@dx
  x_check <- is.array(x) || is.numeric(x)
  dx_check <- ifelse(is_zero(dx), T, match_dim(x, dx))
  x_check && dx_check
}

match_dim <- function(x, dx) {
  # The row number of each element of dx must match the length of x
  all(Map(nrow, dx@matrices) == length(x))
}

#' S4 class "dual"
#'
#' @description This class attaches a dual component to a number / an array.
#' @slot x Any "numeric", i.e. vector, matrix or array.
#' @slot dx 0 or class "matrix_list", see \link{matrix-list-class} for more detail.
#' @import methods
#' @examples
#' a <- new("dual", x = randn(2,2), dx = 0)
#' b <- new("dual", x = randn(2,2), dx = init_dx(list(c(4, 1), c(4, 2)), 2))
setClass(
  "dual",
  representation(
    x = "array_or_numeric",
    dx = "list_matrices_or_zero"
  ),
  validity = check_dual
)

parent_of <- function(x) { x@x }
deriv_of <- function(x) { x@dx }
