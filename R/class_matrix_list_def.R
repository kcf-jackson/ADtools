#=========================================================================
# Class "matrix_list"
#-------------------------------------------------------------------------

is_matrix <- function(x) {
  is.matrix(x) || (attr(class(x), "package") == "Matrix")
}

all_matrices <- function(object) {
  all(purrr::map_lgl(object@matrices, is_matrix))
}

#' S4 class "matrix_list"
#'
#' @name matrix-list-class
#' @slot matrices A list of matrices. Each matrix must be of class 'matrix' or
#' an object from the package "Matrix.
#' @import methods
#' @examples
#' library(ADtools)
#' a <- new("matrix_list", matrices = purrr::map(1:4, ~randn(2,2)))
setClass(
  "matrix_list",
  representation(matrices = "list"),
  validity = all_matrices
)

element_of <- function(x) { x@matrices }
