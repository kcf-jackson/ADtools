#' @include class_dual_def.R
NULL

#' Length of an Object
#' @param x A "dual" object.
setMethod(
  "length",
  signature(x = "dual"),
  function(x) {
    length(x@x)
  }
)


#' Coerce the first component of the dual object into a vector.
#' @param x A "dual" object.
setMethod(
  "as.vector",
  signature(x = "dual"),
  function(x) {
    x@x <- as.vector(x@x)
    x
  }
)


#'  Coerce the first component of the dual object into a matrix.
#' @method as.matrix dual
#' @param x A "dual" object.
#' @param nrow a positive integer; the desired number of rows.
#' @param ncol a positive integer; the desired number of rows.
#' @export
as.matrix.dual <- function(x, nrow, ncol) {
  x@x <- matrix(x@x, nrow = nrow, ncol = ncol)
  x
}

#'  Coerce the first component of the dual object into a matrix.
#' @param x A "dual" object.
#' @param ... Use `nrow`, `ncol` to specify the number of rows and columns.
setMethod("as.matrix", signature(x = "dual"), as.matrix.dual)
