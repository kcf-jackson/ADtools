#' @include class_dual_def.R
NULL

#' Length of an Object
#' @param x A "dual" object.
setMethod("length",
          signature(x = "dual"),
          function(x) { length(parent_of(x)) }
)

#' Coerce the first component of the dual object into a vector.
#' @param x A "dual" object.
setMethod("as.vector",
          signature(x = "dual"),
          function(x) {
            x@x <- as.vector(parent_of(x))
            x
          }
)


#'  Coerce the first component of the dual object into a matrix.
#' @method as.matrix dual
#' @param x A "dual" object.
#' @param ... Use `nrow`, `ncol` to specify the number of rows and columns.
#' @export
as.matrix.dual <- function(x, ...) {
  x@x <- matrix(parent_of(x), ...)
  x
}

#'  Coerce the first component of the dual object into a matrix.
#' @param x A "dual" object.
#' @param ... Use `nrow`, `ncol` to specify the number of rows and columns.
setMethod("as.matrix",
          signature(x = "dual"),
          as.matrix.dual
)
