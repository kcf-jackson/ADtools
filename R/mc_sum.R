#' Sum of matrix elements
#' @param x A "dual" object.
setMethod("sum",
           signature(x = "dual"),
           function(x) {
             x@dx <- t(colSums(deriv_of(x)))
             x@x <- sum(parent_of(x))
             x
           }
)

#' Row sum of a matrix.
#' @param x A "dual" object.
setMethod("rowSums",
           signature(x = "dual"),
           function(x) {
             x@dx <- d_rowSums(x)
             x@x <- rowSums(parent_of(x))
             x
           }
)

#' Column sum of a matrix.
#' @param x A "dual" object.
setMethod("colSums",
           signature(x = "dual"),
           function(x) {
             x@dx <- d_colSums(x)
             x@x <- colSums(parent_of(x))
             x
           }
)

#' Trace of a matrix
#' @param x A square matrix
#' @export
tr <- function(x) {
  if (nrow(x) != ncol(x)) stop("Input must be a square matrix")
  sum(diag(x))
}

#' Trace of a matrix
#' @param x A "dual" object.
setMethod("tr",
           signature(x = "dual"),
           function(x) { sum(diag(x)) }
)
