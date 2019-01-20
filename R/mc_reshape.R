#' Diagonal matrix
#' @param x A matrix, or a vector.
#' @export
diag.dual <- function(x) {
  x@dx = d_diagonal(x)
  x@x = diag(parent_of(x))
  x
}

#' Diagonal matrix
#' @param x A "dual" object.
setMethod("diag", signature(x = "dual"), diag.dual)


#' Vectorisation
#' @param x A matrix.
#' @export
vec <- function(x) {
  as.matrix(as.numeric(x))
}


#' Vectorisation
#' @param x A "dual" object.
setMethod("vec",
          signature(x = "dual"),
          function(x) {
            # vec doesn't change the derivative
            x@x <- as.matrix(as.numeric(parent_of(x)))
            x
          }
)


#' Half-vectorisation
#' @param x A matrix.
#' @export
vech <- function(x) {
  if (nrow(x) != ncol(x)) stop("Input should be a square matrix")
  as.matrix(x[lower.tri(x, diag = T)])
}


#' Half-vectorisation
#' @param x A "dual" object.
setMethod("vech",
          signature(x = "dual"),
          function(x) {
            ind <- which(lower.tri(parent_of(x), diag = T))
            x <- vec(x)
            x[ind, , drop = F]
          }
)
