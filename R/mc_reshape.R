#' Diagonal matrix
setMethod("diag",
          signature(x = "dual"),
          function(x) {
            x@dx = d_diagonal(x)
            x@x = Diagonal0(parent_of(x))
            x
          }
)


#' Vectorisation
#' @param x A matrix.
#' @export
vec <- function(x) {
  as.matrix(as.numeric(x))
}


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


setMethod("vech",
          signature(x = "dual"),
          function(x) {
            ind <- which(lower.tri(parent_of(x), diag = T))
            x <- vec(x)
            x[ind, , drop = F]
          }
)
