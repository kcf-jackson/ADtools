#' Diagonal matrix
#' @param x A matrix, or a vector.
#' @export
diag.dual <- function(x) {
  x@dx <- d_diagonal(x)
  x@x <- diag(x@x)
  x
}

d_diagonal <- function(x) {
  m0 <- x@x
  dx <- x@dx
  if (is.matrix(m0)) {
    diag_ind <- seq(1, length(m0), nrow(m0) + 1)
    return(dx[diag_ind, , drop = F])
  }
  if (is.vector(m0)) {
    new_dx <- zero_matrix0(length(m0)^2, ncol(dx))
    diag_ind <- seq(1, length(m0)^2, length(m0) + 1)
    new_dx[diag_ind, ] <- dx
    return(new_dx)
  }
  stop("The input is not a matrix or a vector.")
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
setMethod(
  "vec",
  signature(x = "dual"),
  function(x) {
    # vec doesn't change the derivative
    x@x <- as.matrix(as.numeric(x@x))
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
setMethod(
  "vech",
  signature(x = "dual"),
  function(x) {
    ind <- which(lower.tri(x@x, diag = T))
    x <- vec(x)
    x[ind, , drop = F]
  }
)
