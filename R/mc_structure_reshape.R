#' Diagonal matrix
#' @inheritParams base::diag
#' @param x A matrix, a vector or a dual number. In the last case, all other parameters
#' are ignored as \code{diag.dual} only takes the first argument.
#' @export
diag <- function(x, ...) {
  UseMethod("diag", x, ...)
}

# diag <- function(x, nrow, ncol, names = TRUE) {
#   if (class(x) == "dual") {
#     return(diag.dual(x))
#   }
#   base::diag(x, nrow, ncol, names)
# }

#' Diagonal matrix
#' @param x A "dual" object.
#' @rdname diag_dual
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

#' @rdname diag_dual
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
  if (nrow(x) != ncol(x)) {
    stop("Input should be a square matrix")
  }
  as.matrix(x[lower.tri(x, diag = T)])
}
# References
# 1. matrixcalc::vech
# 2. https://en.wikipedia.org/wiki/Vectorization_(mathematics)#Half-vectorization


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


#' Coerce the first component of the dual object into a matrix.
#' @name as.matrix.dual
#' @method as.matrix dual
#' @param x A "dual" object.
#' @param ... Unused.
#' @export
as.matrix.dual <- function(x, ...) {
  x@x <- as.matrix(x@x)
  x
}
setMethod("as.matrix", signature(x = "dual"), as.matrix.dual)


#' Coerce the first component of the dual object into a matrix.
#' @method matrix dual
#' @inheritParams matrix
#' @rdname matrix.dual
#' @export
matrix.dual <- function(data, nrow, ncol = 1, byrow = F, dimnames = NULL) {
  if (missing(nrow)) {
    nrow <- round(length(data@x) / ncol)
  }
  ncol <- round(length(data@x) / nrow)
  data@x <- matrix(data@x,
    nrow = nrow, ncol = ncol,
    byrow = byrow, dimnames = dimnames
  )
  if (byrow) {
    ind <- as.numeric(matrix(seq_along(data@x), nrow = nrow, ncol = ncol, byrow = T))
    data@dx <- data@dx[ind, ]
  }
  data
}

#' @rdname matrix.dual
setMethod("matrix", signature(data = "dual"), matrix.dual)


#' Matrices
#' @param data A "dual" object.
#' @param nrow a positive integer; the desired number of rows.
#' @param ncol a positive integer; the desired number of rows.
#' @param byrow T or F; whether to fill the matrix by rows.
#' @param dimnames A dimnames attribute for the matrix: NULL or
#' a list of length 2 giving the row and column names respectively.
#' An empty list is treated as NULL, and a list of length one as
#' row names. The list can be named, and the list names will be
#' used as names for the dimensions.
#' @export
matrix <- function(data, ...) {
  UseMethod("matrix", data, ...)
}
