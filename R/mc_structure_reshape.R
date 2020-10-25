# #' Diagonal matrix
# #' @inheritParams base::diag
# #' @param x A matrix, a vector or a dual number. In the last case, all other parameters
# #' are ignored as \code{diag.dual} only takes the first argument.
# #' @export
# diag <- function(x, ...) {
#   UseMethod("diag", x, ...)
# }

# #' Diagonal matrix
# #' @inheritParams base::diag
# #' @export
# diag <- function(x, nrow, ncol, names = TRUE) {
#   if (class(x) == "dual") {
#     return(diag.dual(x))
#   }
#   base::diag(x, nrow, ncol, names)
# }

#' Diagonal matrix
#' @param x A "dual" object.
#' @rdname diag_dual
#' @note A single column matrix is treated as a vector. If one wants to call
#' \code{diag} on a single column matrix \code{x}, one can call
#' \code{x[1, 1]} instead.
#' @export
diag.dual <- function(x) {
  x@dx <- d_diagonal(x)
  m0 <- x@x
  if (is_matrix(m0)) {
    if (ncol(m0) == 1) {
      m0 <- as.vector(m0)  # then proceeds as the vector case
    } else {
      x@x <- diag(m0)
      return(x)
    }
  }
  if (is_vector(m0)) {
    x@x <- Diagonal0(x = m0)
    return(x)
  }
  if (is_scalar(m0)) {
    x@x <- Diagonal0(n = m0)
    return(x)
  }
  stop("The input must be a (dual) matrix / vector / scalar.")
}

d_diagonal <- function(x) {
  m0 <- x@x
  dx <- x@dx
  if (is_matrix(m0) && (ncol(m0) > 1)) {
    # extract from dx the corresponding diagonal entries
    diag_ind <- seq(1, length(m0), nrow(m0) + 1)
    return(dx[diag_ind, , drop = FALSE])
  }

  m0 <- as.vector(m0)
  if (is_scalar(m0)) {
    return(zero_matrix0(m0^2, ncol(dx)))
  }

  # handle both a single-column matrix and a vector
  new_dx <- zero_matrix0(length(m0)^2, ncol(dx))
  diag_ind <- seq(1, length(m0)^2, length(m0) + 1)
  new_dx[diag_ind, ] <- dx
  return(new_dx)
}

#' @rdname diag_dual
setMethod("diag", signature(x = "dual"), diag.dual)




#' Vectorisation
#' 
#' @param x A matrix.
#' 
#' @examples 
#' A <- randn(3, 3)
#' vec(A)
#' 
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
#' 
#' @param x A matrix.
#' 
#' @examples 
#' A <- randn(3, 3)
#' vech(A)
#' 
#' @export
vech <- function(x) {
  if (nrow(x) != ncol(x)) {
    stop("Input should be a square matrix")
  }
  as.matrix(x[lower.tri(x, diag = TRUE)])
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
    ind <- which(lower.tri(x@x, diag = TRUE))
    x <- vec(x)
    x[ind, , drop = FALSE]
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
#' @rdname as.matrix.dual
setMethod("as.matrix", signature(x = "dual"), as.matrix.dual)




#' Coerce the first component of the dual object into a matrix.
#' @method matrix dual
#' @inheritParams matrix
#' @rdname matrix.dual
#' @export
matrix.dual <- function(data, nrow, ncol = 1, byrow = FALSE, dimnames = NULL) {
  if (missing(nrow)) {
    nrow <- round(length(data@x) / ncol)
  }
  ncol <- round(length(data@x) / nrow)
  data@x <- matrix(data@x,
    nrow = nrow, ncol = ncol,
    byrow = byrow, dimnames = dimnames
  )
  if (byrow) {
    ind <- as.numeric(matrix(seq_along(data@x), nrow = nrow, ncol = ncol, byrow = TRUE))
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
#' @param byrow TRUE or FALSE; whether to fill the matrix by rows.
#' @param dimnames A dimnames attribute for the matrix: NULL or
#' a list of length 2 giving the row and column names respectively.
#' An empty list is treated as NULL, and a list of length one as
#' row names. The list can be named, and the list names will be
#' used as names for the dimensions.
#' @export
matrix <- function(data, ...) {
  UseMethod("matrix", data, ...)
}




#' Construct a lower triangular matrix from a vector
#'
#' @rdname lower-triangular
#'
#' @param data A numeric vector.
#' @param nrow A positive integer; the desired number of rows.
#' @param ncol A positive integer; the desired number of rows.
#' @param diag A logical; should the diagonal be included ?
#'
#' @examples
#' lower_tri_matrix(1:3, 3, 3)
#' lower_tri_matrix(1:6, 3, 3, diag = TRUE)
#' 
#' @export
lower_tri_matrix <- function(data, nrow = 1, ncol = 1, diag = FALSE) {
  y <- matrix(0, nrow = nrow, ncol = ncol)
  y[lower.tri(y, diag = diag)] <- data
  y
}

#' @rdname lower-triangular
setMethod(
  "lower_tri_matrix",
  signature(data = "dual"),
  function(data, nrow = 1, ncol = 1, diag = FALSE) {
    x <- lower_tri_matrix(data@x, nrow = nrow, ncol = ncol, diag)
    dx <- matrix(0, nrow = nrow * ncol, ncol = ncol(data@dx))

    ind <- which(lower.tri(x, diag = diag))
    dx[ind, ] <- as.matrix(data@dx)

    new("dual", x = x, dx = dx)
  }
)




#' Construct an upper triangular matrix from a vector
#' 
#' @rdname upper-triangular
#' 
#' @param data A numeric vector.
#' @param nrow A positive integer; the desired number of rows.
#' @param ncol A positive integer; the desired number of rows.
#' @param diag A logical; should the diagonal be included ?
#' 
#' @examples
#' upper_tri_matrix(1:3, 3, 3)
#' upper_tri_matrix(1:6, 3, 3, diag = TRUE)
#' 
#' @export
upper_tri_matrix <- function(data, nrow = 1, ncol = 1, diag = FALSE) {
  y <- matrix(0, nrow = nrow, ncol = ncol)
  y[upper.tri(y, diag = diag)] <- data
  y
}

#' @rdname upper-triangular
setMethod(
  "upper_tri_matrix",
  signature(data = "dual"),
  function(data, nrow = 1, ncol = 1, diag = FALSE) {
    x <- upper_tri_matrix(data@x, nrow = nrow, ncol = ncol, diag)
    dx <- matrix(0, nrow = nrow * ncol, ncol = ncol(data@dx))

    ind <- which(upper.tri(x, diag = diag))
    dx[ind, ] <- as.matrix(data@dx)

    new("dual", x = x, dx = dx)
  }
)
