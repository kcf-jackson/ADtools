#' @include class_dual_def.R
NULL

#' Inverse of 'dual'-class objects
#' @method solve dual
#' @name solve.dual
#' @param a A "dual" object.
#' @param b A "dual" object.
#' @param ... Other arguments passed to 'base::solve'. See '?solve' for detail.
#' @export
solve.dual <- function(a, b, ...) {
  if (missing(b)) {
    fun <- call_S4("solve", a = "dual", b = "missing")
  } else {
    fun <- solve_dual
  }
  fun(a, b, ...)
}

#' @rdname solve.dual
setMethod("solve",
  signature(a = "dual", b = "missing"),
  function(a, b, ...) {
    inv_a <- solve(a@x, ...)
    a@dx <- d_solve(a, inv_a)
    a@x <- inv_a
    a
  }
)

solve_dual <- function(a, b, ...) {
  solve(a, ...) %*% b
}

#' @rdname solve.dual
setMethod("solve", signature(a = "dual", b = "dual"), solve_dual)

#' @rdname solve.dual
setMethod("solve", signature(a = "ANY", b = "dual"), solve_dual)

#' @rdname solve.dual
setMethod("solve", signature(a = "dual", b = "ANY"), solve_dual)


#' Transpose of 'dual'-class objects
#' @param x A "dual" object.
setMethod("t",
  signature(x = "dual"),
  function(x) {
    x@dx <- d_transpose(x)
    x@x <- t(x@x)
    x
  }
)

#' Crossproduct of 'dual'-class objects
#' @param x A "dual" object.
#' @param y A "dual" object.
setMethod("tcrossprod",
  signature(x = "dual", y = "missing"),
  function(x, y) {
    x@dx <- d_XXT(x)
    x@x <- tcrossprod(x@x)
    x
  }
)


tcrossprod_dual <- function(x, y) { x %*% t(y) }

#' Crossproduct of 'dual'-class objects
#' @param x A "dual" object.
#' @param y Numeric matrix.
setMethod("tcrossprod", signature(x = "dual", y = "ANY"), tcrossprod_dual)

#' Crossproduct of 'dual'-class objects
#' @param x A "dual" object.
#' @param y A "dual" object.
setMethod("tcrossprod", signature(x = "dual", y = "dual"), tcrossprod_dual)


crossprod_dual <- function(x, y) { t(x) %*% y }

#' Crossproduct of 'dual'-class objects
#' @param x A "dual" object.
#' @param y A "dual" object.
setMethod("crossprod",
  signature(x = "dual", y = "missing"),
  function(x, y) { t(x) %*% x }
)

#' Crossproduct of 'dual'-class objects
#' @param x A "dual" object.
#' @param y A "dual" object.
setMethod("crossprod", signature(x = "dual", y = "ANY"), crossprod_dual)

#' Crossproduct of 'dual'-class objects
#' @param x A "dual" object.
#' @param y A "dual" object.
setMethod("crossprod", signature(x = "dual", y = "dual"), crossprod_dual)


#' Cholesky decomposition
#' @param x A numeric matrix.
#' @note This function returns a lower-triangular matrix.
#' @export
chol0 <- function(x) { t(chol(x)) }

#' Cholesky decomposition of 'dual'-class objects
#' @param x A "dual" object.
#' @note The Cholesky decomposition used in this function returns a
#' lower triangular matrix.
setMethod("chol0",
  signature(x = "dual"),
  function(x) {
    L <- chol0(x@x)
    dL <- d_chol(L, x)
    x@x <- L
    x@dx <- dL
    x
  }
)


#' Determinant of a matrix
#' @name matrix_determinant
#' @inherit base::det
#' @export
det <- function(x, ...) {
  UseMethod("det", x)
}

#' @rdname matrix_determinant
#' @export
det.default <- base::det

#' Determinant of a 'dual'-class object
#' @param x A "dual" object.
#' @param ... Other parameters to be passed to `base::det`.
#' @export
det.dual <- function(x, ...) {
  px <- x@x
  det_x <- det(px, ...)
  x@x <- det_x
  x@dx <- det_x * t(as.numeric(t(solve(px))))
  x
}

#' Determinant of a 'dual'-class object
#' @param x A "dual" object.
#' @param ... Other parameters to be passed to `base::det`.
setMethod("det", signature(x = "dual"), det.dual)
