#' @include class_dual_def.R
NULL

#' Inverse of 'dual'-class objects
#' @param a A "dual" object.
setMethod("solve",
  signature(a = "dual"),
  function(a) {
    inv_a <- solve(parent_of(a))
    new("dual", x = inv_a, dx = d_solve(a, inv_a), param = param_of(a))
  }
)

#' Transpose of 'dual'-class objects
#' @param x A "dual" object.
setMethod("t",
  signature(x = "dual"),
  function(x) {
    new("dual", x = t(parent_of(x)), dx = d_transpose(x), param = param_of(x))
  }
)

#' Crossproduct of 'dual'-class objects
#' @param x A "dual" object.
#' @param y A "dual" object.
setMethod("tcrossprod",
  signature(x = "dual", y = "missing"),
  function(x, y) {
    new("dual",
        x = tcrossprod(parent_of(x)),
        dx = d_XXT(x),
        param = param_of(x))
  }
)

#' Crossproduct of 'dual'-class objects
#' @param x A "dual" object.
#' @param y A "dual" object.
setMethod("tcrossprod",
  signature(x = "dual", y = "dual"),
  function(x, y) { x %*% t(y) }
)


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
setMethod("crossprod",
  signature(x = "dual", y = "dual"),
  function(x, y) { t(x) %*% y }
)


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
    L <- chol0(parent_of(x))
    dL <- d_chol(L, x)
    x@x <- L
    x@dx <- dL
    x
  }
)


#' Determinant of a matrix
#' This is a wrapper of `base::det`.
#' @param x A numeric matrix.
#' @param ... Other parameters to be passed to `base::det`.
#' @export
det0 <- function(x, ...) { det(x, ...) }

#' Determinant of a 'dual'-class object
#' @param x A "dual" object.
setMethod("det0",
  signature(x = "dual"),
  function(x) {
    px <- parent_of(x)
    det_x <- det(px)
    x@x <- det_x
    x@dx <- det_x * t(as.numeric(t(solve(px))))
    x
  }
)
