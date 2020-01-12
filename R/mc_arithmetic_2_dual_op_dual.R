#' @include class_dual_def.R
NULL

#' Addition of 'dual'-class objects
#' @param e1 A "dual" object.
#' @param e2 A "dual" object.
setMethod(
  "+",
  signature(e1 = "dual", e2 = "dual"),
  function(e1, e2) {
    x <- e1@x + e2@x
    dx <- d_sum(e1, e2)
    e1@x <- x
    e1@dx <- dx
    e1
  }
)

#' Subtraction of 'dual'-class objects
#' @param e1 A "dual" object.
#' @param e2 A "dual" object.
setMethod(
  "-",
  signature(e1 = "dual", e2 = "dual"),
  function(e1, e2) {
    x <- e1@x - e2@x
    new("dual", x = x, dx = d_minus(e1, e2))
  }
)

#' (Element-wise) Multiplication of 'dual'-class objects
#' @param e1 A "dual" object.
#' @param e2 A "dual" object.
setMethod(
  "*",
  signature(e1 = "dual", e2 = "dual"),
  function(e1, e2) {
    x <- e1@x * e2@x
    dx <- d_scalar_prod(e1, e2)
    e1@x <- x
    e1@dx <- dx
    e1
  }
)

#' (Element-wise) Division of 'dual'-class objects
#' @param e1 A "dual" object.
#' @param e2 A "dual" object.
setMethod(
  "/",
  signature(e1 = "dual", e2 = "dual"),
  function(e1, e2) {
    x <- e1@x / e2@x
    new("dual", x = x, dx = d_divide(e1, e2))
  }
)

#' Matrix multiplication of 'dual'-class objects
#' @param x A "dual" object.
#' @param y A "dual" object.
setMethod(
  "%*%",
  signature(x = "dual", y = "dual"),
  function(x, y) {
    res <- x@x %*% y@x
    new("dual", x = res, dx = d_matrix_prod(x, y))
  }
)

#' Kronecker multiplication of 'dual'-class objects
#' @param X A "dual" object.
#' @param Y A "dual" object.
setMethod(
  "kronecker",
  signature(X = "dual", Y = "dual"),
  function(X, Y) {
    res <- X@x %x% Y@x
    new("dual", x = res, dx = d_kronecker(X, Y))
  }
)

#' Powers of 'dual'-class objects
#' @param e1 A "dual" object.
#' @param e2 A "dual" object.
setMethod(
  "^",
  signature(e1 = "dual", e2 = "dual"),
  function(e1, e2) {
    if (any(e1@x <= 0)) {
      stop("When the exponent is a dual number, the base must be a matrix with positive entries to use e1 ^ e2 = exp(e2 * log(e1))")
    }
    if (!is_scalar(e2@x)) {
      stop("The exponent must be a scalar.")
    }
    px <- e1@x ^ e2@x
    dx <- as.numeric(px) * (e2 * log(e1))@dx
    e1@x <- px
    e1@dx <- dx
    e1
  }
)
