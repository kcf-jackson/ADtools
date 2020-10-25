#' @include class_dual_def.R
NULL

# This file contains functions to handle cross-classes arithmetic.
# Whenever an operation is performed on a dual-number and a scalar / matrix,
# it is assumed to be performed component-wise. Furthermore, all
# non-dual objects are assumed to have zero matrices as Jacobian.

#' Addition of 'dual'-class objects (dual-ANY)
#' @param e1 A "dual" object.
#' @param e2 "ANY" object.
setMethod(
  "+",
  signature(e1 = "dual", e2 = "ANY"),
  function(e1, e2) {
    fail_recycling <- (length(e1@x) %% length(e2) != 0) &&
      (length(e2) %% length(e1@x) != 0)

    if (fail_recycling) {
      stop("Dimensions of e1 and e2 do not match.")
    }

    e1@x <- e1@x + e2
    # If the size has increased after the addition, e.g. when e1
    # has length 1 and e2 has length 3, then dx should expand to
    # the appropriate length. e1@dx is duplicated because e1@x is
    # recycled.
    if (nrow(e1@dx) != length(e1@x)) {
      e1@dx <- mapreduce(numeric(length(e2) / nrow(e1@dx)), ~e1@dx, rbind)
    }
    e1
  }
)

#' Addition of 'dual'-class objects (ANY-dual)
#' @param e1 "ANY" object.
#' @param e2 A "dual" object.
setMethod(
  "+",
  signature(e1 = "ANY", e2 = "dual"),
  function(e1, e2) {
    e2 + e1
  }
)

#' Subtraction of 'dual'-class objects (dual-ANY)
#' @param e1 A "dual" object.
#' @param e2 "ANY" object.
setMethod(
  "-",
  signature(e1 = "dual", e2 = "ANY"),
  function(e1, e2) {
    e1 + (-e2)
  }
)

#' Subtraction of 'dual'-class objects (unary dual)
#' @param e1 A "dual" object.
#' @param e2 "missing" object.
setMethod(
  "-",
  signature(e1 = "dual", e2 = "missing"),
  function(e1, e2) {
    e1@x <- -e1@x
    e1@dx <- -e1@dx
    e1
  }
)

#' Subtraction of 'dual'-class objects (dual-ANY)
#' @param e1 "ANY" object.
#' @param e2 A "dual" object.
setMethod(
  "-",
  signature(e1 = "ANY", e2 = "dual"),
  function(e1, e2) {
    e1 + (-e2)
  }
)

#' (Element-wise) Multiplication of 'dual'-class objects (dual-ANY)
#' @param e1 A "dual" object.
#' @param e2 "ANY" object.
setMethod(
  "*",
  signature(e1 = "dual", e2 = "ANY"),
  function(e1, e2) {
    if (is_scalar(e2)) {
      e1@x <- e1@x * e2
      e1@dx <- e1@dx * e2
    } else {
      if (is_scalar(e1@x)) {
        # e1 is scalar-dual, e2 is matrix
        e1@x <- e1@x * e2
        e1@dx <- as.numeric(e2) %*% e1@dx
      } else {
        # e1 is matrix-dual, e2 is matrix
        e1@x <- e1@x * e2
        e1@dx <- e1@dx * as.numeric(e2)
      }
    }
    e1
  }
)

#' (Element-wise) Multiplication of 'dual'-class objects (ANY-dual)
#' @param e1 "ANY" object.
#' @param e2 A "dual" object.
setMethod(
  "*",
  signature(e1 = "ANY", e2 = "dual"),
  function(e1, e2) {
    if (is_scalar(e1)) {
      e2@x <- e1 * e2@x
      e2@dx <- e1 * e2@dx
    } else {
      if (is_scalar(e2@x)) {
        # e1 is matrix, e2 is dual-scalar
        e2@x <- e1 * e2@x
        e2@dx <- as.numeric(e1) %*% e2@dx
      } else {
        # e1 is matrix, e2 is dual-matrix
        e2@x <- e1 * e2@x
        e2@dx <- as.numeric(e1) * e2@dx
      }
    }
    e2
  }
)

#' (Element-wise) Division of 'dual'-class objects (dual-ANY)
#' @param e1 A "dual" object.
#' @param e2 "ANY" object.
setMethod(
  "/",
  signature(e1 = "dual", e2 = "ANY"),
  function(e1, e2) {
    if (is_scalar(e2)) {
      e1@x <- e1@x / e2
      e1@dx <- e1@dx / e2
    } else {
      if (is_scalar(e1@x)) {
        e1@x <- e1@x / e2
        e1@dx <- (rep(1, length(e2)) %*% e1@dx) / as.numeric(e2)
      } else {
        e1@x <- e1@x / e2
        e1@dx <- e1@dx / as.numeric(e2)
      }
    }
    e1
  }
)

#' (Element-wise) Division of 'dual'-class objects (ANY-dual)
#' @param e1 "ANY" object.
#' @param e2 A "dual" object.
setMethod(
  "/",
  signature(e1 = "ANY", e2 = "dual"),
  function(e1, e2) {
    x <- e2@x
    dx <- e2@dx
    if (is_scalar(x)) {
      inv_a <- 1 / x
      d_inv_a <- -inv_a^2 * dx
      e2@x <- e1 * inv_a
      e2@dx <- as.numeric(e1) %*% d_inv_a
    } else {
      e2@x <- e1 / x
      e2@dx <- - as.numeric(e1 / x^2) * dx
    }
    return(e2)
  }
)


#' Powers of 'dual'-class objects (dual-ANY)
#' @param e1 A "dual" object.
#' @param e2 "ANY" object.
setMethod(
  "^",
  signature(e1 = "dual", e2 = "ANY"),
  function(e1, e2) {
    if (!is_scalar(e2)) {
      stop("The exponent must be a scalar.")
    }

    if (e2 == 0) {
      e1@x <- e1@x^e2
      dx_dim <- dim(e1@dx)
      e1@dx <- zero_matrix0(dx_dim[1], dx_dim[2])
      return(e1)
    }

    px <- e1@x
    dx <- e1@dx
    e1@x <- px^e2
    e1@dx <- e2 * as.numeric(px^(e2 - 1)) * dx
    e1
  }
)


#' Powers of 'dual'-class objects (ANY-dual)
#' @param e1 "ANY" object.
#' @param e2 A "dual" object.
setMethod(
  "^",
  signature(e1 = "ANY", e2 = "dual"),
  function(e1, e2) {
    if (any(e1 <= 0)) {
      stop("When the exponent is a dual number, the base must be a matrix with positive entries to use e1 ^ e2 = exp(e2 * log(e1))")
    }
    if (!is_scalar(e2@x)) {
      stop("The exponent must be a scalar.")
    }
    px <- e1 ^ e2@x
    dx <- as.numeric(px) * (e2 * log(e1))@dx
    e2@x <- px
    e2@dx <- dx
    e2
  }
)


#' Matrix multiplication of 'dual'-class objects (dual-ANY)
#' @param x A "dual" object.
#' @param y "ANY" object.
setMethod(
  "%*%",
  signature(x = "dual", y = "ANY"),
  function(x, y) {
    x %*% empty_dual(y, ncol(x@dx))
  }
)


#' Matrix multiplication of 'dual'-class objects (ANY-dual)
#' @param x "ANY" object.
#' @param y A "dual" object.
setMethod(
  "%*%",
  signature(x = "ANY", y = "dual"),
  function(x, y) {
    empty_dual(x, ncol(y@dx)) %*% y
  }
)

#' Kronecker multiplication of 'dual'-class objects (dual-ANY)
#' @param X A "dual" object.
#' @param Y "ANY" object.
setMethod(
  "kronecker",
  signature(X = "dual", Y = "ANY"),
  function(X, Y) {
    X %x% empty_dual(Y, ncol(X@dx))
  }
)

#' Kronecker multiplication of 'dual'-class objects (ANY-dual)
#' @param X "ANY" object.
#' @param Y A "dual" object.
setMethod(
  "kronecker",
  signature(X = "ANY", Y = "dual"),
  function(X, Y) {
    empty_dual(X, ncol(Y@dx)) %x% Y
  }
)

empty_dual <- function(x, d) {
  new("dual", x = x, dx = zero_matrix0(length(x), d))
}
