#' @include class_dual_def.R
NULL

# This file contains functions to handle cross-classes arithmetic.
# Whenever an operation is performed on a dual-number and a scalar / matrix,
# it is assumed to be performed component-wise. Furthermore, all
# non-dual objects are assumed to have zero matrices as Jacobian.

#' Addition of 'dual'-class objects
#' @param e1 A "dual" object.
#' @param e2 "ANY" object.
setMethod("+",
          signature(e1 = "dual", e2 = "ANY"),
          function(e1, e2) { e1@x <- parent_of(e1) + e2; e1 })

#' Addition of 'dual'-class objects
#' @param e1 "ANY" object.
#' @param e2 A "dual" object.
setMethod("+",
          signature(e1 = "ANY", e2 = "dual"),
          function(e1, e2) { e2@x <- e1 + parent_of(e2); e2 })

#' Subtraction of 'dual'-class objects
#' @param e1 A "dual" object.
#' @param e2 "ANY" object.
setMethod("-",
          signature(e1 = "dual", e2 = "ANY"),
          function(e1, e2) { e1@x <- parent_of(e1) - e2; e1})

#' Subtraction of 'dual'-class objects
#' @param e1 A "dual" object.
#' @param e2 "ANY" object.
setMethod("-",
          signature(e1 = "dual", e2 = "missing"),
          function(e1, e2) {
            e1@x <- - parent_of(e1)
            e1@dx <- - deriv_of(e1)
            e1
          }
)

#' Subtraction of 'dual'-class objects
#' @param e1 "ANY" object.
#' @param e2 A "dual" object.
setMethod("-",
          signature(e1 = "ANY", e2 = "dual"),
          function(e1, e2) {
            e2@x <- e1 - parent_of(e2)
            e2@dx <- - deriv_of(e2)
            e2
          }
)

#' (Element-wise) Multiplication of 'dual'-class objects
#' @param e1 A "dual" object.
#' @param e2 "ANY" object.
setMethod("*",
          signature(e1 = "dual", e2 = "ANY"),
          function(e1, e2) {
            if (is_scalar(e2)) {
              e1@x <- parent_of(e1) * e2
              e1@dx <- deriv_of(e1) * e2
            } else {
              # matrix case
              e1@x <- parent_of(e1) * e2
              e1@dx <- deriv_of(e1) * as.numeric(e2)
            }
            e1
          }
)

#' (Element-wise) Multiplication of 'dual'-class objects
#' @param e1 "ANY" object.
#' @param e2 A "dual" object.
setMethod("*",
          signature(e1 = "ANY", e2 = "dual"),
          function(e1, e2) {
            if (is_scalar(e1)) {
              e2@x <- e1 * parent_of(e2)
              e2@dx <- e1 * deriv_of(e2)
            } else {
              # matrix case
              e2@x <- e1 * parent_of(e2)
              e2@dx <- as.numeric(e1) * deriv_of(e2)
            }
            e2
          }
)

#' (Element-wise) Division of 'dual'-class objects
#' @param e1 A "dual" object.
#' @param e2 "ANY" object.
setMethod("/",
          signature(e1 = "dual", e2 = "ANY"),
          function(e1, e2) {
            if (is_scalar(e2)) {
              e1@x <- parent_of(e1) / e2
              e1@dx <- deriv_of(e1) / e2
            } else {
              e1@x <- parent_of(e1) / e2
              e1@dx <- deriv_of(e1) / as.numeric(e2)
            }
            e1
          }
)

#' (Element-wise) Division of 'dual'-class objects
#' @param e1 "ANY" object.
#' @param e2 A "dual" object.
setMethod("/",
          signature(e1 = "ANY", e2 = "dual"),
          function(e1, e2) {
            x <- e2@x
            dx <- e2@dx
            for (j in 1:ncol(x)) {
              for (i in 1:nrow(x)) {
                deriv_rind <- i + (j-1) * nrow(x)
                inv_a <- 1 / x[i, j, drop = F]
                d_inv_a <- - as.numeric(inv_a)^2 * dx[deriv_rind, , drop = F]
                k <- ifelse(is_scalar(e1), e1, e1[i, j])
                e2@x[i, j] <- k * inv_a
                e2@dx[deriv_rind, ] <- k * d_inv_a
              }
            }
            return(e2)
          }
)

#' Matrix multiplication of 'dual'-class objects
#' @param x A "dual" object.
#' @param y "ANY" object.
setMethod("%*%",
          signature(x = "dual", y = "ANY"),
          function(x, y) { x %*% dual(y, ncol(x@dx), -1) })

#' Matrix multiplication of 'dual'-class objects
#' @param x "ANY" object.
#' @param y A "dual" object.
setMethod("%*%",
          signature(x = "ANY", y = "dual"),
          function(x, y) { dual(x, ncol(y@dx), -1) %*% y })

#' Kronecker multiplication of 'dual'-class objects
#' @param X A "dual" object.
#' @param Y "ANY" object.
setMethod("kronecker",
          signature(X = "dual", Y = "ANY"),
          function(X, Y) { X %x% dual(Y, ncol(X@dx), -1) })

#' Kronecker multiplication of 'dual'-class objects
#' @param X "ANY" object.
#' @param Y A "dual" object.
setMethod("kronecker",
          signature(X = "ANY", Y = "dual"),
          function(X, Y) { dual(X, ncol(Y@dx), -1) %x% Y })
