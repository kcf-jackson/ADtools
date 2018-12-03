#' @include class_dual_def.R
NULL

#' Addition of 'dual'-class objects
#' @param e1 A "dual" object.
#' @param e2 A "dual" object.
setMethod("+",
          signature(e1 = "dual", e2 = "dual"),
          function(e1, e2) {
            x <- parent_of(e1) + parent_of(e2)
            new("dual", x = x, dx = d_sum(e1, e2), param = param_of(e1))
          }
)

#' Subtraction of 'dual'-class objects
#' @param e1 A "dual" object.
#' @param e2 A "dual" object.
setMethod("-",
          signature(e1 = "dual", e2 = "dual"),
          function(e1, e2) {
            x <- parent_of(e1) - parent_of(e2)
            new("dual", x = x, dx = d_minus(e1, e2), param = param_of(e1))
          }
)

#' (Element-wise) Multiplication of 'dual'-class objects
#' @param e1 A "dual" object.
#' @param e2 A "dual" object.
setMethod("*",
          signature(e1 = "dual", e2 = "dual"),
          function(e1, e2) {
            x <- parent_of(e1) * parent_of(e2)
            new("dual", x = x, dx = d_scalar_prod(e1, e2), param = param_of(e1))
          }
)

#' (Element-wise) Division of 'dual'-class objects
#' @param e1 A "dual" object.
#' @param e2 A "dual" object.
setMethod("/",
          signature(e1 = "dual", e2 = "dual"),
          function(e1, e2) {
            x <- parent_of(e1) / parent_of(e2)
            new("dual", x = x, dx = d_divide(e1, e2), param = param_of(e1))
          }
)

#' Matrix multiplication of 'dual'-class objects
#' @param x A "dual" object.
#' @param y A "dual" object.
setMethod("%*%",
          signature(x = "dual", y = "dual"),
          function(x, y) {
            res <- parent_of(x) %*% parent_of(y)
            new("dual", x = res, dx = d_matrix_prod(x, y), param = param_of(x))
          }
)

#' Kronecker multiplication of 'dual'-class objects
#' @param X A "dual" object.
#' @param Y A "dual" object.
setMethod("kronecker",
          signature(X = "dual", Y = "dual"),
          function(X, Y) {
            res <- parent_of(X) %x% parent_of(Y)
            new("dual", x = res, dx = d_kronecker(X, Y), param = param_of(X))
          }
)
