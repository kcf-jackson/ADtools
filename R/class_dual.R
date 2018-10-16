#' @include class_matrix_list_def.R
#' @include class_dual_def.R
NULL

#=========================================================================
# Arithmetic for dual numbers
#-------------------------------------------------------------------------
#' (Scalar) Addition of 'dual'-class objects
#' @param e1 A "dual" object.
#' @param e2 A "dual" object.
setMethod("+",
  signature(e1 = "dual", e2 = "dual"),
  function(e1, e2) {
    x <- parent_of(e1) + parent_of(e2)
    new("dual", x = x, dx = d_sum(e1, e2))
  }
)

#' (Scalar) Subtraction of 'dual'-class objects
#' @param e1 A "dual" object.
#' @param e2 A "dual" object.
setMethod("-",
  signature(e1 = "dual", e2 = "dual"),
  function(e1, e2) {
    x <- parent_of(e1) - parent_of(e2)
    new("dual", x = x, dx = d_minus(e1, e2))
  }
)

#' (Scalar) Multiplication of 'dual'-class objects
#' @param e1 A "dual" object.
#' @param e2 A "dual" object.
setMethod("*",
  signature(e1 = "dual", e2 = "dual"),
  function(e1, e2) {
    x <- parent_of(e1) * parent_of(e2)
    new("dual", x = x, dx = d_scalar_prod(e1, e2))
  }
)

#' (Scalar) Division of 'dual'-class objects
#' @param e1 A "dual" object.
#' @param e2 A "dual" object.
setMethod("/",
  signature(e1 = "dual", e2 = "dual"),
  function(e1, e2) {
    x <- parent_of(e1) / parent_of(e2)
    new("dual", x = x, dx = d_divide(e1, e2))
  }
)

#' Matrix multiplication of 'dual'-class objects
#' @param x A "dual" object.
#' @param y A "dual" object.
setMethod("%*%",
  signature(x = "dual", y = "dual"),
  function(x, y) {
    res <- parent_of(x) %*% parent_of(y)
    new("dual", x = res, dx = d_matrix_prod(x, y))
  }
)

#' Matrix multiplication of 'dual'-class objects
#' @param x A "dual" object.
#' @param y A "dual" object.
setMethod("%x%",
  signature(X = "dual", Y = "dual"),
  function(X, Y) {
    res <- parent_of(X) %x% parent_of(Y)
    new("dual", x = res, dx = d_kronecker(X, Y))
  }
)

#' Inverse of 'dual'-class objects
#' @param a A "dual" object.
setMethod("solve",
  signature(a = "dual"),
  function(a) {
    inv_a <- solve(parent_of(a))
    new("dual", x = inv_a, dx = d_solve(a, inv_a))
  }
)

#' Transpose of 'dual'-class objects
#' @param a A "dual" object.
setMethod("t",
  signature(x = "dual"),
  function(x) {
    new("dual", x = t(parent_of(x)), dx = d_transpose(x))
  }
)

#' Crossproduct of 'dual'-class objects
#' @param a A "dual" object.
setMethod("tcrossprod",
  signature(x = "dual", y = "missing"),
  function(x, y) {
    new("dual", x = tcrossprod(parent_of(x)), dx = d_XXT(x))
  }
)

#' Crossproduct of 'dual'-class objects
#' @param a A "dual" object.
setMethod("tcrossprod",
  signature(x = "dual", y = "dual"),
  function(x, y) { x %*% t(y) }
)


#' Crossproduct of 'dual'-class objects
#' @param a A "dual" object.
setMethod("crossprod",
  signature(x = "dual", y = "missing"),
  function(x, y) { t(x) %*% x }
)

#' Crossproduct of 'dual'-class objects
#' @param a A "dual" object.
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
#' @param a A "dual" object.
#' @note The Cholesky decomposition used in this function returns a
#' lower triangular matrix.
setMethod("chol0",
  signature(x = "dual"),
  function(x) {
    L <- chol0(parent_of(x))
    new("dual", x = L, dx = d_chol(L, x))
  }
)

#' #' Combine 'dual'-class objects by Columns
#' #' @param a A "dual" object.
#' #' @note The Cholesky decomposition used in this function returns a
#' #' lower triangular matrix.
#' setMethod("cbind",
#'   signature(x = "dual"),
#'   function(x, y) {
#'     cbind(parent_of(x), parent_of(y))
#'
#'     new("dual", x = L, dx = d_chol(L, x))
#'   }
#' )
#'

#' # CHECK MADNESS REPO
#' #' Subsetting
#' setMethod("[",
#'   signature(x = "dual", i = "integer", j = "missing"),
#'   function(x, i, j, ...) {
#'     new("dual", x = parent_of(x[i,j,...]), dx = d_subset(x, i, j, ...))
#'   }
#' )
#' setMethod("[",
#'   signature(x = "dual", i = "missing", j = "integer"),
#'   function(x, i, j, ...) {
#'     new("dual", x = parent_of(x[i,j,...]), dx = d_subset(x, i, j, ...))
#'   }
#' )
#'
#' setMethod("length",
#'   signature(x = "dual"),
#'   function(x) { length(parent_of(x)) }
#' )


#' #'
#' setMethod("rnorm",
#'   signature(n = "int", mean = "dual", sd = "dual"),
#'   function(n, mean, sd) {
#'     z <- rnorm(n)
#'     res <- parent_of(mean) + parent_of(sd) * z
#'     new("dual", x = res, dx = d_rnorm(mean, sd, z))
#'   }
#' )
#'
#' #'
#' setMethod("rgamma",
#'   signature(n, shape, scale),
#'   function(n, shape, scale) {
#'
#'   }
#' )
#'
#' #'
#' setMethod("rchisq",
#'   signature(n, df),
#'   function(n, df) {
#'
#'   }
#' )
#'
#' #'
#' setMethod("rexp",
#'   signature(n = "int", rate = "dual"),
#'   function(n, rate) {
#'     rgamma(n, 1)
#'   }
#' )
