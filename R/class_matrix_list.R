#' @include class_matrix_list_def.R
NULL

#=========================================================================
# Methods
#-------------------------------------------------------------------------
setMethod_by_operation <- function(op) {
  op_fun <- match.fun(op)
  setMethod(op,
    signature(e1 = "matrix_list", e2 = "matrix_list"),
    function(e1, e2) {
      m <- purrr::map2(element_of(e1), element_of(e2), op_fun)
      new("matrix_list", matrices = m)
    }
  )
  setMethod(op,
    signature(e1 = "ANY", e2 = "matrix_list"),
    function(e1, e2) {
      m <- purrr::map2(list(e1), element_of(e2), op_fun)
      new("matrix_list", matrices = m)
    }
  )
  setMethod(op,
    signature(e1 = "matrix_list", e2 = "ANY"),
    function(e1, e2) {
      m <- purrr::map2(element_of(e1), list(e2), op_fun)
      new("matrix_list", matrices = m)
    }
  )
}

#' Addition of 'matrix_list'-class objects
#' @param e1 A "matrix_list" object.
#' @param e2 A "matrix_list" object.
setMethod_by_operation("+")

#' Subtraction of 'matrix_list'-class objects
#' @param e1 A "matrix_list" object.
#' @param e2 A "matrix_list" object.
setMethod_by_operation("-")

#' (Element-wise) Multiplication of 'matrix_list'-class objects
#' @param e1 A "matrix_list" object.
#' @param e2 A "matrix_list" object.
setMethod_by_operation("*")

#' (Element-wise) Division of 'matrix_list'-class objects
#' @param e1 A "matrix_list" object.
#' @param e2 A "matrix_list" object.
setMethod_by_operation("/")

# Handling exception
#' Negation of 'matrix_list'-class objects
#' @param e1 A "matrix_list" object.
#' @param e2 Empty
setMethod("-",
  signature(e1 = "matrix_list", e2 = "missing"),
  function(e1, e2) {
    new("matrix_list", matrices = purrr::map(element_of(e1), `-`))
  }
)

#-------------------------------------------------------------------------
#' Matrix multiplication of 'matrix_list'-class objects
#' @param x A "matrix_list" object.
#' @param y A "matrix_list" object.
# %*% has variable names different to the others, so it needs separate handling.
setMethod("%*%",
  signature(x = "matrix_list", y = "matrix_list"),
  function(x, y) {
    m <- purrr::map2(element_of(x), element_of(y), `%*%`)
    new("matrix_list", matrices = m)
  }
)

#' Matrix multiplication of 'matrix_list'-class objects
#' @param x Any object other than "matrix_list" object.
#' @param y A "matrix_list" object.
setMethod("%*%",
  signature(x = "ANY", y = "matrix_list"),
  function(x, y) {
    m <- purrr::map2(list(x), element_of(y), `%*%`)
    new("matrix_list", matrices = m)
  }
)

#' Matrix multiplication of 'matrix_list'-class objects
#' @param x A "matrix_list" object.
#' @param y Any object other than "matrix_list" object.
setMethod("%*%",
  signature(x = "matrix_list", y = "ANY"),
  function(x, y) {
    m <- purrr::map2(element_of(x), list(y), `%*%`)
    new("matrix_list", matrices = m)
  }
)

#-------------------------------------------------------------------------
#' Matrix kronecker multiplication of 'matrix_list'-class objects
#' @param X A "matrix_list" object.
#' @param Y A "matrix_list" object.
# %x% has variable names different to the others, so it needs separate handling.
setMethod("%x%",
  signature(X = "matrix_list", Y = "matrix_list"),
  function(X, Y) {
    m <- purrr::map2(element_of(X), element_of(Y), `%x%`)
    new("matrix_list", matrices = m)
  }
)

#' Matrix kronecker multiplication of 'matrix_list'-class objects
#' @param X Any object other than "matrix_list" object.
#' @param Y A "matrix_list" object.
setMethod("%x%",
  signature(X = "ANY", Y = "matrix_list"),
  function(X, Y) {
    m <- purrr::map2(list(X), element_of(Y), `%x%`)
    new("matrix_list", matrices = m)
  }
)

#' Matrix kronecker multiplication of 'matrix_list'-class objects
#' @param X A "matrix_list" object.
#' @param Y Any object other than "matrix_list" object.
setMethod("%x%",
  signature(X = "matrix_list", Y = "ANY"),
  function(X, Y) {
    m <- purrr::map2(element_of(X), list(Y), `%x%`)
    new("matrix_list", matrices = m)
  }
)

#-------------------------------------------------------------------------
#' Combine two 'matrix_list'-class objects by Columns
#' @param x A "matrix_list" object.
#' @param y A "matrix_list" object.
setMethod("cbind2",
  signature(x = "matrix_list", y = "matrix_list"),
  function(x, y) {
    m <- purrr::map2(element_of(X), element_of(Y), `cbind2`)
    new("matrix_list", matrices = m)
  }
)

#' Combine two 'matrix_list'-class objects by Columns
#' @param x Any object other than "matrix_list" object.
#' @param y A "matrix_list" object.
setMethod("cbind2",
  signature(x = "ANY", y = "matrix_list"),
  function(x, y) {
    m <- purrr::map2(list(x), element_of(y), `cbind2`)
    new("matrix_list", matrices = m)
  }
)

#' Combine two 'matrix_list'-class objects by Columns
#' @param x A "matrix_list" object.
#' @param y Any object other than "matrix_list" object.
setMethod("cbind2",
  signature(x = "matrix_list", y = "ANY"),
  function(x, y) {
    m <- purrr::map2(element_of(x), list(y), `cbind2`)
    new("matrix_list", matrices = m)
  }
)
