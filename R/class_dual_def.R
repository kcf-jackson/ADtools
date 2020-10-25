check_dual <- function(object) {
  x <- object@x
  dx <- object@dx

  x_check <- is_scalar(x) || is_vector(x) || is_matrix(x)
  dx_check <- is_matrix(dx)
  dim_check <- nrow(dx) == length(x)   # match dimension

  x_check && dx_check && dim_check
}

is_scalar <- function(x) {
  is.numeric(x) && (length(x) == 1)
}

is_vector <- function(x) {
  is.numeric(x) && is.vector(x) && (length(x) > 1)
}

is_matrix <- function(x) {
  is.matrix(x) || is_sparse_matrix(x)
}

is_sparse_matrix <- function(x) {
  attr_cls <- attr(class(x), "package")
  !is.null(attr_cls) && (attr_cls == "Matrix")
}


#' S4 class "dual"
#'
#' @description This class attaches a dual component to a number / matrix.
#' @slot x scalar, vector or matrix; also accepts any matrix classes from the "Matrix" package.
#' @slot dx matrix; also accepts any matrix classes from the "Matrix" package.
#' @slot param a named list, containing the column indices each variable occupies.
#'
#' @import methods
#'
#' @note Users should not construct the object directly, instead, use the
#' constructor helper \link{dual} provided.
setClass(
  "dual",
  representation(x = "ANY", dx = "ANY"),
  validity = check_dual
)
