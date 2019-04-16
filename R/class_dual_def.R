check_dual <- function(object) {
  x <- object@x
  dx <- object@dx
  x_check <- is.array(x) || is.numeric(x) || is_matrix(x)
  dim_check <- nrow(dx) == length(x)   # match dimension
  x_check && is_matrix(dx) && dim_check
}
is_matrix <- function(x) {
  (class(x) == "matrix") || is_sparse_matrix(x)
}
is_sparse_matrix <- function(x) {
  attr_cls <- attr(class(x), "package")
  !is.null(attr_cls) && (attr_cls == "Matrix")
}

#' S4 class "dual"
#'
#' @description This class attaches a dual component to a number / an array.
#' @slot x scalar, vector or matrix; also accepts any matrix classes from the "Matrix" package.
#' @slot dx matrix; also accepts any matrix classes from the "Matrix" package.
#' @slot param a named list, containing the column indices each variable occupies.
#' @import methods
#' @note Users should not construct the object directly, instead, use the
#' constructor helper `dual` provided.
setClass(
  "dual",
  representation(x = "ANY", dx = "ANY", param = "list"),
  validity = check_dual
)


#' Dual number extractor
#' @param x A dual number
#' @export
parent_of <- function(x) { x@x }

#' Dual number extractor
#' @param x A dual number
#' @export
deriv_of <- function(x) { x@dx }

#' Dual number extractor
#' @param x A dual number
#' @export
param_of <- function(x) { x@param }
