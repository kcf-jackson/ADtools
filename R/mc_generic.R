#' @include class_dual_def.R
NULL

#' Length of an Object
#' @param x A "dual" object.
setMethod("length", signature(x = "dual"), function(x) length(x@x))


#' Dimension of an Object
#' @param x A "dual" object.
setMethod("dim", signature(x = "dual"), function(x) dim(x@x))


#' Number of rows
#' @param x A "dual" object.
setMethod("nrow", signature(x = "dual"), function(x) nrow(x@x))


#' Number of columns
#' @param x A "dual" object.
setMethod("ncol", signature(x = "dual"), function(x) ncol(x@x))


#' Rounding of Numbers
#' @param x A "dual" object.
#' @param digits An integer indicating the number of decimal places.
#' @note The function 'round' does not have a derivative over the real line.
#' The derivative will be kept unchanged. The reason of not dropping it is
#' that sometimes one may need to round a matrix to correct floating-point
#' errors. This is often used when an apparent symmetric matrix does not pass
#' the symmetry test.
#' @export
round.dual <- function(x, digits = 0) {
  x@x <- round(x@x, digits = digits)
  return(x)
}

#' Rounding of Numbers
#' @param x A "dual" object.
#' @param digits integer indicating the number of decimal places.
setMethod("round", signature(x = "dual", digits = "integer"), round.dual)
