#' @include class_dual_def.R
NULL


subset_fun <- function(x, i, j, drop = F) {
  new("dual",
      x = parent_of(x)[i, j, drop = drop],
      dx = d_subset(x, i, j),
      param = x@param)
}
setMethod("[", signature(x = "dual", i = "numeric", j = "missing", drop = "logical"), subset_fun)
setMethod("[", signature(x = "dual", i = "missing", j = "numeric", drop = "logical"), subset_fun)
setMethod("[", signature(x = "dual", i = "numeric", j = "numeric", drop = "logical"), subset_fun)


#' Return the First or Last Part of an Object
#' @param x A "dual" object.
#' @param n A single integer.
setMethod("head",
  signature(x = "dual"),
  function(x, n = 6) {
    px <- parent_of(x)
    if (is.vector(px)) return(x[seq(n)])
    x[seq(n), , drop = F]
  }
)


#' Return the First or Last Part of an Object
#' @param x A "dual" object.
#' @param n A single integer.
setMethod("tail",
  signature(x = "dual"),
  function(x, n = 6) {
    px <- parent_of(x)
    if (is.vector(px)) {
      L <- length(px)
      ind <- L - max(n-1, 0):0
      return(x[ind])
    }
    L <- nrow(px)
    ind <- L - max(n-1, 0):0
    x[ind, , drop = F]
  }
)
