#' @include class_dual_def.R
NULL


subset_fun <- function(x, i, j, drop = F) {
  new("dual",
      x = parent_of(x)[i, j, drop = drop],
      dx = d_subset(x, i, j),
      param = param_of(x))
}
setMethod("[", signature(x = "dual", i = "numeric", j = "missing"), subset_fun)
setMethod("[", signature(x = "dual", i = "missing", j = "numeric"), subset_fun)
setMethod("[", signature(x = "dual", i = "numeric", j = "numeric"), subset_fun)


#' Return the First or Last Part of an Object
#' @param x A "dual" object.
#' @param n A single integer.
setMethod("head",
  signature(x = "dual"),
  function(x, n = 6) {
    assertthat::assert_that(length(n) == 1)

    px <- parent_of(x)
    # if (is.vector(px)) {
    #   L <- length(px)
    #   assertthat::assert_that((-L < n) && (n != 0))
    #
    #   n <- ifelse(n < 0, max(length(x) + n, 0), min(n, L))
    #   return(x[seq_len(n)])
    # }
    # matrix case
    L <- nrow(px)
    assertthat::assert_that((-L < n) && (n != 0))

    n <- ifelse(n < 0, max(L + n, 0), min(n, L))
    x[seq_len(n), , drop = F]
  }
)


#' Return the First or Last Part of an Object
#' @param x A "dual" object.
#' @param n A single integer.
setMethod("tail",
  signature(x = "dual"),
  function(x, n = 6) {
    assertthat::assert_that(n != 0)

    px <- parent_of(x)
    # if (is.vector(px)) {
    #   L <- length(px)
    #   assertthat::assert_that((-L < n) && (n != 0))
    #
    #   n <- ifelse(n < 0, max(L + n, 0), min(L, n))
    #   ind <- seq.int(to = L, length.out = n)
    #   return(x[ind])
    # }
    # matrix case
    L <- nrow(px)
    assertthat::assert_that((-L < n) && (n != 0))

    n <- ifelse(n < 0, max(L + n, 0), min(L, n))
    ind <- seq.int(to = L, length.out = n)
    x[ind, , drop = F]
  }
)
