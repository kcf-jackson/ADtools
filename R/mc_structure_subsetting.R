#' @include class_dual_def.R
NULL


subset_fun <- function(x, i, j, drop = F) {
  x@dx <- d_subset(x, i, j)
  x@x <- x@x[i, j, drop = drop]
  x
}

d_subset <- function(a, i, j) {
  x <- a@x
  dx <- a@dx
  nr <- nrow(x)
  nc <- ncol(x)
  if (missing(i) && missing(j)) {
    stop("At least one of index i and index j should be present.") # nocov
  } else if (missing(i) && !missing(j)) {
    ind <- mapreduce(j, ~ seq(nr) + nr * (.x - 1), c)
  } else if (!missing(i) && missing(j)) {
    ind <- mapreduce(i, ~ .x + nr * (seq(nc) - 1), c)
  } else {
    ind <- map2reduce(i, j, ~ .x + nr * (.y - 1), c)
  }
  dx[sort(ind), , drop = FALSE]
}

#' Extract parts of an object
#' @param x A "dual" object.
#' @param i integer; the row index.
#' @param j integer; the column index.
#' @param drop T or F; if T, returns a vector when one dimension of the matrix is 1.
#' @rdname index-subset
setMethod("[", signature(x = "dual", i = "numeric", j = "missing"), subset_fun)

#' @rdname index-subset
setMethod("[", signature(x = "dual", i = "missing", j = "numeric"), subset_fun)

#' @rdname index-subset
setMethod("[", signature(x = "dual", i = "numeric", j = "numeric"), subset_fun)


#' Return the First or Last Part of an Object
#' @param x A "dual" object.
#' @param n A single integer.
setMethod("head",
  signature(x = "dual"),
  function(x, n = 6) {
    assertthat::assert_that(length(n) == 1)

    px <- x@x
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

    px <- x@x
    L <- nrow(px)
    assertthat::assert_that((-L < n) && (n != 0))

    n <- ifelse(n < 0, max(L + n, 0), min(L, n))
    ind <- seq.int(to = L, length.out = n)
    x[ind, , drop = F]
  }
)
