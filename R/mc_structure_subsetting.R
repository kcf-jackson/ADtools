#' @include class_dual_def.R
NULL


#' Extract parts of an object
#' @param x A "dual" object.
#' @param i integer; the row index.
#' @param j integer; the column index.
#' @param drop TRUE or FALSE; if TRUE, returns a vector when one
#' dimension of the matrix is 1.
#' @rdname index-subset
Extract.dual <- function(x, i, j, drop = FALSE) {
  x@dx <- d_subset(x, i, j)
  if (missing(i) && missing(j)) {
    x@x <- x@x[,, drop = drop]
  } else if (missing(i)) {
    x@x <- x@x[, j, drop = drop]
  } else if (missing(j)) {
    x@x <- x@x[i, , drop = drop]
  } else {
    x@x <- x@x[i, j, drop = drop]
  }
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

#' @rdname index-subset
setMethod("[", signature(x = "dual", i = "numeric", j = "missing"), Extract.dual)

#' @rdname index-subset
setMethod("[", signature(x = "dual", i = "missing", j = "numeric"), Extract.dual)

#' @rdname index-subset
setMethod("[", signature(x = "dual", i = "numeric", j = "numeric"), Extract.dual)


#' Return the First or Last Part of an Object
#' @rdname head.dual
#' @param x A "dual" object.
#' @param n A single integer.
head.dual <- function(x, n = 6) {
  assertthat::assert_that(length(n) == 1)

  px <- x@x
  L <- nrow(px)
  assertthat::assert_that((-L < n) && (n != 0))

  n <- ifelse(n < 0, max(L + n, 0), min(n, L))
  x[seq_len(n), , drop = FALSE]
}

#' @rdname head.dual
setMethod("head", signature(x = "dual"), head.dual)


#' Return the First or Last Part of an Object
#' @rdname tail.dual
#' @param x A "dual" object.
#' @param n A single integer.
tail.dual <- function(x, n = 6) {
  assertthat::assert_that(length(n) == 1)

  px <- x@x
  L <- nrow(px)
  assertthat::assert_that((-L < n) && (n != 0))

  n <- ifelse(n < 0, max(L + n, 0), min(L, n))
  ind <- seq.int(to = L, length.out = n)
  x[ind, , drop = FALSE]
}

#' @rdname tail.dual
setMethod("tail", signature(x = "dual"), tail.dual)
