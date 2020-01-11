#' Sum of matrix elements
#' @param x A "dual" object.
setMethod(
  "sum",
  signature(x = "dual"),
  function(x) {
    x@dx <- t(colSums(x@dx))
    x@x <- sum(x@x)
    x
  }
)


#' Row sum of a matrix.
#' @param x A "dual" object.
setMethod(
  "rowSums",
  signature(x = "dual"),
  function(x) {
    x@dx <- d_rowSums(x)
    x@x <- rowSums(x@x)
    x
  }
)

d_rowSums <- function(x) {
  px <- x@x
  px_len <- length(px)
  px_nr <- nrow(px)

  dx <- x@dx
  purrr::map(
    1:px_nr,
    ~ colSums(dx[seq(.x, px_len, px_nr), ])
  ) %>%
    do.call(rbind, .)
}


#' Column sum of a matrix.
#' @param x A "dual" object.
setMethod(
  "colSums",
  signature(x = "dual"),
  function(x) {
    x@dx <- d_colSums(x)
    x@x <- colSums(x@x)
    x
  }
)

d_colSums <- function(x) {
  px <- x@x
  px_len <- length(px)
  px_nr <- nrow(px)

  dx <- x@dx
  purrr::map2(
    seq(1, px_len, px_nr), seq(px_nr, px_len, px_nr),
    ~ colSums(dx[.x:.y, ])
  ) %>%
    do.call(rbind, .)
}


#' Trace of a matrix
#' @param x A square matrix
#' @export
tr <- function(x) {
  if (nrow(x) != ncol(x)) {
    stop("Input must be a square matrix")
  }
  sum(diag(x))
}

#' Trace of a matrix
#' @param x A "dual" object.
setMethod(
  "tr",
  signature(x = "dual"),
  function(x) {
    sum(diag(x))
  }
)

#' Mean of vector or matrix elements
#' @param x A "dual" object.
setMethod(
  "mean",
  signature(x = "dual"),
  function(x) {
    x@dx <- t(colMeans(x@dx))
    x@x <- mean(x@x)
    x
  }
)
