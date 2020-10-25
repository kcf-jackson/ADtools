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
  px_nc <- ncol(px)
  dx <- x@dx

  m0 <- mapreduce(1:px_nr, ~cbind(.x, seq(.x, px_len, px_nr)), rbind)
  Matrix::sparseMatrix(m0[,1], m0[,2], x = 1) %*% dx
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
  px_nc <- ncol(px)

  dx <- x@dx
  s <- seq(1, px_len, px_nr)
  purrr::map2(s, s + px_nr - 1, ~colSums(dx[.x:.y, ])) %>%
    do.call(rbind, .)
}


#' Trace of a matrix
#' 
#' @param x A square numeric matrix
#' 
#' @examples 
#' A <- randn(2, 2)
#' tr(A)
#' 
#' B <- randn(3, 3)
#' tr(B)
#' 
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

#' Row mean of a matrix.
#' @param x A "dual" object.
setMethod(
  "rowMeans",
  signature(x = "dual"),
  function(x) {
    x@dx <- d_rowSums(x) / NCOL(x@x)
    x@x <- rowMeans(x@x)
    x
  }
)


#' Column mean of a matrix.
#' @param x A "dual" object.
setMethod(
  "colMeans",
  signature(x = "dual"),
  function(x) {
    x@dx <- d_colSums(x) / NROW(x@x)
    x@x <- colMeans(x@x)
    x
  }
)
