# Partition a matrix by rows into a list of matrices
partition_by_rows <- function(Z, nrow_block, nblock) {
  if (missing(nrow_block)) {
    nrow_block <- nrow(Z) / nblock
  }
  assertthat::assert_that(
    nrow(Z) %% nrow_block == 0,
    msg = "The number of rows is not divisible by the block size."
  )

  purrr::map(
    seq(1, nrow(Z), nrow_block),
    ~Z[.x:(.x + nrow_block - 1), ]
  )
}

`%L-circledast%` <- function(x, y) {
  B <- x
  Z <- partition_by_rows(y, nblock = ncol(x))
  mapreduce(1:nrow(B), function(k) {
    purrr::reduce(purrr::map2(B[k, ], Z, `*`), `+`)
  }, rbind2)
}

`%L-boxdot%` <- function(x, y) {
  A <- x
  V <- partition_by_rows(y, nrow_block = ncol(x))
  mapreduce(V, ~A %*% .x, rbind2)
}

#' Post-multiplying a kronecker product
#' @description Compute (B \%x\% A) Z
#' @param B numeric matrix.
#' @param A numeric matrix.
#' @param Z numeric matrix.
#' @export
BxAZ <- function(B, A, Z) {
  A %L-boxdot% (B %L-circledast% Z)
}

# optim_BxAZ <- function(B, A, Z) {
#   Z <- partition_by_rows(Z, nblock = ncol(B))
#   1:nrow(B) %>%
#     purrr::map(function(k) {
#       A %*% purrr::reduce(purrr::map2(B[k, ], Z, `*`), `+`)
#     }) %>%
#     do.call(rbind, .)
# }


# Partition a matrix by columns into a list of matrices
partition_by_columns <- function(Z, ncol_block, nblock) {
  if (missing(ncol_block)) {
    ncol_block <- ncol(Z) / nblock
  }
  assertthat::assert_that(
    ncol(Z) %% ncol_block == 0,
    msg = "The number of columns is not divisible by the block size."
  )

  purrr::map(
    seq(1, ncol(Z), ncol_block),
    ~Z[, .x:(.x + ncol_block - 1)]
  )
}

`%R-circledast%` <- function(x, y) {
  B <- y
  Z <- partition_by_columns(x, nblock = nrow(y))
  do.call(cbind, purrr::map(1:ncol(B), function(k) {
    purrr::reduce(purrr::map2(B[, k], Z, `*`), `+`)
  }))
}

`%R-boxdot%` <- function(x, y) {
  A <- y
  V <- partition_by_columns(x, ncol_block = nrow(y))
  do.call(cbind, purrr::map(V, ~.x %*% A))
}


#' Pre-multiplying a kronecker product
#' @description Compute X (B \%x\% A)
#' @param X numeric matrix.
#' @param B numeric matrix.
#' @param A numeric matrix.
#' @export
XBxA <- function(X, B, A) {
  (X %R-circledast% B) %R-boxdot% A
}


#' Compute (B \%x\% I) D
#' @param B numeric matrix.
#' @param D numeric matrix.
#' @export
BxID <- function(B, D) {
  B %L-circledast% D
}

#' Compute (I \%x\% C) D
#' @param C numeric matrix.
#' @param D numeric matrix.
#' @export
IxCD <- function(C, D) {
  C %L-boxdot% D
}

#' Compute A (B \%x\% I)
#' @param A numeric matrix.
#' @param B numeric matrix.
#' @export
ABxI <- function(A, B) {
  A %R-circledast% B
}

#' Compute A (I \%x\% C)
#' @param A numeric matrix.
#' @param C numeric matrix.
#' @export
AIxC <- function(A, C) {
  A %R-boxdot% C
}
