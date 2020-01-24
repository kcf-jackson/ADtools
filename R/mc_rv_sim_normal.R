#' @include class_dual_def.R
NULL

# ===============================================================
# Gaussian distribution
# ===============================================================
#' Simulate univariate normal random variates
#' @name rnorm0
#' @param n Positive integer; the number of samples.
#' @param mean A dual number or a numeric vector; the mean of the normal distribution.
#' @param sd A dual number or a numeric vector; the standard deviation of the normal distribution.
#' @export
rnorm0 <- function(n, mean = 0, sd = 1) {
  rnorm(n, mean, sd)
}

rnorm0_dual <- function(n, mean, sd) {
  len_mean <- length(mean)
  assertthat::assert_that(len_mean == 1 || len_mean == n)
  len_sd <- length(sd)
  assertthat::assert_that(len_sd == 1 || len_sd == n)
  mean + as.matrix(sd * rnorm(n))
}

#' @rdname rnorm0
setMethod(
  "rnorm0",
  signature(n = "numeric", mean = "dual", sd = "dual"),
  rnorm0_dual
)

#' @rdname rnorm0
setMethod(
  "rnorm0",
  signature(n = "numeric", mean = "ANY", sd = "dual"),
  rnorm0_dual
)

#' @rdname rnorm0
setMethod(
  "rnorm0",
  signature(n = "numeric", mean = "dual", sd = "ANY"),
  rnorm0_dual
)


#' Simulate multivariate normal random variates
#' @param n Positive integer; the number of samples.
#' @param mean mean vector of the normal distribution.
#' @param sigma covariance matrix of the normal distribution.
#' @return A numeric matrix, where every column is a sample.
#' @export
rmvnorm0 <- function(n, mean, sigma) {
  Z <- t(mvtnorm::rmvnorm(n, numeric(length(mean))))
  t(mean + chol0(sigma) %*% Z)
}

#' Simulate multivariate normal random variates
#' @rdname rmvnorm0_dual
#' @param n Positive integer; the number of samples.
#' @param mean A numeric vector or a dual number; the mean of the normal distribution.
#' @param sigma A numeric matrix or a dual number; the standard deviation of the normal distribution.
rmvnorm0_dual <- function(n, mean, sigma) {
  U <- t(chol0(sigma))
  Z <- mvtnorm::rmvnorm(n, numeric(length(mean)))
  add_vector_to_matrix_row(mean, Z %*% U)
}


#' @rdname rmvnorm0_dual
setMethod(
  "rmvnorm0",
  signature(n = "numeric", mean = "dual", sigma = "dual"),
  rmvnorm0_dual
)

#' @rdname rmvnorm0_dual
setMethod(
  "rmvnorm0",
  signature(n = "numeric", mean = "dual", sigma = "ANY"),
  rmvnorm0_dual
)

#' @rdname rmvnorm0_dual
setMethod(
  "rmvnorm0",
  signature(n = "numeric", mean = "ANY", sigma = "dual"),
  rmvnorm0_dual
)
