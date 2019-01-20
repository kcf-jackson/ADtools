#' @include class_dual_def.R
NULL

#===============================================================
# Gaussian distribution
#===============================================================
#' Simulate univariate normal random variates
#' @inherit stats::rnorm
#' @export
rnorm0 <- function(n, mean = 0, sd = 1) {
  rnorm(n, mean, sd)
}


#' Simulate univariate normal random variates
#' @param n Positive integer; the number of samples.
#' @param mean A dual number; the mean of the normal distribution.
#' @param sd A dual number; the standard deviation of the normal distribution.
setMethod("rnorm0",
  signature(n = "numeric", mean = "dual", sd = "dual"),
  function(n, mean, sd) {
    len_mean <- length(mean)
    assertthat::assert_that(len_mean == 1 || len_mean == n)

    len_sd <- length(sd)
    assertthat::assert_that(len_sd == 1 || len_sd == n)

    mean + sd * rnorm(n)
  }
)

#' Simulate multivariate normal random variates
#' @param n Positive integer; the number of samples.
#' @param mean mean vector of the normal distribution.
#' @param sigma covariance matrix of the normal distribution.
# #' @param ... Other parameters to be passed to `mvtnorm::rmvnorm`
#' @export
rmvnorm0 <- function(n, mean, sigma) {
  mvtnorm::rmvnorm(n, mean, sigma)
}

#' Simulate multivariate normal random variates
#' @param n Positive integer; the number of samples.
#' @param mean A dual number; the mean of the normal distribution.
#' @param sigma A dual number; the standard deviation of the normal distribution.
setMethod("rmvnorm0",
  signature(n = "numeric", mean = "dual", sigma = "dual"),
  function(n, mean, sigma) {
    mapreduce(
      seq(n),
      ~mean + chol0(sigma) %*%
        t(mvtnorm::rmvnorm(1, numeric(length(mean)))),
      cbind2
    )
  }
)
