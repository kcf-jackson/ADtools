#' The density of the multivariate t distribution
#' @param x vector or matrix of quantiles. If x is a matrix, each row is taken to be a quantile.
#' @param delta the vector of noncentrality parameters.
#' @param sigma numeric matrix; scale matrix.
#' @param df degress of freedom.
#' @param log logical; whether to return log density value.
#' @examples
#' sigma <- crossprod(randn(3, 3))
#' data <- rmvt0(2, sigma = sigma, df = 2)
#' dmvt0(data, delta = rep(0, 3), sigma = sigma, df = 2)
#' @export
dmvt0 <- function(x, delta, sigma, df = 1, log = TRUE) {
  p <- nrow(sigma)
  a <- log(gamma(0.5 * (df + p))) - log(gamma(df / 2))
  inv_sigma <- solve(sigma)
  b <- 0.5 * p * log(df) + 0.5 * p * log(pi)
  b2 <- 0.5 * log(det(sigma))
  b3 <- 0.5 * (df + p)

  # cx <- t(t(x) - delta)
  # cx <- add_vector_to_matrix_row(-delta, x)
  cx <- x + one_matrix0(NROW(x), 1) %*% t(-delta)
  b4 <- log(1 + (rowSums(cx %*% inv_sigma * cx) / df))
  logretval <- (a - b - b2) - b3 * b4
  if (log) logretval else exp(logretval)
}


# Compute t(x) %*% A %*% x
quadratic_form <- function(x, A) {
  rowSums(x %*% A * x)
}


#' The density of the multivariate normal distribution
#' @param x vector or matrix of quantiles. If x is a matrix, each row is taken to be a quantile.
#' @param mean numeric vector; the mean vector.
#' @param sigma numeric matrix; the covariance matrix.
#' @param log logical; if TRUE, returns the log value.
#' @export
dmvnorm0 <- function(x, mean, sigma, log = FALSE) {
  d <- length(mean)
  cx <- x + one_matrix0(NROW(x), 1) %*% t(-mean)
  logretval <- - 0.5 * quadratic_form(cx, solve(sigma)) -
    (d / 2) * log(2 * pi) - 0.5 * log(det(sigma))
  if (log) logretval else exp(logretval)
}


#' The density of the normal distribution
#' @param x vector or matrix of quantiles. If x is a matrix, each row is taken to be a quantile.
#' @param mean numeric vector; the mean vector.
#' @param sd numeric vector; the sd vector.
#' @param log logical; if TRUE, returns the log value.
#' @export
dnorm0 <- function(x, mean, sd, log = FALSE) {
  cx <- x + one_matrix0(NROW(x), 1) %*% t(-mean)
  logretval <- as.vector(-0.5 * (cx / sd)^2 - log(sd) - 0.5 * log(2 * pi))
  if (log) logretval else exp(logretval)
}


#' The density of the student-t distribution
#' @param x vector or matrix of quantiles. If x is a matrix, each row is taken to be a quantile.
#' @param df degrees of freedom (> 0, maybe non-integer).
#' @param log logical; if TRUE, returns the log value.
#' @export
dt0 <- function(x, df, log = FALSE) {
  logretval <- log(gamma((df + 1) / 2)) - 0.5 * log(df * pi) -
    log(gamma(df / 2)) - (df + 1) / 2 * log(1 + x^2 / df)
  if (log) logretval else exp(logretval)
}


#' The density of the gamma distribution
#' @param x vector or matrix of quantiles. If x is a matrix, each row is taken to be a quantile.
#' @param shape positive number; the shape parameter.
#' @param rate positive number; the rate parameter.
#' @param scale positive number; the scale parameter.
#' @param log logical; if TRUE, returns the log value.
#' @export
dgamma0 <- function(x, shape, rate, scale = 1 / rate, log = FALSE) {
  if (!missing(rate) && !missing(scale)) {
    if (abs(rate * scale - 1) < 1e-15)
      warning("specify 'rate' or 'scale' but not both")
    else
      stop("specify 'rate' or 'scale' but not both")
  }
  k <- shape
  theta <- scale
  logretval <- (k-1) * log(x) - x / theta - log(gamma(k)) - k * log(theta)
  if (log) logretval else exp(logretval)
}


#' The density of the chi-squared distribution
#' @param x vector or matrix of quantiles. If x is a matrix, each row is taken to be a quantile.
#' @param df degrees of freedom (> 0, maybe non-integer).
#' @param log logical; if TRUE, returns the log value.
#' @export
dchisq0 <- function(x, df, log = FALSE) {
  half_k <- df / 2
  logretval <- (half_k - 1) * log(x) - x / 2 - half_k * log(2) - log(gamma(half_k))
  if (log) logretval else exp(logretval)
}


#' The density of the exponential distribution
#' @param x vector or matrix of quantiles. If x is a matrix, each row is taken to be a quantile.
#' @param rate positive number; the rate parameter.
#' @param log logical; if TRUE, returns the log value.
#' @export
dexp0 <- function(x, rate = 1, log = FALSE) {
  retval <- rate * exp(-rate * x)
  if (log) log(retval) else retval
}


#' The density of the Wishart distribution
#' @param X numeric matrix.
#' @param df degrees of freedom (> 0, maybe non-integer).
#' @param Scale numeric matrix; the scale matrix.
#' @param log logical; if TRUE, returns the log value.
#' @export
dWishart0 <- function(X, df, Scale, log = FALSE) {
  p <- nrow(X)
  n <- df
  # assertthat::assert_that(n > p - 1)
  V <- Scale

  logretval <- (n - p - 1) / 2 * log(det(X)) - tr(solve(V, X)) / 2 -
    (n * p / 2) * log(2) - n / 2 * log(det(V)) - lmvgamma(p, n / 2)

  if (log) logretval else exp(logretval)
}

lmvgamma <- function(p, a) {
  if (p == 1) {
    log(gamma(a))
  } else {
    0.5 * (p - 1) * log(pi) + lmvgamma(p - 1, a) + log(gamma(a + 0.5 * (1 - p)))
  }
}


#' The density of the inverse Wishart distribution
#' @param X numeric matrix.
#' @param df degrees of freedom (> 0, maybe non-integer).
#' @param Scale numeric matrix; the scale matrix.
#' @param log logical; if TRUE, returns the log value.
#' @export
dinvWishart0 <- function(X, df, Scale, log = FALSE) {
  p <- nrow(X)
  v <- df
  # assertthat::assert_that(v > p - 1)
  Phi <- Scale

  logretval <- v / 2 * log(det(Phi)) - (v + p + 1) / 2 * log(det(X)) -
    0.5 * tr(Phi %*% solve(X)) - (v * p / 2) * log(2) - lmvgamma(p, v / 2)

  if (log) logretval else exp(logretval)
}
