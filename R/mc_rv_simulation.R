#' @include class_dual_def.R
NULL

setClassUnion("dual_or_numeric", c("dual", "numeric"))

# A note about the implementation
# The 'parameter-recycling' mechanism is not implemented for random variable
# simulation involving 'dual' number since it penalises performance and it
# justifies the practice of supplying the wrong number of parameters.
# The mechanism is implemented in `rgamma1` (which does not involve 'dual'
# number) to match `rgamma` from the base 'stats' package.

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


#===============================================================
# Gamma distribution
#===============================================================
#' Simulate gamma random variates
#' @param n Positive integer; the number of samples.
#' @param shape Positive number; the shape of the gamma distribution.
#' @param scale Positive number; the scale of the gamma distribution.
#' @param method base or inv_tf; base refers to `stats::rgamma` while
#' inv_tf refers to inverse transform.
#' @note Inverse transform is slower, but it is provided to remedy that
#' base R uses two algorithms to simulate gamma random variables
#' in different parameter regions, which creates a discontinuity in the
#' pathwise derivative.
#' @export
rgamma0 <- function(n, shape, scale, method = "base") {
  if (method == "base")   return(rgamma(n, shape, scale = scale))
  if (method == "inv_tf") return(rgamma1(n, shape, scale = scale))
  stop("method must be one of 'base' and 'inv_tf'.")
}


# Inverse-transform by root-finding
rgamma1 <- function(n, shape, scale = 1) {
  len_shape <- length(shape)
  len_scale <- length(scale)
  if ((len_shape == 1) && (len_scale == 1)) {
    return(purrr::map_dbl(runif(n), gamma_inv_tf, shape = shape, scale = scale))
  } else {
    shape <- rep(shape, ceiling(n / len_shape))
    scale <- rep(scale, ceiling(n / len_scale))
    return(purrr::pmap_dbl(
      list(u = runif(n), shape = shape[1:n], scale = scale[1:n]),
      gamma_inv_tf
    ))
  }
}
gamma_inv_tf <- function(u, shape, scale) {
  interval <- c(0, shape + 2 * sqrt(shape / scale^2))
  .f <- function(x) { pgamma(x, shape = shape, scale = scale) - u }
  uniroot(
    .f, interval, tol = .Machine$double.eps^0.75, maxiter = 10000,
    extendInt = "upX"
  )$root
}


#' Simulate gamma random variates
#' @param n Positive integer; the number of samples.
#' @param shape A dual number or a scalar; the shape of the gamma distribution.
#' @param scale A dual number or a scalar; the scale of the gamma distribution.
#' @param method base or inv_tf; base refers to `stats::rgamma` while
#' inv_tf refers to inverse transform.
#' @note At least on of `shape` and `scale` should be a dual number.
#' Otherwise, `stats::rgamma` is called instead.
#' @rdname gamma_rv
setMethod("rgamma0",
  signature(n = "numeric", shape = "dual", scale = "dual_or_numeric"),
  function(n, shape, scale, method = "base") {
    len_shape <- length(shape)
    assertthat::assert_that(len_shape == 1 || len_shape == n)

    len_scale <- length(scale)
    assertthat::assert_that(len_scale == 1 || len_scale == n)

    rgamma_single <- function(shape, scale) {
      g <- rgamma0(1, parent_of(shape), scale = 1, method = method)
      d_g <- d_rgamma_num_dual(g, shape)
      g_dual <- new("dual", x = g, dx = deriv_of(d_g), param = shape@param)
      g_dual * scale
    }
    # handle parameters of different length
    if ((len_shape == 1) && (len_scale == 1)) {
      return(mapreduce(seq(n), ~rgamma_single(shape, scale), rbind2))
    }
    if ((len_shape == 1) && (len_scale == n)) {
      return(mapreduce(seq(n), ~rgamma_single(shape, scale[.x]), rbind2))
    }
    if ((len_shape == n) && (len_scale == 1)) {
      return(mapreduce(seq(n), ~rgamma_single(shape[.x], scale), rbind2))
    }
    if ((len_shape == n) && (len_scale == n)) {
      return(mapreduce(seq(n), ~rgamma_single(shape[.x], scale[.x]), rbind2))
    }
  }
)

#' Simulate gamma random variates
#' @rdname gamma_rv
setMethod("rgamma0",
  signature(n = "numeric", shape = "numeric", scale = "dual"),
  function(n, shape, scale, method = "base") {
    len_shape <- length(shape)
    assertthat::assert_that(len_shape == 1 || len_shape == n)

    len_scale <- length(scale)
    assertthat::assert_that(len_scale == 1 || len_scale == n)

    rgamma_single <- function(shape, scale) {
      g <- rgamma0(1, shape, scale = 1, method = method)
      g * scale
    }

    # handle parameters of different length
    if ((len_shape == 1) && (len_scale == 1)) {
      return(mapreduce(seq(n), ~rgamma_single(shape, scale), rbind2))
    }
    if ((len_shape == 1) && (len_scale == n)) {
      return(mapreduce(seq(n), ~rgamma_single(shape, scale[.x]), rbind2))
    }
    if ((len_shape == n) && (len_scale == 1)) {
      return(mapreduce(seq(n), ~rgamma_single(shape[.x], scale), rbind2))
    }
    if ((len_shape == n) && (len_scale == n)) {
      return(mapreduce(seq(n), ~rgamma_single(shape[.x], scale[.x]), rbind2))
    }
  }
)

# d_rgamma <- function(g, alpha) {
#   # takes an simulated value from gamma distribution and the corresponding
#   # parameter alpha, returns the derivative "d G(alpha, 1) / d alpha"
#   f <- function(t) { log(t) * dgamma(t, alpha, 1) }
#   num_1 <- integrate(f, 0, g)$value
#   num_2 <- digamma(alpha) * pgamma(g, alpha, 1)
#   - (num_1 - num_2) / dgamma(g, alpha, 1)
# }

d_rgamma_num_dual <- function(g, alpha) {
  # g is a scalar, alpha is a dual number.
  alpha0 <- parent_of(alpha)
  f <- function(t) { log(t) * dgamma(t, alpha0, 1) }
  num_1 <- integrate(f, 0, g)$value
  num_2 <- digamma(alpha0) * pgamma(g, alpha0, 1)
  - (num_1 - num_2) / dgamma(g, alpha0, 1) * alpha  # correct for derivative
}

# Info: Inverse-Gamma
# X ~ Gamma(alpha, beta) <=> X^{-1} ~ IGamma(alpha, beta)
# (alpha, beta) ~ (shape, rate) parametrisation
# https://en.wikipedia.org/wiki/Inverse-gamma_distribution


#===============================================================
# Wishart / Inverse-Wishart distribution
#===============================================================
#' Simulate Wishart random variates
#' @param v A scalar; degrees of freedom.
#' @param M A matrix; the matrix parameter of the Wishart distribution.
#' @param method base or inv_tf; base refers to the function in the
#' `stats` package while inv_tf refers to inverse transform.
#' @export
rWishart0 <- function(v, M, method = "base") {
  if (method == "base") return(stats::rWishart(1, v, M)[,,1])

  chi_samples <- seq(nrow(M)) %>%
    purrr::map_dbl(~rchisq0(1, v - .x + 1, method = method))
  A <- diag(sqrt(chi_samples))
  A[lower.tri(A)] <- rnorm(sum(lower.tri(A)))
  L <- chol0(M)
  tcrossprod(L %*% A)
}

#' Simulating from Wishart distribution using Bartlett decomposition
#' @param v A dual number or a scalar; degrees of freedom.
#' @param M A dual number or a matrix; the matrix parameter of the Wishart distribution.
#' @param method base or inv_tf; base refers to the function in the
#' `stats` package while inv_tf refers to inverse transform.
#' @rdname wishart_rv
setMethod("rWishart0",
  signature(v = "dual", M = "dual"),
  function(v, M, method = "base") {
    n <- nrow(parent_of(M))
    # Implement: diag(A) <- sqrt(rchisq(n, df = v - seq(n) + 1))
    chisq_samples <- seq(n) %>%
      purrr::map(~rchisq0(1, v - .x + 1, method = method)) %>%
      purrr::reduce(rbind2)
    A <- diag(as.vector(sqrt(chisq_samples)))
    # Changing lower triangular part of the matrix; this doesn't change dA
    A_x <- parent_of(A)
    A_x[lower.tri(A_x)] <- rnorm(sum(lower.tri(A_x)))
    A@x <- A_x

    L <- chol0(M)
    tcrossprod(L %*% A)
  }
)

#' @rdname wishart_rv
setMethod("rWishart0",
  signature(v = "ANY", M = "dual"),
  function(v, M, method = "base") {
    n <- nrow(parent_of(M))
    # Implement: diag(A) <- sqrt(rchisq(n, df = v - seq(n) + 1))
    chisq_samples <- seq(n) %>%
      purrr::map_dbl(~rchisq0(1, v - .x + 1, method = method))
    A <- diag(sqrt(chisq_samples))
    A[lower.tri(A)] <- rnorm(sum(lower.tri(A)))
    L <- chol0(M)
    tcrossprod(L %*% A)
  }
)

#' @rdname wishart_rv
setMethod("rWishart0",
  signature(v = "dual", M = "ANY"),
  function(v, M, method = "base") {
    n <- nrow(M)
    # Implement: diag(A) <- sqrt(rchisq(n, df = v - seq(n) + 1))
    chisq_samples <- seq(n) %>%
      purrr::map(~rchisq0(1, v - .x + 1, method = method)) %>%
      purrr::reduce(rbind2)
    A <- diag(as.vector(sqrt(chisq_samples)))
    # Changing lower triangular part of the matrix; this doesn't change dA
    A_x <- parent_of(A)
    A_x[lower.tri(A_x)] <- rnorm(sum(lower.tri(A_x)))
    A@x <- A_x

    L <- chol0(M)
    tcrossprod(L %*% A)
  }
)


# Info: Inverse-Wishart
# X ~ W^{-1}(Sigma, v) <=> X^{-1} ~ W(Sigma^{-1}, v)
# i.e. to draw X from Inverse-Wishart, one draws Y from Wishart
# with inverse Sigma, then set X = Y^{-1}.
# https://en.wikipedia.org/wiki/Inverse-Wishart_distribution


#===============================================================
# Exponential and Chi-squared distribution
#===============================================================
#' Simulate Chi-square random variates
#' @param n Positive integer; number of observations.
#' @param df Degrees of freedom (can be non-integer).
#' @param method base or inv_tf; base refers to `stats::rchisq` while
#' inv_tf refers to inverse transform.
#' @export
rchisq0 <- function(n, df, method = "base") {
  if (method == "base") return(stats::rchisq(n, df))
  rgamma0(n, df / 2, scale = 2, method = method)
}

#' Simulate Chi-square random variates
#' @param n Positive integer; number of observations.
#' @param df Degrees of freedom (can be non-integer).
#' @param method base or inv_tf; base refers to the function in the
#' `stats` package while inv_tf refers to inverse transform.
setMethod("rchisq0",
  signature(n = "numeric", df = "dual"),
  function(n, df, method = "base") {
    rgamma0(n, df / 2, scale = 2, method = method)
  }
)


#' Simulate exponential random variates
#' @inherit stats::rexp
#' @export
rexp0 <- function(n, rate = 1) {
  stats::rexp(n, rate)
}


#' Simulate exponential random variates
#' @param n Positive integer; the number of samples.
#' @param rate A dual number; the rate of the exponential distribution.
setMethod("rexp0",
  signature(n = "numeric", rate = "dual"),
  function(n, rate) {
    len_rate <- length(rate)
    assertthat::assert_that(len_rate == 1 || len_rate == n)
    as.matrix(stats::rexp(n)) / rate
  }
)
