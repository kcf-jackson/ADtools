#' @include class_dual_def.R
NULL

# A note about the implementation
# The 'parameter-recycling' mechanism is not implemented for random variable
# simulation involving 'dual' number since it penalises performance and it
# justifies the practice of supplying the wrong number of parameters.
# The mechanism is implemented in `rgamma1` (which does not involve 'dual'
# number) to match `rgamma` from the base 'stats' package.

#===============================================================
# Gamma distribution
#---------------------------------------------------------------
#' Simulate gamma random variates
#'
#' @param n Positive integer; the number of samples.
#' @param shape Positive number; the shape of the gamma distribution.
#' @param scale Positive number; the scale of the gamma distribution.
#' @param method 'base' or 'inv_tf'; 'base' refers to `stats::rgamma` while
#' 'inv_tf' refers to inverse transform.
#'
#' @note Inverse transform is slower, but it is provided to remedy that
#' base R uses two algorithms to simulate gamma random variables
#' in different parameter regions, which creates a discontinuity in the
#' pathwise derivative.
#'
#' @examples
#' n <- 10
#' rgamma0(n, shape = 1, scale = 1)
#'
#' @export
rgamma0 <- function(n, shape, scale, method = "inv_tf") {
  if (method == "base")   return(rgamma(n, shape, scale = scale))
  if (method == "inv_tf") return(rgamma1(n, shape, scale = scale))
  stop("method must be one of 'base' and 'inv_tf'.")
}

# Inverse-transform by root-finding
rgamma_invt <- inverse_transform(pgamma)
rgamma1 <- function(n, shape, scale = 1) {
  len_shape <- length(shape)
  len_scale <- length(scale)
  assertthat::assert_that(len_shape == 1 || len_shape == n)
  assertthat::assert_that(len_scale == 1 || len_scale == n)

  param <- cbind(1:n, shape, scale)
  purrr::map2_dbl(param[, 2], param[, 3], ~rgamma_invt(1, shape = .x, scale = .y))
}


rgamma0_dual <- function(n, shape, scale, method = "inv_tf") {
  len_shape <- length(shape)
  assertthat::assert_that(len_shape == 1 || len_shape == n)

  len_scale <- length(scale)
  assertthat::assert_that(len_scale == 1 || len_scale == n)

  rgamma_single <- function(shape) {
    g_dual <- shape
    g <- rgamma0(1, shape@x, scale = 1, method = method)
    g_dual@x <- g
    g_dual@dx <- d_rgamma(g, shape@x) * shape@dx
    g_dual
  }
  if (len_shape == 1) {
    return(mapreduce(seq(n), ~rgamma_single(shape), rbind2) * scale)
  } else {
    return(mapreduce(seq(n), ~rgamma_single(shape[.x]), rbind2) * scale)
  }
  # rgamma_single <- function(shape, scale) {
  #   g_dual <- shape
  #   g <- rgamma0(1, shape@x, scale = 1, method = method)
  #   g_dual@x <- g
  #   g_dual@dx <- d_rgamma(g, shape@x) * shape@dx
  #   g_dual * scale
  # }
  # handle parameters of different length
  # if ((len_shape == 1) && (len_scale == 1)) {
  #   return(mapreduce(seq(n), ~rgamma_single(shape, scale), rbind2))
  # }
  # if ((len_shape == 1) && (len_scale == n)) {
  #   return(mapreduce(seq(n), ~rgamma_single(shape, scale[.x]), rbind2))
  # }
  # if ((len_shape == n) && (len_scale == 1)) {
  #   return(mapreduce(seq(n), ~rgamma_single(shape[.x], scale), rbind2))
  # }
  # if ((len_shape == n) && (len_scale == n)) {
  #   return(mapreduce(seq(n), ~rgamma_single(shape[.x], scale[.x]), rbind2))
  # }
}


#' Simulate gamma random variates
#'
#' @param n Positive integer; the number of samples.
#' @param shape A dual number or a scalar; the shape of the gamma distribution.
#' @param scale A dual number or a scalar; the scale of the gamma distribution.
#' @param method 'base' or 'inv_tf'; 'base' refers to `stats::rgamma` while
#' 'inv_tf' refers to inverse transform.
#'
#' @note At least one of `shape` and `scale` should be a dual number.
#'
#' @rdname gamma-rv
setMethod("rgamma0",
  signature(n = "numeric", shape = "dual", scale = "dual"),
  rgamma0_dual
)

#' @rdname gamma-rv
setMethod("rgamma0",
  signature(n = "numeric", shape = "dual", scale = "numeric"),
  rgamma0_dual
)

d_rgamma <- function(g, alpha0) {
  h <- 1e-8
  numerator <- -(pgamma(g, shape = alpha0 + h) - pgamma(g, shape = alpha0)) / h
  denominator <- dgamma(g, shape = alpha0)
  numerator / denominator
}

#' @rdname gamma-rv
setMethod("rgamma0",
  signature(n = "numeric", shape = "numeric", scale = "dual"),
  function(n, shape, scale, method = "inv_tf") {
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


# Info: Inverse-Gamma
# X ~ Gamma(alpha, beta) <=> X^{-1} ~ IGamma(alpha, beta)
# (alpha, beta) ~ (shape, rate) parametrisation
# https://en.wikipedia.org/wiki/Inverse-gamma_distribution


#===============================================================
# Wishart / Inverse-Wishart distribution
#---------------------------------------------------------------
#' Simulate Wishart random variates
#'
#' @param v A scalar; degrees of freedom.
#' @param M A matrix; the matrix parameter of the Wishart distribution.
#' @param method base or inv_tf; base refers to the function in the
#' `stats` package while inv_tf refers to inverse transform.
#'
#' @examples
#' d <- 4
#' rWishart0(3, crossprod(randn(d, d)))
#'
#' @export
rWishart0 <- function(v, M, method = "inv_tf") {
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
  function(v, M, method = "inv_tf") {
    n <- nrow(M@x)
    # Implement: diag(A) <- sqrt(rchisq(n, df = v - seq(n) + 1))
    chisq_samples <- mapreduce(
      seq(n), ~rchisq0(1, v - .x + 1, method = method), rbind2
    )
    A <- diag(as.vector(sqrt(chisq_samples)))
    # Changing lower triangular part of the matrix; this doesn't change dA
    A_x <- A@x
    A_x[lower.tri(A_x)] <- rnorm(sum(lower.tri(A_x)))
    A@x <- A_x

    L <- chol0(M)
    tcrossprod(L %*% A)
  }
)

#' @rdname wishart_rv
setMethod("rWishart0",
  signature(v = "ANY", M = "dual"),
  function(v, M, method = "inv_tf") {
    n <- nrow(M@x)
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
  function(v, M, method = "inv_tf") {
    n <- nrow(M)
    # Implement: diag(A) <- sqrt(rchisq(n, df = v - seq(n) + 1))
    chisq_samples <- mapreduce(
      seq(n), ~rchisq0(1, v - .x + 1, method = method), rbind2
    )
    A <- diag(as.vector(sqrt(chisq_samples)))
    # Changing lower triangular part of the matrix; this doesn't change dA
    A_x <- A@x
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
#---------------------------------------------------------------
#' Simulate Chi-square random variates
#'
#' @param n Positive integer; number of observations.
#' @param df Degrees of freedom (can be non-integer).
#' @param method base or inv_tf; base refers to `stats::rchisq` while
#' inv_tf refers to inverse transform.
#'
#' @examples
#' n <- 10
#' rchisq0(n, df = 3)
#'
#' @export
rchisq0 <- function(n, df, method = "inv_tf") {
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
  function(n, df, method = "inv_tf") {
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
