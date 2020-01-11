#' The density of the multivariate t distribution
#' @param x vector or matrix of quantiles. If x is a matrix, each row is taken to be a quantile.
#' @param delta the vector of noncentrality parameters.
#' @param sigma scale matrix.
#' @param df degress of freedom.
#' @param log logical; whether to return log density value.
#' @export
#' @examples
#' sigma <- crossprod(randn(3, 3))
#' data <- rmvt0(2, sigma = sigma, df = 2)
#' dmvt0(data, delta = rep(0, 3), sigma = sigma, df = 2)
dmvt0 <- function(x, delta, sigma, df = 1, log = TRUE) {
  p <- nrow(sigma)
  a <- log(gamma(0.5 * (df + p))) - log(gamma(df / 2))
  inv_sigma <- solve(sigma)
  b <- 0.5 * p * log(df) + 0.5 * p * log(pi)
  b2 <- 0.5 * log(det(sigma))
  b3 <- 0.5 * (df + p)

  # datum_LL <- function(x) {
  #   b4 <- log(1 + as.vector(t(x - delta) %*% inv_sigma %*% (x - delta)) / df)
  #   logretval <- a - b - b2 - b3 * b4
  #   if (log) logretval else exp(logretval)
  # }
  # mapreduce(1:nrow(x), ~datum_LL(x[.x, ]), rbind2)

  # cx <- t(t(x) - delta)
  delta_mat <- t(mapreduce(1:nrow(x), ~delta, cbind2))
  cx <- x - delta_mat
  b4 <- log(1 + (rowSums(cx %*% inv_sigma * cx) / df))
  logretval <- a - b - b2 - b3 * b4
  if (log) logretval else exp(logretval)
}

# dmvnorm0
# dnorm0
# dt0
# dmvt0
# dgamma0
# dWishart0
# dchisq0
# dexp0
#
