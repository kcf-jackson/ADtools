#' @include class_dual_def.R
NULL


#' Simulation by inverse transform
#' @param cdf A function; the distribution function of a random variable.
inverse_transform <- function(cdf) {
  # Root-finding requires continuous cdf
  function(n, ...) {                                         # nocov
    param <- list(...)
    purrr::map_dbl(runif(n), ~ single_sim(.x, param, cdf))
  }                                                          # nocov
}

single_sim <- function(u, param, cdf) {
  interval <- c(0, 10)
  .f <- function(x) do.call(cdf, append(list(x), param)) - u
  res <- uniroot(.f, interval,
    tol = .Machine$double.eps^0.75,
    maxiter = 10000, extendInt = "upX"
  )
  res$root
}

rt_invt <- inverse_transform(pt)


#' Simulate random variates from the student-t distribution
#'
#' @param n positive integer; the number of samples.
#' @param df degrees of freedom (> 0, maybe non-integer). df = Inf is allowed.
#'
#' @examples
#' n <- 10
#' rt0(10, df = 3)
#'
#' @export
rt0 <- function(n, df) {
  assertthat::assert_that(length(df) == 1 || length(df) == n)
  if (length(df) == 1) {
    return(rt_invt(n = n, df = df))
  } else if (length(df) == n) {
    return(purrr::map_dbl(df, ~ rt_invt(n = 1, df = .x)))
  }
}

#' Simulate random variates from the student-t distribution
#' @param n positive integer; the number of samples.
#' @param df A dual number; the degree of freedom.
setMethod(
  "rt0",
  signature(n = "numeric", df = "dual"),
  function(n, df) {
    assertthat::assert_that(length(df) == 1 || length(df) == n)

    simulate_single <- function(df) {
      px <- rt0(1, df = df@x)
      dx <- d_student_t(px, df@x) * df@dx
      df@x <- px
      df@dx <- dx
      df
    }

    if (length(df) == 1) {
      res <- mapreduce(1:n, ~ simulate_single(df), rbind2)
    } else if (length(df) == n) {
      res <- mapreduce(1:n, ~ simulate_single(df[.x]), rbind2)
    }
    res
  }
)

d_student_t <- function(x, df) {
  h <- 1e-8
  numerator <- -(pt(x, df = df + h) - pt(x, df = df)) / h
  denominator <- dt(x, df = df)
  numerator / denominator
}


#' Simulate random variates from the multivariate t distribution
#'
#' @param n positive integer; number of observations.
#' @param sigma scale matrix.
#' @param df degree of freedom.
#' @param delta non-centrality parameters.
#'
#' @examples
#' n <- 10
#' d <- 3
#' rmvt0(n, sigma = crossprod(randn(d, d)), df = 2)
#'
#' @export
rmvt0 <- function(n, sigma, df, delta = 0) {
  sim <- function(v, mu, Sigma) {
    sqrt_W <- sqrt(v / rchisq0(1, df = v))
    A <- chol0(Sigma)
    p <- nrow(Sigma)
    Z <- rmvnorm0(1, mean = rep(0, p), sigma = diag(p))
    t(mu + sqrt_W * A %*% t(Z))
  }

  mapreduce(1:n, ~ sim(v = df, mu = delta, Sigma = sigma), rbind2)
}
