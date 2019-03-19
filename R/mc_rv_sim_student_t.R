#' @include class_dual_def.R
NULL


#' Simulation by inverse transform
#' @param cdf A function; the distribution function of a random variable.
inverse_transform <- function(cdf) {
  # Root-finding requires continuous cdf
  single_sim <- function(u, param) {
    interval <- c(0, 10)
    .f <- function(x) do.call(cdf, append(list(x), param)) - u
    res <- uniroot(.f, interval, tol = .Machine$double.eps^0.75,
                   maxiter = 10000, extendInt = "upX")
    res$root
  }
  function(n, ...) {
    param <- list(...)
    purrr::map_dbl(runif(n), ~single_sim(.x, param))
  }
}

rt_invt <- inverse_transform(pt)

#' Simulate random variates from the student-t distribution
#' @param n positive integer; number of observations.
#' @param df degrees of freedom (> 0, maybe non-integer). df = Inf is allowed.
#' @export
rt0 <- function(n, df) {
  assertthat::assert_that(length(df) == 1 || length(df) == n)
  if (length(df) == 1) {
    return(rt_invt(n = n, df = df))
  } else if (length(df) == n) {
    return(purrr::map_dbl(1:n, ~rt_invt(n = 1, df = df[.x])))
  }
}

#' Simulate random variates from the student-t distribution
#' @param n A scalar; the random sample.
#' @param df A dual number; the degree of freedom.
setMethod("rt0",
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
       res <- mapreduce(1:n, ~simulate_single(df), rbind2)
     } else if (length(df) == n) {
       res <- mapreduce(1:n, ~simulate_single(df[.x]), rbind2)
     }
     res
   }
)

d_student_t <- function(x, df) {
  # x is a scalar, dgf0 is a dual number
  h <- 1e-8
  numerator <- - (pt(x, df = df + h) - pt(x, df = df)) / h
  denominator <- dt(x, df = df)
  numerator / denominator
}
