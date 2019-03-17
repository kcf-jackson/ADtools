#' @include class_dual_def.R
NULL

#' Simulate random variates from the student-t distribution
#' @param n positive integer; number of observations.
#' @param df degrees of freedom (> 0, maybe non-integer). df = Inf is allowed.
#' @export
rt0 <- function(n, df) {
  stats::rt(n, df)
}

#' Simulate random variates from the student-t distribution
#' @param n A scalar; the random sample.
#' @param df A dual number; the degree of freedom.
setMethod("rt0",
   signature(n = "numeric", df = "dual"),
   function(n, df) {
     assertthat::assert_that(length(df) == 1 || length(df) == n)

     simulate_single <- function(df) {
       px <- df@x  # the degree of freedom
       s <- rt0(1, px)
       new("dual", x = s,
           dx = d_student_t(s, px) * deriv_of(df),
           param = param_of(df))
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
