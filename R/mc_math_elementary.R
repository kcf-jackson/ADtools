element_wise_derivative <- function(f, df) {
  function(x) {                          # nocov
    px <- x@x
    dx <- x@dx

    x@x <- f(px)
    for (i in 1:nrow(dx)) {
      dx[i, ] <- df(px[i]) * dx[i, ]
    }                                    # nocov
    x@dx <- dx

    return(x)
  }                                      # nocov
}

lambda <- pryr::f

#' Element-wise square-root of a dual object
#' @param x A "dual" object.
setMethod(
  "sqrt",
  signature(x = "dual"),
  element_wise_derivative(sqrt, lambda(0.5 / sqrt(x)))
)

#' Element-wise sine of a dual object
#' @param x A "dual" object.
setMethod(
  "sin",
  signature(x = "dual"),
  element_wise_derivative(sin, cos)
)

#' Element-wise cosine of a dual object
#' @param x A "dual" object.
setMethod(
  "cos",
  signature(x = "dual"),
  element_wise_derivative(cos, lambda(-sin(x)))
)

#' Element-wise tangent of a dual object
#' @param x A "dual" object.
setMethod(
  "tan",
  signature(x = "dual"),
  element_wise_derivative(tan, lambda(cos(x)^{-2}))
)

#' Element-wise exponential of a dual object
#' @param x A "dual" object.
setMethod(
  "exp",
  signature(x = "dual"),
  element_wise_derivative(exp, exp)
)

#' Element-wise logarithm of a dual object
#' @param x A "dual" object.
setMethod(
  "log",
  signature(x = "dual"),
  element_wise_derivative(log, lambda(1 / x))
)

#' Element-wise gamma of a dual object
#' @param x A "dual" object.
setMethod(
  "gamma",
  signature(x = "dual"),
  element_wise_derivative(gamma, lambda(digamma(x) * gamma(x)))
)
