element_wise_derivative <- function(f, df) {
  function(x) {
    px <- parent_of(x)
    dx <- deriv_of(x)

    x@x <- f(px)
    for (i in 1:nrow(dx)) {
      dx[i, ] <- df(px[i]) * dx[i, ]
    }
    x@dx <- dx

    return(x)
  }
}

d_sqrt <- function(x) { 0.5 / sqrt(x) }

#' Element-wise square-root of a dual object
#' @param x A "dual" object.
setMethod("sqrt",
          signature(x = "dual"),
          element_wise_derivative(sqrt, d_sqrt)
)
