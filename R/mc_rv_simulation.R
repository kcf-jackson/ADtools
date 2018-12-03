#' #' @include class_dual_def.R
#' NULL
#'
#' setMethod("rnorm",
#'           signature(n = "integer", mean = "dual", sd = "dual"),
#'           function(n, mean, sd) { mean + sd * rnorm(n) }
#' )
#'
#' setGeneric("rmvnorm", function(n, mean, sigma) {
#'   standardGeneric("rmvnorm")
#' })
#'
#' setMethod("rmvnorm",
#'           signature(n = "integer", mean = "dual", sigma = "dual"),
#'           function(n, mean, sigma) {
#'             mean + chol0(sigma) %*% mvtnorm::rmvnorm(n, numeric(length(mean)))
#'           }
#' )
#'
#' # setMethod("rgamma0",
#' #   signature(n = "integer", shape = "dual", scale = "dual"),
#' #   function(n, shape, scale) {
#' #     g <- rgamma(1, shape@x, scale = 1)
#' #     d_g <- d_rgamma(g, shape@x)
#' #
#' #     new("dual", x = g * scale@x, dx = d_g * scale@x, param = shape@param)
#' #
#' #     # {shape, scale} are {k, theta}, {shape, rate} are {alpha, beta}
#' #     g <- rgamma(1, shape)
#' #     scale * rgamma(n, shape)
#' #   }
#' # )
#'
#' # setMethod("rgamma0",
#' #           signature(n = "integer", shape = "dual", scale = "ANY"),
#' #           function() {}
#' # )
#' #
#' # setMethod("rgamma0",
#' #           signature(n = "integer", shape = "ANY", scale = "dual"),
#' #           function() {}
#' # )
#'
#' d_rgamma <- function(g, alpha) {
#'   # takes an simulated value from gamma distribution and the corresponding
#'   # parameter alpha, returns the derivative "d G(alpha, 1) / d alpha"
#'   f <- function(t) { log(t) * dgamma(t, alpha, 1) }
#'   num_1 <- integrate(f, 0, g)$value
#'   num_2 <- digamma(alpha) * pgamma(g, alpha, 1)
#'   - (num_1 - num_2) / dgamma(g, alpha, 1)
#' }
#'
#'
#' setMethod("rchisq",
#'           signature(n = "integer", df = "dual"),
#'           function(n, df) { rgamma0(n, df / 2, 2) }
#' )
#'
#' setMethod("rexp",
#'           signature(n = "integer", rate = "dual"),
#'           function(n, rate) {
#'             # Note that exp(rate = Lambda)
#'             # = gamma(shape = 1, scale = 1 / Lambda)
#'             # = gamma(shape = 1, rate = Lambda)
#'             rgamma(n, shape = 1, scale = 1 / rate)
#'           }
#' )
