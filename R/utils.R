#' Generate a matrix randomly from a normal distribution
#' @param ... Any number of dimensions.
#' @param mean vector of means.
#' @param sd vector of standard deviations.
#' @export
randn <- function(..., mean = 0, sd = 1) {
  structure(rnorm(prod(...), mean = mean, sd = sd), dim = c(...))
}


#' Generate a matrix randomly from a uniform distribution
#' @param ... Any number of dimensions.
#' @param min lower limit of the distribution. Must be finite.
#' @param max upper limit of the distribution. Must be finite.
#' @export
randu <- function(..., min = 0, max = 1) {
  structure(runif(prod(...), min = min, max = max), dim = c(...))
}


mapreduce <- function(x, f, g, ...) {
  x %>% purrr::map(f, ...) %>% purrr::reduce(g)
}


map2reduce <- function(x, y, f, g, ...) {
  purrr::reduce(purrr::map2(x, y, f, ...), g)
}


map_row <- function(m0, f) {
  purrr::map(1:nrow(m0), ~f(m0[.x, ]))
}


call_S4 <- function(f_name, ...) {
  getMethod(f_name, list(...))@.Data
}


# map_col <- function(m0) {
#   purrr::map(1:ncol(m0), ~f(m0[, .x]))
# }
