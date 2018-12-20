#' Generate a matrix randomly
#' @param ... Any number of dimensions.
#' @export
randn <- function(...) {
  structure(rnorm(prod(...)), dim = c(...))
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


# map_col <- function(m0) {
#   purrr::map(1:ncol(m0), ~f(m0[, .x]))
# }
