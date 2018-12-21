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


#' Single-level destructuring
#' @description Strict usage: c(x, y) %<-% list(a, b, ...).
#' LHS must start with 'c'; RHS must be a list;
#' RHS must have at least as many elements as LHS.
#' @note Multiple-level is avoided intentionally, i.e. nested destructuring
#' is not allowed.
#' @examples
#' c(x, y) %<-% list(1, 2)
#' c(x, y) %<-% list(1, 2, 3)
`%<-%` <- function(lhs, rhs) {
  lhs_chr <- Map(deparse, substitute(lhs))  # Parse LHS into characters
  variables <- lhs_chr[-1]

  if (lhs_chr[[1]] != "c") stop("LHS must start with 'c'.")
  if (!is.list(rhs)) stop("RHS must be a list.")
  if (length(rhs) < length(variables))
    stop("RHS must have at least as many elements as LHS.")

  for (i in seq_along(variables)) {
    assign(variables[[i]], rhs[[i]], envir = parent.frame())
  }
  invisible(NULL)
}
