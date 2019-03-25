#' Finite differencing
#' @param f A function of which the derivative is seeked.
#' @param vary A named list of variables; the variables to be varied.
#' @param fix A named list of variables; the variables to be fixed.
#' @param h The finite differencing parameter; the size of perturbation.
#' @param seed Seed; for pathwise derivative only.
#' @examples
#' f <- function(y, X, beta) { y - X %*% beta }
#' finite_diff(f,
#'   vary = list(beta = c(5, 6)),
#'   fix = list(X = matrix(1:4, 2, 2), y = c(2, 3))
#' )
#'
#' g <- function(X, Y) { X %*% Y }
#' finite_diff(g, vary = list(X = randn(2, 2), Y = randn(2, 2)))
#'
#' @export
finite_diff <- function(f, vary, fix = NULL, h = 1e-8, seed) {
  has_seed <- !missing(seed)
  x <- unlist(vary)
  f_vec <- function(vec0) {
    if (has_seed) set.seed(seed)
    do.call(f, append(relist(vec0, vary), fix))
  }
  vec_finite_diff(f_vec, x, h) %>%
    name_matrix(vary)
}

vec_finite_diff <- function(f, x, h = 1e-8) {
  ufx <- as.numeric(unlist(f(x)))  # avoid duplicate evaluation
  finite_deriv <- function(x) {
    (as.numeric(unlist(f(x))) - ufx) / h
  }
  perturbate <- function(v, h) {
    purrr::map(seq_along(v), function(i) { v[i] <- v[i] + h; v })
  }
  perturbate(x, h) %>%
    purrr::map(finite_deriv) %>%
    do.call(cbind, .)
}

name_matrix <- function(x, input) {
  colnames(x) <- paste("d", names(unlist(input)), sep = "_")
  if (is.null(rownames(x))) {
    rownames(x) <- paste("d_output", seq(nrow(x)), sep = "_")
  }
  x
}
