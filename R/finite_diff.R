#' Finite difference method
#'
#' @param f A function of which the derivative is sought.
#' @param wrt Character vector; the name of the variables to differentiate with respect to.
#' @param at A named list of variables; the point at which the derivative is evaluated.
#' @param h The finite differencing parameter; the size of perturbation.
#' @param seed Seed; for pathwise derivative only.
#' @param method "forward" for forward differencing or "central" for central differencing.
#'
#' @return A numeric matrix; the Jacobian matrix.
#'
#' @examples
#' f <- function(y, X, beta) { y - X %*% beta }
#' finite_diff(
#'   f, wrt = "beta",
#'   at = list(X = matrix(1:4, 2, 2), y = c(2, 3), beta = c(5, 6))
#' )
#'
#' g <- function(X, Y) { X %*% Y }
#' finite_diff(g, at = list(X = randn(2, 2), Y = randn(2, 2)))
#'
#' h <- function(x) { exp(-x^2) }
#' finite_diff(h, at = list(x = 0.01), wrt = "x")
#' finite_diff(h, at = list(x = 0.01), wrt = "x", method = "central")
#'
#' @export
finite_diff <- function(f, wrt = NULL, at, h = 1e-8, seed,
                        method = "forward") {
  f_args <- formalArgs(f)
  wrt <- wrt %||% f_args
  vary <- at[wrt]
  fix <- at[setdiff(f_args, wrt)]

  has_seed <- !missing(seed)
  x <- unlist(vary)
  f_vec <- function(vec0) {
    if (has_seed) set.seed(seed)
    do.call(f, append(relist(vec0, vary), fix))
  }

  if (method == "forward") {
    vec_finite_diff(f_vec, x, h) %>%
      name_matrix(vary)
  } else if (method == "central") {
    vec_central_diff(f_vec, x, h) %>%
      name_matrix(vary)
  } else {
    stop("The argument 'method' must be 'forward' or 'central'.")
  }
}

vec_finite_diff <- function(f, x, h = 1e-8) {
  to_vec <- function(x) as.numeric(unlist(x))

  ufx <- to_vec(f(x))  # avoid duplicate evaluation
  finite_deriv <- function(x) {
    (to_vec(f(x)) - ufx) / h
  }

  perturb <- function(v, h) {
    purrr::map(seq_along(v), function(i) { v[i] <- v[i] + h; v })
  }

  perturb(x, h) %>%
    purrr::map(finite_deriv) %>%
    do.call(cbind, .)
}

vec_central_diff <- function(f, x, h = 1e-8) {
  to_vec <- function(x) as.numeric(unlist(x))

  finite_deriv <- function(x, i) {
    y <- x
    x[i] <- x[i] + h
    y[i] <- y[i] - h
    (to_vec(f(x)) - to_vec(f(y))) / (2 * h)
  }

  seq_along(x) %>%
    purrr::map(~finite_deriv(x, i = .x)) %>%
    do.call(cbind, .)
}

# Alternative implementation (where perturbation is done on-the-fly)
# Note: It turns out this implementation does not necessarily reduce
# memory footprint and computation time.
# vec_finite_diff_2 <- function(f, x, h = 1e-8) {
#   to_vec <- function(x) as.numeric(unlist(x))
#
#   ufx <- to_vec(f(x))  # avoid duplicate evaluation
#   finite_deriv <- function(x, i) {
#     x[i] <- x[i] + h
#     (to_vec(f(x)) - ufx) / h
#   }
#
#   seq_along(x) %>%
#     purrr::map(~finite_deriv(x, i = .x)) %>%
#     dplyr::bind_cols()
# }

name_matrix <- function(x, input) {
  x <- as.matrix(x)
  colnames(x) <- paste("d", names(unlist(input)), sep = "_")
  rownames(x) <- rownames(x) %||% paste("d_output", seq(nrow(x)), sep = "_")
  x
}
