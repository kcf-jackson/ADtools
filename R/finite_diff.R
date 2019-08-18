#' Finite differencing
#' @param f A function of which the derivative is seeked.
#' @param vary A named list of variables; the variables to be varied.
#' @param fix A named list of variables; the variables to be fixed.
#' @param h The finite differencing parameter; the size of perturbation.
#' @param seed Seed; for pathwise derivative only.
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
#' @export
finite_diff <- function(f, wrt = NULL, at, h = 1e-8, seed) {
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
  
  vec_finite_diff(f_vec, x, h) %>%
    name_matrix(vary)
}

vec_finite_diff <- function(f, x, h = 1e-8) {
  to_vec <- function(x) as.numeric(unlist(x))
  
  ufx <- to_vec(f(x))  # avoid duplicate evaluation
  finite_deriv <- function(x) {
    (to_vec(f(x)) - ufx) / h
  }

  perturbate <- function(v, h) {
    purrr::map(seq_along(v), function(i) { v[i] <- v[i] + h; v })
  }

  perturbate(x, h) %>%
    purrr::map(finite_deriv) %>%
    purrr::bind_cols()
}

name_matrix <- function(x, input) {
  colnames(x) <- paste("d", names(unlist(input)), sep = "_")
  rownames(x) <- rownames(x) %||% paste("d_output", seq(nrow(x)), sep = "_")
  x
}
