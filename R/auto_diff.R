#' Automatic differentiation
#'
#' @param f A function of which the derivative is sought.
#' @param wrt A character vector; the name of the variables to differentiate with respect to.
#' @param at A named list of variables; the point at which the derivative is evaluated.
#'
#' @return A dual number with components "x" and "dx". The first gives the value of `f`, and the
#' second gives the derivative of `f`.
#'
#' @examples
#' f <- function(y, X, beta) { y - X %*% beta }
#' auto_diff(
#'   f, wrt = "beta",
#'   at = list(beta = c(5,6), X = matrix(1:4, 2, 2), y = c(2,3))
#' )
#'
#' g <- function(X, Y) { X %*% Y }
#' X <- randn(2, 2)
#' Y <- randn(2, 2)
#' auto_diff(g, at = list(X = X, Y = Y))
#'
#' @export
auto_diff <- function(f, wrt = NULL, at) {
  f_args <- formalArgs(f)
  wrt <- wrt %||% f_args
  vary <- at[wrt]
  fix <- at[setdiff(f_args, wrt)]

  dual_args <- append(duals(vary), fix)
  AD_result <- do.call(f, dual_args)
  recursive_tidy(AD_result, vary)
}


#' Converting a list of parameters into a list of dual numbers
#'
#' @param vary A named list of parameters
#'
#' @return A named list of dual numbers.
#'
#' @examples
#' X <- randn(2, 2)
#' y <- rnorm(2)
#' duals(list(X = X, y = y))
#'
#' @export
duals <- function(vary) {
  dims <- purrr::map(vary, length)
  ind <- seq_along(vary)
  purrr::map2(vary, ind, ~dual(.x, dims, .y))
}


# An interface to call tidy_dx recursively
recursive_tidy <- function(x, vary) {
  if (is.list(x)) {
    purrr::map(x, recursive_tidy, vary = vary)
  } else {
    if (class(x) == "dual") tidy_dx(x, vary) else x
  }
}


#' Add rownames and colnames to the dual component
#'
#' @param x_dual A dual number
#' @param vary A named list of parameters
tidy_dx <- function(x_dual, vary) {
  make_colnames <- function(x) {
  	magrittr::set_colnames(
  		x = x,
  		unlist(purrr::map2(
  			.x = names(vary),
  			.y = purrr::map(vary, seq_along),
  			.f = ~paste("d_", .x, .y, sep = "")
  		))
  	)
  }
  make_rownames <- function(x) {
  	magrittr::set_rownames(
  		x = x,
  		paste("d_output", seq(nrow(x)), sep = "_")
  	)
  }

  x_dual@dx <- as.matrix(x_dual@dx) %>%
    make_colnames() %>%
    make_rownames()
  x_dual
}
