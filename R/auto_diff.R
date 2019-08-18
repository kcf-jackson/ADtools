#' Automatic differentiation
#' @param f A function of which the derivative is seeked.
#' @param vary A named list of variables; the variables to be varied.
#' @param fix A named list of variables; the variables to be fixed.
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
  tidy_dx(AD_result, vary)
}


#' Converting a list of parameters into a list of dual numbers
#' @param vary A named list of parameters
#' @examples
#' \dontrun{
#' X <- randn(2, 2)
#' y <- randn(2)
#' dual_list(list(X = X, y = y))
#' }
#' @export
duals <- function(vary) {
  dims <- purrr::map(vary, length)
  ind <- seq_along(vary)
  purrr::map2(vary, ind, ~dual(.x, dims, .y))
}


#' Add rownames and colnames to the dual component
#' @param x_dual A dual number
#' @param vary A named list of parameters
tidy_dx <- function(x_dual, vary) {
  make_colnames <- function(x) {
	magrittr::set_colnames(
		x = x, 
		purrr::map2_chr(
			x = names(vary),
			y = purrr::map(vary, length),
			f = ~paste("d_", .x, .y, sep = "")
		)
	)
  }
  make_rownames <- function(x) {
	magrittr::set_rownames(
		x = x, 
		paste("d_output", seq(nrow(x)), sep = "_")
	)
  }

  x_dual@dx %>%
    magrittr::make_colnames()
    magrittr::make_rownames()
}
