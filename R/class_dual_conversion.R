#' Automatic differentiation
#' @param f A function of which the derivative is seeked.
#' @param vary A named list of variables; the variables to be varied.
#' @param fix A named list of variables; the variables to be fixed.
#' @examples
#' f <- function(y, X, beta) { y - X %*% beta }
#' auto_diff(f,
#'   vary = list(beta = c(5,6)),
#'   fix = list(X = matrix(1:4, 2, 2), y = c(2,3))
#' )
#'
#' g <- function(X, Y) { X %*% Y }
#' X <- randn(2, 2)
#' Y <- randn(2, 2)
#' auto_diff(g, vary = list(X = X, Y = Y))
#'
#' @export
auto_diff <- function(f, vary, fix = NULL) {
  do.call(f, append(dual_list(vary), fix))
}


#' Converting a list of parameters into a list of dual numbers
#' @param ls_params A named list of parameters
#' @examples
#' \dontrun{
#' X <- randn(2, 2)
#' y <- randn(2)
#' dual_list(list(X = X, y = y))
#' }
#' @export
dual_list <- function(ls_params) {
  list0_len <- purrr::map_dbl(ls_params, length)
  purrr::map(
    seq_along(ls_params),
    ~dual(ls_params[[.x]], list0_len, .x)
  ) %>%
    setNames(names(ls_params))
}


#' Single-level destructuring
#' @description Strict usage: list(x, y) %<-% list(a, b, ...).
#' LHS must start with 'list'; RHS must be a list;
#' RHS must have at least as many elements as LHS.
#' @param lhs An expression; variables to be assigned.
#' @param value A list; values to assign.
#' @note Multiple-level is avoided intentionally, i.e. nested destructuring
#' is not allowed.
#' @examples
#' \dontrun{
#' list(x, y) %<-% list(1, 2)
#' list(x, y) %<-% list(1, 2, 3)
#' }
#' @export
`%<-%` <- function(lhs, value) {
  lhs_chr <- Map(deparse, substitute(lhs))  # Parse LHS into characters
  variables <- lhs_chr[-1]

  if (lhs_chr[[1]] != "list") stop("LHS must start with 'list'.")
  if (!is.list(value)) stop("RHS must be a list.")
  if (length(value) < length(variables))
    stop("RHS must have at least as many elements as LHS.")

  for (i in seq_along(variables)) {
    assign(variables[[i]], value[[i]], envir = parent.frame())
  }
  invisible(NULL)
}
