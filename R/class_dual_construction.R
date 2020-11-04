#' Dual number constructor
#'
#' @param x The object to be converted to a dual number.
#' @param param_dim A named list. The dimension of the dual component to be attached.
#' @param ind Integer; the index in `param_dim` corresponding to `x`. Use `-1` if it is not applicable.
#'
#' @return A dual number with components "x" and "dx". The first gives the value of `f`, and the
#' second gives the derivative of `f`.
#'
#' @examples
#' # Suppose X is a 2 x 2 matrix, Y is a 3 x 1 vector, Z is a 2 x 3 matrix, and
#' # we wish to attach dual components {dX, dY, dZ} to X.
#' dual(randn(2, 2), list(X = 4, Y = 3, Z = 6), ind = 1)
#'
#' @export
dual <- function(x, param_dim, ind) {
  x <- cast_vector_into_matrix(x)
  dx <- init_dx(length(x), param_dim, ind)
  new("dual", x = x, dx = dx)
}

cast_vector_into_matrix <- function(x) {
  if (is_vector(x)) as.matrix(x) else x
}


#' Initialise the dual component
#'
#' @param num_dim A number; the length of the numerator of a derivative.
#' @param denom_dim Numeric vector; the length of the denominators of a derivative.
#' @param num_ind An integer that indicates the location of the derivative
#' corresponding to x; input -1 if none.
#'
#' @return A numeric matrix; the dual component.
#'
# #' @examples
# #' init_dx(4, c(2, 5, 1), -1)
# #' init_dx(4, c(2, 5, 1), 1)
# #' init_dx(4, c(2, 5, 1), 2)
# #' init_dx(4, c(2, 5, 1), 3)
# #'
#' @keywords internal
init_dx <- function(num_dim, denom_dim, num_ind) {
  deriv <- purrr::map(denom_dim, ~zero_matrix0(num_dim, .x))
  if (num_ind != -1) {
    deriv[[num_ind]] <- band_matrix0(num_dim, denom_dim[[num_ind]])
  }
  do.call(cbind, deriv)
}
