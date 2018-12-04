#=========================================================================
# Class "dual"
#-------------------------------------------------------------------------
setClassUnion("array_or_numeric", c("array", "numeric"))

is_matrix <- function(x) {
  (attr(class(x), "package") == "Matrix") || class(x) == "matrix"
}

check_dual <- function(object) {
  x <- object@x
  dx <- object@dx
  # x_check <- is.array(x) || is.numeric(x)
  dx_check <- nrow(dx) == length(x)   # match dimension
  is_matrix(x) && is_matrix(dx) && dx_check
}


#' S4 class "dual"
#'
#' @description This class attaches a dual component to a number / an array.
#' @slot x Numeric matrix.
#' @slot dx matrix; also accepts any matrix classes from the "Matrix" package.
#' @import methods
#' @examples
#' a <- new("dual", x = randn(2,2), dx = 0)
#' b <- new("dual", x = randn(2,2), dx = init_dx(list(c(4, 1), c(4, 2)), 2))
setClass(
  "dual",
  representation(
    x = "ANY",
    dx = "ANY",
    param = "list"
  ),
  validity = check_dual
)


#' Dual number extractor
#' @param x A dual number
#' @export
parent_of <- function(x) { x@x }

#' Dual number extractor
#' @param x A dual number
#' @export
deriv_of <- function(x) { x@dx }

#' Dual number extractor
#' @param x A dual number
#' @export
param_of <- function(x) { x@param }


#' Dual number constructor
#' @param x The object to be converted to a dual number.
#' @param param_dim The dimension of the dual component to be attached.
#' @param ind Integer; the index of the parameter identical to x.
#' @export
dual <- function(x, param_dim, ind) {
  param_name <- names(param_dim)
  if (is.null(param_name)) {
    param_name <- paste0("V", seq_along(param_dim))
  }
  param_dim <- unlist(param_dim)
  new("dual",
      x = x,
      dx = init_dx(length(x), param_dim, ind),
      param = setNames(dim_to_col_range(param_dim), param_name))
}

dim_to_col_range <- function(dim0) {
  end <- cumsum(dim0)
  start <- c(1, head(end, -1) + 1)
  map_row(cbind(start, end), identity)
}


#' Extract derivative
#' @param x A dual number
#' @param wrt A character string; the parameter name.
#' @export
get_deriv <- function(x, wrt) {
  rng <- param_of(x)[[wrt]]
  as.matrix(deriv_of(x))[, rng[1]:rng[2], drop = F]
}


#' Initialise the dual component, i.e. the storage matrices for derivatives
#'
#' @param num_dim A number; the length of the numerator of a derivative.
#' @param denom_dim Numeric vector; the length of the denominators of a derivative.
#' @param num_ind An integer that indicates the location of the derivative
#' corresponding to x; input -1 if none.
#'
#' @examples
#' init_dx(4, c(2,5,1), -1)
#' init_dx(4, c(2,5,1), 1)
#' init_dx(4, c(2,5,1), 2)
#' init_dx(4, c(2,5,1), 3)
#'
#' @export
init_dx <- function(num_dim, denom_dim, num_ind) {
  deriv <- purrr::map(denom_dim, ~zero_matrix0(num_dim, .x))
  if (num_ind != -1) {
    deriv[[num_ind]] <- one_matrix0(num_dim, denom_dim[num_ind])
  }
  do.call(cbind, deriv)
}
