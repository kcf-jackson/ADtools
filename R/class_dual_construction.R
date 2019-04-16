#' Dual number constructor
#' @param x The object to be converted to a dual number.
#' @param param_dim A named list. The dimension of the dual component to be attached.
#' @param ind Integer; the index in `param_dim` corresponding to `x`. Use `-1` if it is not applicable.
#' @examples
#' # Suppose X is a 2 x 2 matrix, Y is a 3 x 1 vector, Z is a 2 x 3 matrix, and
#' # we wish to attach dual components {dX, dY, dZ} to X.
#' dual(randn(2, 2), list(X = 4, Y = 3, Z = 6), ind = 1)
#' @export
dual <- function(x, param_dim, ind) {
  x <- cast_vector_into_matrix(x)
  dim0 <- unlist(param_dim)
  param_name <- if_null_then(names(param_dim), paste0("V", seq_along(param_dim)))
  new("dual",
    x = x,
    dx = init_dx(length(x), dim0, ind),
    param = setNames(dim_to_col_range(dim0), param_name)
  )
}

cast_vector_into_matrix <- function(x) {
  if (is.vector(x) && !is_scalar(x)) {
    return(as.matrix(x))
  }
  return(x)
}

dim_to_col_range <- function(dim0) {
  end <- cumsum(dim0)
  start <- c(1, head(end, -1) + 1)
  purrr::map2(start, end, ~c(start = .x, end = .y))
}

col_range_to_dim <- function(dim0) {
  purrr::map_dbl(dim0, ~diff(.x) + 1)
}


#' Extract derivative
#' @param x A dual number
#' @param wrt A character vector; the parameter name.
#' @export
get_deriv <- function(x, wrt) {
  if (missing(wrt)) wrt <- names(param_of(x))
  get_range <- function(list0) {
    mapreduce(x = list0, f = ~seq(.x[1], .x[2]), g = c)
  }
  make_colnames <- function(list0) {
    map2reduce(
      x = names(list0),
      y = purrr::map(list0, ~ seq(diff(.x) + 1)),
      f = ~paste("d_", .x, .y, sep = ""),
      g = c
    )
  }
  set_rownames <- function(x) {
    rownames(x) <- paste("d_output", seq(nrow(x)), sep = "_")
    x
  }

  var_of_interest <- param_of(x)[wrt]
  rng <- get_range(var_of_interest)
  d_var <- make_colnames(var_of_interest)

  as.matrix(x@dx)[, rng, drop = F] %>%
    magrittr::set_colnames(d_var) %>%
    set_rownames()
}


#' Initialise the dual component, i.e. the storage matrices for derivatives
#' @param num_dim A number; the length of the numerator of a derivative.
#' @param denom_dim Numeric vector; the length of the denominators of a derivative.
#' @param num_ind An integer that indicates the location of the derivative
#' corresponding to x; input -1 if none.
#' @examples
#' \dontrun{
#' init_dx(4, c(2, 5, 1), -1)
#' init_dx(4, c(2, 5, 1), 1)
#' init_dx(4, c(2, 5, 1), 2)
#' init_dx(4, c(2, 5, 1), 3)
#' }
#' @keywords internal
init_dx <- function(num_dim, denom_dim, num_ind) {
  deriv <- purrr::map(denom_dim, ~zero_matrix0(num_dim, .x))
  if (num_ind != -1) {
    deriv[[num_ind]] <- one_matrix0(num_dim, denom_dim[num_ind])
  }
  do.call(cbind, deriv)
}
