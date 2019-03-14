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

  if (is.null(dim(x)) && length(x) > 1) x <- as.matrix(x)
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

col_range_to_dim <- function(dim0) {
  dim0 %>% purrr::map_dbl(diff) %>% magrittr::add(1)
}


#' Extract derivative
#' @param x A dual number
#' @param wrt A character vector; the parameter name.
#' @export
get_deriv <- function(x, wrt) {
  if (missing(wrt)) wrt <- names(param_of(x))
  get_range <- function(list0) {
    list0 %>%
      purrr::map(~seq(.x[1], .x[2])) %>%
      do.call(c, .)
  }
  make_colnames <- function(list0) {
    purrr::map2(
      .x = names(list0),
      .y = purrr::map(list0, ~seq(diff(.x) + 1)),
      ~paste("d_", .x, .y, sep = "")
    ) %>%
      do.call(c, .)
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
#' init_dx(4, c(2,5,1), -1)
#' init_dx(4, c(2,5,1), 1)
#' init_dx(4, c(2,5,1), 2)
#' init_dx(4, c(2,5,1), 3)
#' }
#' @keywords interal
init_dx <- function(num_dim, denom_dim, num_ind) {
  deriv <- purrr::map(denom_dim, ~zero_matrix0(num_dim, .x))
  if (num_ind != -1) {
    deriv[[num_ind]] <- one_matrix0(num_dim, denom_dim[num_ind])
  }
  do.call(cbind, deriv)
}
