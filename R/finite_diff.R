#' Finite differencing
#' @param f A function of which the derivative is seeked.
#' @param vary A named list of variables; the variables to be varied.
#' @param fix A named list of variables; the variables to be fixed.
#' @param h The finite differencing parameter; the size of perturbation.
#' @param seed Seed; for pathwise derivative only.
#' @examples
#' f <- function(y, X, beta) { y - X %*% beta }
#' finite_diff(f,
#'   vary = list(beta = c(5,6)),
#'   fix = list(X = matrix(1:4, 2, 2), y = c(2,3))
#' )
#'
#' g <- function(X, Y) { X %*% Y }
#' finite_diff(g, vary = list(X = randn(2, 2), Y = randn(2, 2)))
#'
#' @export
finite_diff <- function(f, vary, fix = NULL, h = 1e-8, seed) {  # add fix
  dim_len <- if_null_then(dim, length)
  dim_ls <- purrr::map(vary, dim_len)
  has_seed <- !missing(seed)

  x <- list_to_vec(vary)
  f_vec <- function(vec0) {
    if (has_seed) set.seed(seed)
    do.call(f, append(vec_to_list(vec0, dim_ls), fix))    # add append
  }

  vec_finite_diff(f_vec, x, h) %>%
    magrittr::set_colnames(flatten_name(vary)) %>%
    add_rownames()
}


#' Finite differencing (Core)
#' @param f A function of which the derivative is seeked.
#' @param x The point at which the derivative is seeked.
#' @param h The finite differencing parameter; the size of perturbation.
#' @examples
#' \dontrun{
#' f <- function(x) { x[1] + 2* x[2] }
#' vec_finite_diff(f, c(2,3))
#' }
#' @keywords internal
vec_finite_diff <- function(f, x, h = 1e-8) {
  fx <- f(x)  # avoid duplicate evaluation
  finite_deriv <- . %>% { as.numeric((f(.) - fx) / h) }
  perturbate <- function(v, h) {
    purrr::map(seq_along(v), function(i) { v[i] <- v[i] + h; v })
  }

  map_then_call(perturbate(x, h), finite_deriv, cbind)
}


# ===== Converting between list and vector =====
list_to_vec <- function(list0) {
  map_then_call(list0, as.numeric, c)
}

vec_to_list <- function(vec0, ls_dim) {
  structure0 <- function(x, dim) {
    # R deprecates recycling of length 1 array, so extra handling is needed.
    if (isTRUE(dim == 1)) return(x)
    structure(x, dim = dim)
  }
  # extract elements according to ls_dim
  vec0_index <- cumsum(purrr::map_dbl(ls_dim, prod))
  ls_vec <- map_diff(c(0, vec0_index), ~vec0[seq(.x + 1, .y)])
  # reshape elements according to ls_dim
  res <- purrr::map2(ls_vec, ls_dim, ~structure0(.x, dim = .y))
  setNames(res, names(ls_dim))
}

flatten_name <- function(list0) {
  paste("d", names(list_to_vec(list0)), sep = "_")
}

add_rownames <- function(x) {
  rownames(x) <- paste("d_output", seq(nrow(x)), sep = "_")
  x
}

# ===== High-order helper functions =====
map_then_call <- function(x, f, g) {
  x %>% purrr::map(f) %>% do.call(g, .)
}

# Apply f1, if null, then apply f2
if_null_then <- function(f1, f2) {
  function(x) {
    f1x <- f1(x)
    if (is.null(f1x)) return(f2(x))
    return(f1x)
  }
}

# Apply f to every two consecutive items in a vector
map_diff <- function(vec0, f) {
  purrr::map2(head(vec0, -1), tail(vec0, -1), f)
}
