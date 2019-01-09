finite_diff <- function(f, vary, fix = NULL, h = 1e-8, seed) {  # add fix
  has_seed <- !missing(seed)
  x <- unlist(vary)
  f_vec <- function(vec0) {
    if (has_seed) set.seed(seed)
    do.call(f, append(relist(vec0, vary), fix))    # add append
  }
  vec_finite_diff(f_vec, x, h) %>%
    magrittr::set_colnames(flatten_name(vary))
}

vec_finite_diff <- function(f, x, h = 1e-8) {
  ufx <- unlist(f(x))  # avoid duplicate evaluation
  finite_deriv <- . %>% { (unlist(f(.)) - ufx) / h }
  perturbate <- function(v, h) {
    purrr::map(seq_along(v), function(i) { v[i] <- v[i] + h; v })
  }
  res <- map_then_call(perturbate(x, h), finite_deriv, cbind)
  if (is.null(rownames(res))) res <- add_rownames(res)
  res
}

flatten_name <- function(list0) {
  paste("d", names(unlist(list0)), sep = "_")
}

add_rownames <- function(x) {
  rownames(x) <- paste("d_output", seq(nrow(x)), sep = "_")
  x
}

map_then_call <- function(x, f, g) {
  x %>% purrr::map(f) %>% do.call(g, .)
}
