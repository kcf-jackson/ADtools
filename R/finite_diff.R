# Finite differencing
finite_diff <- function(f, x, h = 1e-8) {
  perturbate <- function(v, h) {
    purrr::map(seq_along(v), function(i) { v[i] <- v[i] + h; v })
  }
  diff_fun <- . %>% { as.numeric((f(.) - f(x)) / h) }

  perturbate(x, h) %>%
    purrr::map(diff_fun) %>%
    do.call(cbind, .)
}
