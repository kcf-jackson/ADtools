testthat::context("Test solve, t, tcrossprod, crossprod, chol0, det")

set.seed(123)
# One-argument case (dual)
solve_f <- function(x) solve(x)  # change the argument name to match the other functions
fs <- list(solve_f, t, tcrossprod, crossprod, det)
inputs <- generate_inputs(2:5, lambda(list(x = 1 + randn(i, i))))
ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-5)
test_fs(fs, inputs, ctrl)


fs <- list(chol0)  # requires positive semi-definite input
inputs <- generate_inputs(2:5, lambda(list(x = 10 + crossprod(randn(i, i)))))
ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-5)
test_fs(fs, inputs, ctrl)


# Two-argument case (dual, dual) for tcrossprod, crossprod, solve
fs <- list(tcrossprod, crossprod, lambda(solve(x, y)))
inputs <- generate_inputs(
  2:5, lambda(list(x = 10 + randn(i, i), y = randn(i, i)))
)
ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-5)
test_fs(fs, inputs, ctrl)


# Two-argument case (dual, matrix) for tcrossprod, crossprod, solve
purrr::map(2:5, function(i) {
  m0 <- randn(i, i)
  fs <- list(
    lambda(x, solve(a = x, b = m0)),
    lambda(x, solve(a = 10 + m0, b = x)),
    lambda(x, tcrossprod(x, y = m0)),
    lambda(x, crossprod(x, y = m0))
  )
  inputs <- list( list(x = randn(i, i)) )
  ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-5)
  test_fs(fs, inputs, ctrl)
})
