# Test suite
lambda <- pryr::f


# Generate example inputs
generate_inputs <- function(config_ls, generator) {
  Map(generator, config_ls)
}


# Test a list of functions with a list of examples
test_fs <- function(fs, inputs, ctrl) {
  purrr::map(fs, function(f) {
    purrr::map(inputs, function(input) {
      test_AD_with(f, input, ctrl)
    })
  })
}


# Compare Finite differencing with Automatic Differentiation
test_AD_with <- function(f, input, ctrl) {
  FD <- list(
    x = do.call(f, input),
    dx = finite_diff(f, vary = input)
  )
  AD <- auto_diff(f, vary = input) %>% {
    list(x = parent_of(.), dx = get_deriv(.))
  }
  compare_ls(FD, AD, ctrl)
}


# Compare two list of numeric arrays
compare_ls <- function(ls_1, ls_2, ctrl) {
  purrr::map2(
    ls_1, ls_2,
    ~compare(.x, .y, ctrl$display, ctrl$err_fun, ctrl$epsilon)
  )
}


# Compare two numeric arrays
compare <- function(x, y, display = T, err_fun = abs_err, epsilon) {
  err <- err_fun(x, y)
  if (display) {
    x_vec <- as.numeric(x)
    y_vec <- as.numeric(y)
    print(cbind(x_vec, y_vec, err)[err > epsilon, ])
    cat("Maximum error: ", max(err), "\n")
  }
  testthat::expect_true(max(err) < epsilon)
}


abs_err <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  abs(x - y)
}


rel_err <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  scale <- max(abs(x), abs(y))
  if (scale == 0) return(0)
  abs(x - y) / scale
}
