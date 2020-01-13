# Test suite
lambda <- pryr::f


# Generate example inputs
generate_inputs <- function(config_ls, generator) {
  Map(generator, config_ls)
}


# Test a list of functions with a list of examples
test_fs <- function(fs, inputs, ctrl) {
  invisible(purrr::map(fs, function(f) {
    purrr::map(inputs, function(input) {
      test_AD_with(f, input, ctrl)
    })
  }))
}


# Compare Finite differencing with Automatic Differentiation
test_AD_with <- function(f, input, ctrl) {
  FD <- list(
    x = do.call(f, input),
    dx = finite_diff(f, at = input)
  )
  AD <- auto_diff(f, at = input) %>% {
    list(x = .@x, dx = .@dx)
  }
  compare_ls(FD, AD, ctrl)
}


# Compare two list of numeric arrays
compare_ls <- function(ls_1, ls_2, ctrl) {
  purrr::map2(
    ls_1, ls_2,
    ~do.call(compare, append(list(x = .x, y = .y), ctrl))
  )
}


# Compare two numeric arrays
compare <- function(x, y, display = T, err_fun = abs_err, epsilon,
                    summary_fun = max, note = "Maximum Error") {
  err <- err_fun(x, y)
  if (display) {
    x_vec <- as.numeric(x)
    y_vec <- as.numeric(y)
    print(cbind(x_vec, y_vec, err)[err > epsilon, ])
    cat(note, ": ", max(err), "\n")
  }
  testthat::expect_true(summary_fun(err) < epsilon)
}


abs_err <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  abs(x - y)
}


rel_err <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  scale <- pmax(abs(x), abs(y), 1e-8)  # avoid problem like: rel_err(0, 1e-18) = 1
  ifelse(scale == 0, 0, abs(x - y) / scale)
}
