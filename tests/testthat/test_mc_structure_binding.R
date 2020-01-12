testthat::context("Test binding two dual numbers")

# Functions without S3 declarations need to be wrapped
wrap2 <- function(fs) purrr::map(fs, ~lambda(.x(x, y)))


testthat::test_that("test cbind2, rbind2", {
  g <- function(x, y) rbind2(x * y, x + y)  # additional check
  fs <- list(cbind2, rbind2, g)
  inputs <- generate_inputs(
    config_ls = 2:10,
    generator = lambda(list(x = randn(i,i), y = randn(i,i)))
  )
  ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-7)
  test_fs(wrap2(fs), inputs, ctrl)
})


testthat::test_that("test cbind2, rbind2", {
  ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-8)
  purrr::map(2:10, function(i) {
  	fs <- list(cbind2, rbind2)
  	m0 <- randn(i, i)
  	# dual - numeric
  	gs <- purrr::map(fs, ~lambda(x, .x(x, m0)))
  	inputs <- list(list(x = randn(i, i)))
  	test_fs(gs, inputs, ctrl)
  	# numeric - dual
  	gs <- purrr::map(fs, ~lambda(y, .x(m0, y)))
  	inputs <- list(list(y = randn(i, i)))
  	test_fs(gs, inputs, ctrl)
  })
})
