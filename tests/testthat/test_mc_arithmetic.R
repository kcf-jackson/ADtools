testthat::context("Test dual number arithmetic")

set.seed(123)
`%++%` <- append

#=======================================================================
# dual matrix AND dual matrix
f2 <- function(e1, e2) { e1 %x% e2 }
fs <- list(`+`, `-`, `*`, `/`, f2)
inputs <- generate_inputs(
  2:12,
  lambda(list(e1 = randn(i, i), e2 = 10 + randn(i, i))) # add 10 to avoid dividing by near-zero
)
ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-6)
test_fs(fs, inputs, ctrl)

# Matrix product is like taking squares; precision is often lost with
# finite-difference method, so it is tested separately.
f <- function(e1, e2) { e1 %*% e2 }
fs <- list(f)
inputs <- generate_inputs(
  2:12,
  lambda(list(e1 = 1 + randn(i, i), e2 = 1 + randn(i, i)))
)
ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-6)
test_fs(fs, inputs, ctrl)


#=======================================================================
# dual scalar AND dual matrix
fs <- list(`+`, `-`, `*`, `/`)
inputs <- generate_inputs(2:15, lambda(list(e1 = i, e2 = 10 + randn(i, i)))) %++%
    generate_inputs(2:15, lambda(list(e1 = 10 + randn(i, i), e2 = i))) %++%
    generate_inputs(2:15, lambda(list(e1 = i, e2 = i)))
ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-6)
test_fs(fs, inputs, ctrl)


#=======================================================================
# numeric matrix AND dual matrix
f <- function(e1, e2) { e1 %*% e2 }
f2 <- function(e1, e2) { e1 %x% e2 }
fs <- list(`+`, `-`, `*`, `/`, f2)
ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-6)
purrr::map(2:15, function(i) {
  m0 <- randn(i, i)
  gs <- purrr::map(fs, function(f) lambda(e2, f(m0, e2)))
  inputs <- list( list(e2 = 10 + randn(i, i)) )
  test_fs(gs, inputs, ctrl)
})
purrr::map(2:15, function(i) {
  m0 <- 10 + randn(i, i)
  gs <- purrr::map(fs, function(f) lambda(e1, f(e1, m0)))
  inputs <- list( list(e1 = randn(i, i)) )
  test_fs(gs, inputs, ctrl)
})


# Matrix multiplication
fs <- list(f)
ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-6)
purrr::map(2:15, function(i) {
  m0 <- randn(i, i)
  gs <- purrr::map(fs, function(f) lambda(e2, f(m0, e2)))
  inputs <- list( list(e2 = randn(i, i)) )
  test_fs(gs, inputs, ctrl)
})
purrr::map(2:15, function(i) {
  m0 <- randn(i, i)
  gs <- purrr::map(fs, function(f) lambda(e1, f(e1, m0)))
  inputs <- list( list(e1 = randn(i, i)) )
  test_fs(gs, inputs, ctrl)
})


#=======================================================================
# numeric scalar AND dual matrix
fs <- list(`+`, `-`, `*`, `/`)
ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-6)
purrr::map(2:15, function(i) {
  m0 <- i
  gs <- purrr::map(fs, function(f) lambda(e2, f(m0, e2)))
  inputs <- list( list(e2 = 10 + randn(i, i)) )
  test_fs(gs, inputs, ctrl)
})
purrr::map(2:15, function(i) {
  m0 <- 10 + i
  gs <- purrr::map(fs, function(f) lambda(e1, f(e1, m0)))
  inputs <- list( list(e1 = randn(i, i)) )
  test_fs(gs, inputs, ctrl)
})


#=======================================================================
# numeric matrix AND dual scalar
fs <- list(`+`, `-`, `*`, `/`)
ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-6)
purrr::map(2:15, function(i) {
  m0 <- randn(i, i)
  gs <- purrr::map(fs, function(f) lambda(e2, f(m0, e2)))
  inputs <- list( list(e2 = 10 + runif(1, max = i)) )
  test_fs(gs, inputs, ctrl)
})
purrr::map(2:15, function(i) {
  m0 <- 10 + randn(i, i)
  gs <- purrr::map(fs, function(f) lambda(e1, f(e1, m0)))
  inputs <- list( list(e1 = runif(1, max = i)) )
  test_fs(gs, inputs, ctrl)
})


#=======================================================================
# Power
ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-5) # lower precision for high exponents
purrr::map(2:15, function(i) {
  # dual matrix - dual scalar
  # base must have only positive entries
  gs <- list(`^`)
  inputs <- list( list(e1 = randu(i, i, min = 0.1), e2 = i) )
  test_fs(gs, inputs, ctrl)

  # # numeric matrix - dual scalar
  m0 <- randu(i, i, min = 0.1)   # base must have only positive entries
  gs <- list( lambda(e2, m0 ^ e2) )
  inputs <- list( list(e2 = i) )
  test_fs(gs, inputs, ctrl)

  # dual matrix - numeric scalar
  gs <- list(lambda(e1, e1 ^ i))
  inputs <- list( list(e1 = randu(i, i, min = 0.1)) )
  test_fs(gs, inputs, ctrl)
})
