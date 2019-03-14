# Frames
# Scalar by scalar: +, -, *, /
d_scalar_op_scalar <- function(a, b, d_op) {
  pa <- a@x
  pb <- b@x
  da <- a@dx
  db <- b@dx

  # Make sure we are indeed dealing with scalars
  assertthat::assert_that(length(pa) == 1)
  assertthat::assert_that(length(pb) == 1)

  d_op(da, db, pa, pb)
}

# Scalar by Matrix: +, -, *, /
d_scalar_op_matrix <- function(a, b, d_op) {
  pa <- a@x
  pb <- b@x
  da <- a@dx
  db <- b@dx
  assertthat::assert_that(length(pa) == 1)

  da <- matrix(rep(da, nrow(db)), nrow = nrow(db), byrow = T)
  d_op(da, db, pa, as.numeric(pb))
}

# Matrix by Scalar: +, -, *, /
d_matrix_op_scalar <- function(a, b, d_op) {
  pa <- a@x
  pb <- b@x
  da <- a@dx
  db <- b@dx
  assertthat::assert_that(length(pb) == 1)

  db <- matrix(rep(db, nrow(da)), nrow = nrow(da), byrow = T)
  d_op(da, db, as.numeric(pa), pb)
}

# Components
plus_fun <- function(dx, dy, ...) { dx + dy }
minus_fun <- function(dx, dy, ...) { dx - dy }
multiply_fun <- function(dx, dy, x, y) { dx * y + dy * x }
divide_fun <- function(dx, dy, x, y) { (y * dx - x * dy) / y^2 }

# Factory
d_scalar_plus_scalar <- purrr::partial(d_scalar_op_scalar, d_op = plus_fun)
d_scalar_minus_scalar <- purrr::partial(d_scalar_op_scalar, d_op = minus_fun)
d_scalar_multiply_scalar <- purrr::partial(d_scalar_op_scalar, d_op = multiply_fun)
d_scalar_divide_scalar <- purrr::partial(d_scalar_op_scalar, d_op = divide_fun)

d_scalar_plus_matrix <- purrr::partial(d_scalar_op_matrix, d_op = plus_fun)
d_scalar_minus_matrix <- purrr::partial(d_scalar_op_matrix, d_op = minus_fun)
d_scalar_multiply_matrix <- purrr::partial(d_scalar_op_matrix, d_op = multiply_fun)
d_scalar_divide_matrix <- purrr::partial(d_scalar_op_matrix, d_op = divide_fun)

d_matrix_plus_scalar <- purrr::partial(d_matrix_op_scalar, d_op = plus_fun)
d_matrix_minus_scalar <- purrr::partial(d_matrix_op_scalar, d_op = minus_fun)
d_matrix_multiply_scalar <- purrr::partial(d_matrix_op_scalar, d_op = multiply_fun)
d_matrix_divide_scalar <- purrr::partial(d_matrix_op_scalar, d_op = divide_fun)
