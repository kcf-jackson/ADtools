d_matrix_plus_matrix <- function(a, b) {
  a@dx + b@dx
}

d_matrix_minus_matrix <- function(a, b) {
  a@dx - b@dx
}

d_matrix_multiply_matrix <- function(a, b) {
  # Elementwise multiplication
  pa <- as.numeric(a@x)
  pb <- as.numeric(b@x)
  assertthat::are_equal(length(pa), length(pb))
  da <- a@dx
  db <- b@dx
  diag_v0_times_m0(pa, db) + diag_v0_times_m0(pb, da)
}

d_matrix_divide_matrix <- function(a, b) {
  # Elementwise division
  assertthat::assert_that(identical(dim(a@x), dim(b@x)))
  pa <- as.numeric(a@x)
  pb <- as.numeric(b@x)
  da <- a@dx
  db <- b@dx
  diag_v0_times_m0(-pa / pb^2, db) + diag_v0_times_m0(1 / pb, da)
}


# Factory -----------------------------------------------------------------
d_scalar_plus_scalar     <- function(a, b) d_scalar_op_scalar(a, b, d_op = plus_fun)
d_scalar_minus_scalar    <- function(a, b) d_scalar_op_scalar(a, b, d_op = minus_fun)
d_scalar_multiply_scalar <- function(a, b) d_scalar_op_scalar(a, b, d_op = multiply_fun)
d_scalar_divide_scalar   <- function(a, b) d_scalar_op_scalar(a, b, d_op = divide_fun)

d_scalar_plus_matrix     <- function(a, b) d_scalar_op_matrix(a, b, d_op = plus_fun)
d_scalar_minus_matrix    <- function(a, b) d_scalar_op_matrix(a, b, d_op = minus_fun)
d_scalar_multiply_matrix <- function(a, b) d_scalar_op_matrix(a, b, d_op = multiply_fun)
d_scalar_divide_matrix   <- function(a, b) d_scalar_op_matrix(a, b, d_op = divide_fun)

d_matrix_plus_scalar     <- function(a, b) d_matrix_op_scalar(a, b, d_op = plus_fun)
d_matrix_minus_scalar    <- function(a, b) d_matrix_op_scalar(a, b, d_op = minus_fun)
d_matrix_multiply_scalar <- function(a, b) d_matrix_op_scalar(a, b, d_op = multiply_fun)
d_matrix_divide_scalar   <- function(a, b) d_matrix_op_scalar(a, b, d_op = divide_fun)

# Frames ------------------------------------------------------------------
# Scalar by Scalar: +, -, *, /
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

  da <- matrix(rep(da, NROW(db)), nrow = NROW(db), byrow = TRUE)
  d_op(da, db, pa, as.numeric(pb))
}

# Matrix by Scalar: +, -, *, /
d_matrix_op_scalar <- function(a, b, d_op) {
  pa <- a@x
  pb <- b@x
  da <- a@dx
  db <- b@dx
  assertthat::assert_that(length(pb) == 1)

  db <- matrix(rep(db, NROW(da)), nrow = NROW(da), byrow = TRUE)
  d_op(da, db, as.numeric(pa), pb)
}

# Components --------------------------------------------------------------
plus_fun     <- function(dx, dy, ...) { dx + dy }
minus_fun    <- function(dx, dy, ...) { dx - dy }
multiply_fun <- function(dx, dy, x, y) { dx * y + dy * x }
divide_fun   <- function(dx, dy, x, y) { (y * dx - x * dy) / y^2 }
