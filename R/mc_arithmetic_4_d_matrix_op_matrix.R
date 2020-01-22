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
