d_matrix_plus_matrix <- function(a, b) {
  a@dx + b@dx
}

d_matrix_minus_matrix <- function(a, b) {
  a@dx - b@dx
}

d_matrix_multiply_matrix <- function(a, b) {
  # Elementwise multiplication
  pa <- a@x
  pb <- b@x
  assertthat::assert_that(identical(dim(pa), dim(pb)))
  da <- a@dx
  db <- b@dx
  scale_rows_by_vector(db, as.numeric(pa)) +
    scale_rows_by_vector(da, as.numeric(pb))
}


d_matrix_divide_matrix <- function(a, b) {
  # Elementwise division
  assertthat::assert_that(identical(dim(a@x), dim(b@x)))
  pa <- a@x
  pb <- b@x
  da <- a@dx
  db <- b@dx
  scale_rows_by_vector(db, - as.numeric(pa) / as.numeric(pb)^2) +
    scale_rows_by_vector(da, 1 / as.numeric(pb))
}

scale_columns_by_vector <- function(m0, v0) {
  t(v0 * t(m0))
}

scale_rows_by_vector <- function(m0, v0) {
  v0 * m0
}
