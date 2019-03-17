d_matrix_plus_matrix <- function(a, b) {
  a@dx + b@dx
}

d_matrix_minus_matrix <- function(a, b) {
  a@dx - b@dx
}

d_matrix_multiply_matrix <- function(a, b) {
  # Elementwise division
  pa <- a@x
  pb <- b@x
  assertthat::assert_that(identical(dim(pa), dim(pb)))
  da <- a@dx
  db <- b@dx
  d_res <- da
  for (i in 1:nrow(da)) {
    d_res[i, ] <- multiply_fun(da[i, ], db[i, ], pa[i], pb[i])
  }
  d_res
}

d_matrix_divide_matrix <- function(a, b) {
  # Elementwise division
  assertthat::assert_that(identical(dim(a@x), dim(b@x)))
  pa <- a@x
  pb <- b@x
  da <- a@dx
  db <- b@dx
  d_res <- da
  for (i in 1:nrow(da)) {
    d_res[i, ] <- divide_fun(da[i, ], db[i, ], pa[i], pb[i])
  }
  d_res
}
