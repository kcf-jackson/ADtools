d_matrix_plus_matrix <- function(a, b) {
  deriv_of(a) + deriv_of(b)
}

d_matrix_minus_matrix <- function(a, b) {
  deriv_of(a) - deriv_of(b)
}

d_matrix_multiply_matrix <- function(a, b) {
  # Elementwise division
  assertthat::assert_that(identical(dim(parent_of(a)), dim(parent_of(b))))
  pa <- parent_of(a)
  pb <- parent_of(b)
  da <- deriv_of(a)
  db <- deriv_of(b)
  d_res <- da
  for (i in 1:nrow(da)) {
    d_res[i, ] <- multiply_fun(da[i, ], db[i, ], pa[i], pb[i])
  }
  d_res
}

d_matrix_divide_matrix <- function(a, b) {
  # Elementwise division
  assertthat::assert_that(identical(dim(parent_of(a)), dim(parent_of(b))))
  pa <- parent_of(a)
  pb <- parent_of(b)
  da <- deriv_of(a)
  db <- deriv_of(b)
  d_res <- da
  for (i in 1:nrow(da)) {
    d_res[i, ] <- divide_fun(da[i, ], db[i, ], pa[i], pb[i])
  }
  d_res
}

d_matrix_prod <- function(a, b) {
  A <- parent_of(a)
  B <- parent_of(b)
  dA <- deriv_of(a)
  dB <- deriv_of(b)
  I_A <- Diagonal0(nrow(A))
  I_B <- Diagonal0(ncol(B))
  (I_B %x% A) %*% dB + (t(B) %x% I_A) %*% dA
  # I_x_B_times_C(A, dB) + A_x_I_times_C(t(B), dA)
}
