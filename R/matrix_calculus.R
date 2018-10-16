#=========================================================================
# Differential calculus
#-------------------------------------------------------------------------
d_sum <- function(a, b) {
  da <- deriv_of(a)
  db <- deriv_of(b)
  if (is_zero(da)) return(db)
  if (is_zero(db)) return(da)
  da + db
}

d_minus <- function(a, b) {
  da <- deriv_of(a)
  db <- deriv_of(b)
  if (is_zero(da)) return(-db)
  if (is_zero(db)) return(da)
  da - db
}

# d_scalar_prod <- function(a, b) {
#   x <- parent_of(a)
#   y <- parent_of(b)
#   dx <- deriv_of(a)
#   dy <- deriv_of(b)
#   if (is_zero(dx) || is_zero(dy)) return(0)
#   x * dy + y * dx
# }
# d_divide <- function(a, b) {  }

d_matrix_prod <- function(a, b) {
  A <- parent_of(a)
  B <- parent_of(b)
  dA <- deriv_of(a)
  dB <- deriv_of(b)
  if (is_zero(dA) || is_zero(dB)) return(0)
  I_A <- Matrix::Diagonal(nrow(A))
  I_B <- Matrix::Diagonal(ncol(B))
  (I_B %x% A) %*% dB + (t(B) %x% I_A) %*% dA
  # I_x_B_times_C(A, dB) + A_x_I_times_C(t(B), dA)
}

d_solve <- function(a, inv_a) {
  da <- deriv_of(a)
  - (t(inv_a) %x% inv_a) %*% da
}

d_transpose <- function(a) {
  X <- parent_of(a)
  dX <- deriv_of(a)
  K_nq <- commutation_matrix(nrow(X), ncol(X))
  K_nq %*% dX
}

d_XXT <- function(a) {
  X <- parent_of(a)
  dX <- deriv_of(a)

  n <- nrow(X)
  I_n <- Matrix::Diagonal(n)
  I_nn <- Matrix::Diagonal(n^2)
  K_nn <- commutation_matrix(n, n)
  ((I_nn + K_nn) %*% (X %x% I_n)) %*% dX
}

d_chol <- function(L, a) {
  # LL^T = A
  dA <- deriv_of(a)

  n <- nrow(L)
  I_n <- Matrix::Diagonal(n)
  I_nn <- Matrix::Diagonal(n^2)
  K_nn <- commutation_matrix(n, n)
  E_n <- elimination_matrix(n)
  D_n <- Matrix::t(E_n)

  (D_n %*% solve(E_n %*% (I_nn + K_nn) %*% (L %x% I_n) %*% D_n) %*% E_n) %*% dA
}

d_kronecker <- function(a, b) {
  A <- parent_of(a)
  B <- parent_of(b)
  dA <- deriv_of(a)
  dB <- deriv_of(b)

  m <- nrow(A)
  n <- ncol(A)
  p <- nrow(B)
  q <- ncol(B)
  I_n <- Matrix::Diagonal(n)
  K_qm <- commutation_matrix(q, m)
  I_p <- Matrix::Diagonal(p)
  I_mn <- Matrix::Diagonal(m*n)
  I_pq <- Matrix::Diagonal(p*q)

  (I_n %x% K_qm %x% I_p) %*%
    (as.numeric(A) %x% dB + dA %x% as.numeric(B))
}


d_normal <- function() {

}
# d_gamma

is_zero <- function(x) { is.numeric(x) && (x == 0) }
