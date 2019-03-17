# Other matrix functions
d_solve <- function(a, inv_a) {
  da <- a@dx
  -(t(inv_a) %x% inv_a) %*% da
}


d_transpose <- function(a) {
  X <- a@x
  dX <- a@dx
  K_nq <- commutation_matrix0(nrow(X), ncol(X))
  K_nq %*% dX
}


d_XXT <- function(a) {
  X <- a@x
  dX <- a@dx

  n <- nrow(X)
  I_n <- Diagonal0(n)
  I_nn <- Diagonal0(n^2)
  K_nn <- commutation_matrix0(n, n)
  ((I_nn + K_nn) %*% (X %x% I_n)) %*% dX
}


d_chol <- function(L, a) {
  # LL^T = A
  dA <- a@dx

  n <- nrow(L)
  I_n <- Diagonal0(n)
  I_nn <- Diagonal0(n^2)
  K_nn <- commutation_matrix0(n, n)
  E_n <- elimination_matrix0(n)
  D_n <- Matrix::t(E_n)

  (D_n %*% solve(E_n %*% (I_nn + K_nn) %*% (L %x% I_n) %*% D_n, E_n)) %*% dA
}
