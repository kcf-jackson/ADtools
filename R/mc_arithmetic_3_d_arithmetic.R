# Matrix +, -, *, /
d_op_dispatch <- function(a, b, funs) {
  if (is_scalar(a@x)) {
    if (is_scalar(b@x)) {
      return(funs[[1]](a, b))
    } else {
      return(funs[[2]](a, b))
    }
  } else {
    if (is_scalar(b@x)) {
      return(funs[[3]](a, b))
    } else {
      return(funs[[4]](a, b))
    }
  }
}

d_sum <- function(a, b) {
  d_op_dispatch(a, b, list(
    d_scalar_plus_scalar, d_scalar_plus_matrix,
    d_matrix_plus_scalar, d_matrix_plus_matrix
  ))
}

d_minus <- function(a, b) {
  d_op_dispatch(a, b, list(
    d_scalar_minus_scalar, d_scalar_minus_matrix,
    d_matrix_minus_scalar, d_matrix_minus_matrix
  ))
}

d_scalar_prod <- function(a, b) {
  d_op_dispatch(a, b, list(
    d_scalar_multiply_scalar, d_scalar_multiply_matrix,
    d_matrix_multiply_scalar, d_matrix_multiply_matrix
  ))
}

d_divide <- function(a, b) {
  d_op_dispatch(a, b, list(
    d_scalar_divide_scalar, d_scalar_divide_matrix,
    d_matrix_divide_scalar, d_matrix_divide_matrix
  ))
}

is_sparse <- function(x, threshold = 0.25) {
  is_sparse_matrix(x) && (nnzero(x) / length(x) < threshold)
}

d_matrix_prod <- function(a, b) {
  A <- a@x
  B <- b@x
  dA <- a@dx
  dB <- b@dx

  sum_1 <- if (is_sparse(dB)) {
    I_B <- Diagonal0(NCOL(B))
    (I_B %x% A) %*% dB
  } else {
    if (!is.matrix(A)) A <- as.matrix(A)
    if (!is.matrix(dB)) dB <- as.matrix(dB)
    IxCD(A, dB)
  }

  sum_2 <- if (is_sparse(dA)) {
    I_A <- Diagonal0(NROW(A))
    (t(B) %x% I_A) %*% dA
  } else {
    if (!is.matrix(B)) B <- as.matrix(B)
    if (!is.matrix(dA)) dA <- as.matrix(dA)
    BxID(t(B), dA)
  }

  sum_1 + sum_2
}

d_kronecker <- function(a, b) {
  A <- a@x
  B <- b@x
  dA <- a@dx
  dB <- b@dx

  m <- NROW(A)
  n <- NCOL(A)
  p <- NROW(B)
  q <- NCOL(B)
  I_n <- Diagonal0(n)
  K_qm <- commutation_matrix0(q, m)
  I_p <- Diagonal0(p)
  I_mn <- Diagonal0(m * n)
  I_pq <- Diagonal0(p * q)

  (I_n %x% K_qm %x% I_p) %*%
    ((I_mn %x% as.numeric(B)) %*% dA + (as.numeric(A) %x% I_pq) %*% dB)
  # (as.numeric(A) %x% dB + dA %x% as.numeric(B))  # double-checked for correctness
}
