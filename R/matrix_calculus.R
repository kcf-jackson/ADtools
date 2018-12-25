# Helper functions
is_scalar <- function(x) { length(x) == 1 }


# Matrix +, -, *, /
d_op <- function(a, b, funs) {
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
  d_op(a, b, list(
    d_scalar_plus_scalar, d_scalar_plus_matrix,
    d_matrix_plus_scalar, d_matrix_plus_matrix
  ))
}

d_minus <- function(a, b) {
  d_op(a, b, list(
    d_scalar_minus_scalar, d_scalar_minus_matrix,
    d_matrix_minus_scalar, d_matrix_minus_matrix
  ))
}

d_scalar_prod <- function(a, b) {
  d_op(a, b, list(
    d_scalar_multiply_scalar, d_scalar_multiply_matrix,
    d_matrix_multiply_scalar, d_matrix_multiply_matrix
  ))
}

d_divide <- function(a, b) {
  d_op(a, b, list(
    d_scalar_divide_scalar, d_scalar_divide_matrix,
    d_matrix_divide_scalar, d_matrix_divide_matrix
  ))
}


# Other matrix functions
d_subset <- function(a, i, j) {
  x <- parent_of(a)
  dx <- deriv_of(a)
  nr <- nrow(x)
  nc <- ncol(x)
  if (missing(i) && missing(j)) {
    stop("At least one of index i and index j should be present.")  # nocov
  } else if (missing(i) && !missing(j)) {
    ind <- mapreduce(j, ~seq(nr) + nr * (.x - 1), c)
  } else if (!missing(i) && missing(j)) {
    ind <- mapreduce(i, ~.x + nr * (seq(nc) - 1), c)
  } else {
    ind <- map2reduce(i, j, ~.x + nr * (.y - 1), c)
  }
  dx[ind, , drop = FALSE]
}


d_diagonal <- function(x) {
  m0 <- parent_of(x)
  dx <- deriv_of(x)

  if (is.matrix(m0)) {
    diag_ind <- seq(1, length(m0), nrow(m0) + 1)
    return(dx[diag_ind, , drop = F])
  }

  if (is.vector(m0)) {
    new_dx <- zero_matrix0(length(m0)^2, ncol(dx))
    diag_ind <- seq(1, length(m0)^2, length(m0) + 1)
    new_dx[diag_ind, ] <- dx
    return(new_dx)
  }

  stop("The input is not a matrix or a vector.")
}


d_rowSums <- function(x) {
  px <- parent_of(x)
  px_len <- length(px)
  px_nr <- nrow(px)

  dx <- deriv_of(x)
  purrr::map(
    1:px_nr,
    ~colSums(dx[seq(.x, px_len, px_nr), ])
  ) %>%
    do.call(rbind, .)
}


d_colSums <- function(x) {
  px <- parent_of(x)
  px_len <- length(px)
  px_nr <- nrow(px)

  dx <- deriv_of(x)
  purrr::map2(
    seq(1, px_len, px_nr), seq(px_nr, px_len, px_nr),
    ~colSums(dx[.x:.y, ])
  ) %>%
    do.call(rbind, .)
}


d_solve <- function(a, inv_a) {
  da <- deriv_of(a)
  - (t(inv_a) %x% inv_a) %*% da
}


d_transpose <- function(a) {
  X <- parent_of(a)
  dX <- deriv_of(a)
  K_nq <- commutation_matrix0(nrow(X), ncol(X))
  K_nq %*% dX
}


d_XXT <- function(a) {
  X <- parent_of(a)
  dX <- deriv_of(a)

  n <- nrow(X)
  I_n <- Diagonal0(n)
  I_nn <- Diagonal0(n^2)
  K_nn <- commutation_matrix0(n, n)
  ((I_nn + K_nn) %*% (X %x% I_n)) %*% dX
}


d_chol <- function(L, a) {
  # LL^T = A
  dA <- deriv_of(a)

  n <- nrow(L)
  I_n <- Diagonal0(n)
  I_nn <- Diagonal0(n^2)
  K_nn <- commutation_matrix0(n, n)
  E_n <- elimination_matrix0(n)
  D_n <- Matrix::t(E_n)

  (D_n %*% solve(E_n %*% (I_nn + K_nn) %*% (L %x% I_n) %*% D_n, E_n)) %*% dA
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
  I_n <- Diagonal0(n)
  K_qm <- commutation_matrix0(q, m)
  I_p <- Diagonal0(p)
  I_mn <- Diagonal0(m*n)
  I_pq <- Diagonal0(p*q)

  (I_n %x% K_qm %x% I_p) %*%
    ((I_mn %x% as.numeric(B)) %*% dA + (as.numeric(A) %x% I_pq) %*% dB)
    #(as.numeric(A) %x% dB + dA %x% as.numeric(B))
}
