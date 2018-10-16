commutation_matrix <- function(r, c) {
  if (missing(c)) c <- r
  entries <- expand.grid(1:r, 1:c)
  src <- r * (entries[,2] - 1) + entries[,1]
  tgt <- c * (entries[,1] - 1) + entries[,2]
  Matrix::sparseMatrix(tgt, src, x = 1)
}

elimination_matrix <- function(n) {
  entries <- expand.grid(1:n, 1:n) %>%
    cbind(src = 1:(n^2)) %>%
    dplyr::filter(Var1 >= Var2) %>%
    cbind(tgt = 1:(0.5*n*(n+1)))
  Matrix::sparseMatrix(entries[,"tgt"], entries[,"src"], x = 1)
}
