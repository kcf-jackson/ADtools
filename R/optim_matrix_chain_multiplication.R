#' Interface for optimal matrix chain multiplication
#'
#' @param ... Numeric matrices.
#' @param method "optimal" order or "natural" left-to-right order.
#'
#' @examples
#' A <- randn(20, 5)
#' B <- randn(5, 40)
#' C <- randn(40, 2)
#' system.time({ matrix_prod(A, B, C, method = "optimal") })
#' system.time({ matrix_prod(A, B, C, method = "natural") })
#'
#' @export
matrix_prod <- function(..., method = "optimal") {
  matrix_ls <- list(...)
  if (length(matrix_ls) == 2) {
    matrix_ls[[1]] %*% matrix_ls[[2]]
  } else if (method != "optimal") {
    purrr::reduce(matrix_ls, `%*%`)
  } else {
    optim_prod(matrix_ls)
  }
}

mprod <- matrix_prod   # Alias for internal use


#' Executing the matrix multiplication given the optimal order
#' @param matrix_ls A list of matrices to be multiplied
optim_prod <- function(matrix_ls) {
  d <- c(NROW(matrix_ls[[1]]), Map_dbl(NCOL, matrix_ls))
  S <- mcm_optimal_order(d)$S
  n <- length(d) - 1

  optimal_product <- function(i, j) {
    if (i == j) {
      matrix_ls[[i]]
    } else {
      k <- S[i, j]
      optimal_product(i, k) %*% optimal_product(k + 1, j)
    }
  }

  optimal_product(1, n)
}


#' Find the optimal order of multiplying a matrix chain
#' @param x A numeric vector of matrix dimensions
mcm_optimal_order <- function(x) {
  n <- length(x) - 1       # Number of matrices
  e <- new.env()
  e$S <- matrix(0, n, n)   # Keep track of the optimal order

  recursion <- function(i, j, k) {
    optim_mcm(i, k) + optim_mcm(k+1, j) + x[i] * x[k+1] * x[j+1]
  }

  optim_mcm <- function(i, j) {
    if (i == j) {
      return(0)
    } else {
      cost <- Map_dbl(recursion, i = i, j = j, k = i:(j - 1))
      e$S[i, j] <- i + which.min(cost) - 1
      min(cost)
    }
  }

  list(cost = optim_mcm(1, n), S = e$S)
}


Map_dbl <- function(...) as.numeric(Map(...))   # Type casting
