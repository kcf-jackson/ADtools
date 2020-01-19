#' Post-multiplying a kronecker product
#' @description Compute (B \%x\% A) Z
#' @param B numeric matrix.
#' @param A numeric matrix.
#' @param Z numeric matrix.
#' @export
BxAZ <- function(B, A, Z) lboxdot(A, lcircledast(B, Z))

#' Pre-multiplying a kronecker product
#' @description Compute X (B \%x\% A)
#' @param X numeric matrix.
#' @param B numeric matrix.
#' @param A numeric matrix.
#' @export
XBxA <- function(X, B, A) rboxdot(rcircledast(X, B), A)

#' Compute (B \%x\% I) D
#' @param B numeric matrix.
#' @param D numeric matrix.
#' @export
BxID <- function(B, D) lcircledast(B, D)

#' Compute (I \%x\% C) D
#' @param C numeric matrix.
#' @param D numeric matrix.
#' @export
IxCD <- function(C, D) lboxdot(C, D)

#' Compute A (B \%x\% I)
#' @param A numeric matrix.
#' @param B numeric matrix.
#' @export
ABxI <- function(A, B) rcircledast(A, B)

#' Compute A (I \%x\% C)
#' @param A numeric matrix.
#' @param C numeric matrix.
#' @export
AIxC <- function(A, C) rboxdot(A, C)
