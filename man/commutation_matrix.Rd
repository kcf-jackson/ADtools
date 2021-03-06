% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optim_special_matrix.R
\name{commutation_matrix}
\alias{commutation_matrix}
\title{Commutation matrix}
\usage{
commutation_matrix(r, c)
}
\arguments{
\item{r}{integer; row dimension.}

\item{c}{integer; column dimension.}
}
\description{
Suppose A is a (r x c) matrix, then a (rc x rc) matrix K
is a commutation matrix if K %*% vec(A) = vec(t(A)).
}
\examples{
A <- randn(3, 4)
K <- commutation_matrix(3, 4)
compare <- function(x, y) {
    cbind(x, y, err = abs(x - y))
}
compare(K \%*\% vec(A), vec(t(A)))

}
