#' ADtools: Automatic Differentiation
#' 
#' @importFrom stats rnorm setNames dgamma integrate pgamma rgamma runif uniroot dt pt
#' @importFrom magrittr %>%
#' @importFrom utils relist
#' @import Matrix
#' @import Rcpp
#' @keywords internal
#' @useDynLib ADtools
"_PACKAGE"

utils::globalVariables(c("."))
