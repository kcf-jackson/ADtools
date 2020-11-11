#' ADtools: Automatic Differentiation
#'
#' @importFrom stats rnorm setNames dgamma integrate pgamma rgamma runif uniroot dt pt
#' @importFrom magrittr %>%
#' @importFrom utils relist
#' @import Matrix
#' @import Rcpp
#' @importFrom RcppParallel RcppParallelLibs
#' @keywords internal
#' @useDynLib ADtools, .registration = TRUE
"_PACKAGE"

utils::globalVariables(c("."))
