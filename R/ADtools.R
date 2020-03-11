#' @keywords internal
#' @useDynLib ADtools
"_PACKAGE"

# Suppress R CMD check note
#' @importFrom stats rnorm setNames dgamma integrate pgamma rgamma runif uniroot dt pt
#' @importFrom magrittr %>%
#' @importFrom utils relist
#' @import Matrix
#' @import Rcpp
utils::globalVariables(c("."))
