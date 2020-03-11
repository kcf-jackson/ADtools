#' Reshape a list of dual objects into dual object
#' @param dual_ls A list of dual objects.
#' @param tidy_x A function to join the x component of the dual object.
#' @param tidy_dx A function to join the dx component of the dual object.
# tidy_dual <- function(dual_ls, tidy_x, tidy_dx) {
#   assertthat::assert_that(is_dual(dual_ls[[1]]))
#   new(
#     "dual",
#     x = map_call(dual_ls, ~.x@x, tidy_x),
#     dx = map_call(dual_ls, ~.x@dx, tidy_dx)
#   )
# }
#


#' #' Apply a function to each element of a vector
#' #' @param .x A list or atomic vector.
#' #' @param .f A function or a formula.
#' #' @param ... Additional arguments passed on to the mapped function.
#' map_dbl <- function(.x, .f, ...) {
#'   res <- purrr::map(.x, .f, ...)
#'   if (is_dual(res[[1]])) {
#'     assertthat::assert_that(
#'       all(purrr::map_lgl(res, ~length(.x@x) == 1)),
#'       msg = "Every result must be a single double."
#'     )
#'     new(
#'       "dual",
#'       x = map_call(res, ~.x@x, c),
#'       dx = map_call(res, ~.x@dx, rbind),
#'     )
#'   } else {
#'     assertthat::assert_that(
#'       all(purrr::map_lgl(res, ~length(.x) == 1)),
#'       msg = "Every result must be a single double."
#'     )
#'     do.call(c, res)
#'   }
#' }


#' is_dual <- function(x) class(x) == "dual"
#'
#'
#' map_call <- function(x, f, g, ...) {
#'   do.call(g, purrr::map(x, f, ...))
#' }
