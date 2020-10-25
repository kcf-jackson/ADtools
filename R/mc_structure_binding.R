#' @include class_dual_def.R
NULL

#' Combine 'dual'-class objects by Columns (dual-dual)
#' @param x A "dual" object.
#' @param y A "dual" object.
setMethod(
  "cbind2",
  signature(x = "dual", y = "dual"),
  function(x, y) {
    x@x <- cbind(x@x, y@x)
    x@dx <- rbind(x@dx, y@dx)
    x
  }
)

#' Combine 'dual'-class objects by Columns (dual-ANY)
#' @param x A "dual" object.
#' @param y ANY object.
setMethod(
  "cbind2",
  signature(x = "dual", y = "ANY"),
  function(x, y) {
    y <- empty_dual(y, ncol(x@dx))
    cbind2(x, y)
  }
)

#' Combine 'dual'-class objects by Columns (ANY-dual)
#' @param x ANY object.
#' @param y A "dual" object.
setMethod(
  "cbind2",
  signature(x = "ANY", y = "dual"),
  function(x, y) {
    x <- empty_dual(x, ncol(y@dx))
    cbind2(x, y)
  }
)


#' Combine 'dual'-class objects by Rows (dual-dual)
#' @param x A "dual" object.
#' @param y A "dual" object.
setMethod(
  "rbind2",
  signature(x = "dual", y = "dual"),
  function(x, y) {
    res_x <- rbind(x@x, y@x)
    res_dx <- rbind(x@dx, y@dx)
    r1 <- if_null_then(nrow(x@x), length(x@x))
    r2 <- if_null_then(nrow(y@x), length(y@x))
    nc <- if_null_then(ncol(x@x), length(x@x))
    g <- seq(nc)

    x_ind <- seq_along(x@x)
    res_x_ind <- x_ind + rearrange(rep((g - 1) * r2, r1), nc)
    # res_ind <- x_ind + rearrange(rep((g - 1) * r2, r1), nc)
    # res_dx[res_ind, ] <- x@dx[x_ind, ]

    y_ind <- seq_along(y@x)
    res_y_ind <- y_ind + rearrange(rep(g * r1, r2), nc)
    # res_ind <- y_ind + rearrange(rep(g * r1, r2), nc)
    # res_dx[res_ind, ] <- y@dx[y_ind, ]

    res_dx <- res_dx[order(c(res_x_ind, res_y_ind)), , drop = FALSE]

    x@x <- res_x
    x@dx <- res_dx
    x
  }
)

rearrange <- function(vec0, group_size) {
  as.numeric(t(matrix(vec0, nrow = group_size)))
}

#' Combine 'dual'-class objects by Rows (dual-ANY)
#' @param x A "dual" object.
#' @param y ANY object.
setMethod(
  "rbind2",
  signature(x = "dual", y = "ANY"),
  function(x, y) {
    y <- empty_dual(y, ncol(x@dx))
    rbind2(x, y)
  }
)

#' Combine 'dual'-class objects by Rows (ANY-dual)
#' @param x ANY object.
#' @param y A "dual" object.
setMethod(
  "rbind2",
  signature(x = "ANY", y = "dual"),
  function(x, y) {
    x <- empty_dual(x, ncol(y@dx))
    rbind2(x, y)
  }
)
