#' @include class_dual_def.R
NULL

#' Combine 'dual'-class objects by Columns
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

#' Combine 'dual'-class objects by Columns
#' @param x A "dual" object.
#' @param y ANY object.
setMethod(
  "cbind2",
  signature(x = "dual", y = "ANY"),
  function(x, y) {
    y <- dual(y, get_param_dim(x), -1)
    cbind2(x, y)
  }
)

#' Combine 'dual'-class objects by Columns
#' @param x ANY object.
#' @param y A "dual" object.
setMethod(
  "cbind2",
  signature(x = "ANY", y = "dual"),
  function(x, y) {
    x <- dual(x, get_param_dim(y), -1)
    cbind2(x, y)
  }
)

#' Combine 'dual'-class objects by Rows.
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
    res_ind <- x_ind + rearrange(rep((g - 1) * r2, r1), nc)
    res_dx[res_ind, ] <- x@dx[x_ind, ]

    y_ind <- seq_along(y@x)
    res_ind <- y_ind + rearrange(rep(g * r1, r2), nc)
    res_dx[res_ind, ] <- y@dx[y_ind, ]
    x@x <- res_x
    x@dx <- res_dx
    x
  }
)

if_null_then <- function(x, y) {
  if (is.null(x)) return(y)
  return(x)
}

rearrange <- function(vec0, group_size) {
  as.numeric(t(matrix(vec0, nrow = group_size)))
}


#' Combine 'dual'-class objects by Rows
#' @param x A "dual" object.
#' @param y ANY object.
setMethod(
  "rbind2",
  signature(x = "dual", y = "ANY"),
  function(x, y) {
    y <- dual(y, get_param_dim(x), -1)
    rbind2(x, y)
  }
)

#' Combine 'dual'-class objects by Rows
#' @param x ANY object.
#' @param y A "dual" object.
setMethod(
  "rbind2",
  signature(x = "ANY", y = "dual"),
  function(x, y) {
    x <- dual(x, get_param_dim(y), -1)
    rbind2(x, y)
  }
)
