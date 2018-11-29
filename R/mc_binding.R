#' @include class_dual_def.R
NULL

#' Combine 'dual'-class objects by Columns
#' @param x A "dual" object.
#' @param y A "dual" object.
setMethod("cbind2",
          signature(x = "dual", y = "dual"),
          function(x, y) {
            x@x <- cbind(parent_of(x), parent_of(y))
            x@dx <- rbind(deriv_of(x), deriv_of(y))
            x
          }
)

#' Combine 'dual'-class objects by Rows.
#' @param x A "dual" object.
#' @param y A "dual" object.
setMethod("rbind2",
          signature(x = "dual", y = "dual"),
          function(x, y) {
            x@x <- rbind(parent_of(x), parent_of(y))
            x@dx <- rbind(deriv_of(x), deriv_of(y))
            x
          }
)
