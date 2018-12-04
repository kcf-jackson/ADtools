#' @include class_dual_def.R
NULL

#' Length of an Object
#' @param x A "dual" object.
setMethod("length",
          signature(x = "dual"),
          function(x) { length(parent_of(x)) }
)
