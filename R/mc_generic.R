#' @include class_dual_def.R
NULL

setMethod("length",
          signature(x = "dual"),
          function(x) { length(parent_of(x)) }
)
