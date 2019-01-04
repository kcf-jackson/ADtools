relative_diff <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  abs(x - y) / max(abs(x), abs(y))
}
