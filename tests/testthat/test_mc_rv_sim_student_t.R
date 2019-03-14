# context("Test random variable simulation")
#
# # Helper functions
# compare_FD_and_AD <- function(FD, AD, show = F) {
#   AD <- AD[seq(nrow(FD)), seq(ncol(FD))]
#   rel_err <- relative_diff(FD, AD)
#   if (show) {
#     print(cbind(FD = as.numeric(FD), AD = as.numeric(AD)))
#     print(glue::glue("Maximum relative error over all entries: {max(rel_err)} (entry: {which.max(rel_err)})"))
#   }
#   max(rel_err)
# }
#
# eq_transform <- function(x, y, f, show = F) {
#   f <- purrr::compose(as.numeric, f)
#   if (show) print(cbind(f(x), f(y)))
#   testthat::expect_equal(f(x), f(y))
# }
#
# aeq_transform <- function(x, y, f, show = F) {
#   f <- purrr::compose(as.numeric, f)
#   diff0 <- sum(relative_diff(f(x), f(y)))
#   if (show) print(cbind(f(x), f(y)))
#   testthat::expect_lt(diff0, 1e-8)
# }
#
# # Tests
# test_that("Test student-t simulation", {
#   df <- 3
#   # Check derivative
#   for (n in 1:10) {
#     # AD
#     df0 <- dual(df, list(A = 1, B = 4), 1)
#     set.seed(123)
#     AD_res <- rt0(n, df0)@dx
#     f <- function(param) {
#       set.seed(123)
#       rt0(n, df = param)
#     }
#     FD_res <- finite_diff_test(f, df)
#     testthat::expect_lt(compare_FD_and_AD(FD_res, AD_res), 1e-6)
#     # Check sample
#     set.seed(123)
#     AD_sample <- rt0(n, df0)@x
#     FD_sample <- f(df)
#     testthat::expect_lt(sum(abs(AD_sample - FD_sample)), 1e-8)
#   }
#   # Test cases where df has length > 1
#   param_dim <- list(A = 1, B = 2, C = 1)
#   s <- sample(10, 10)
#   df0 <- dual(s, param_dim, -1)
#   df0@dx[,1] <- 1
#
#   set.seed(123)
#   AD_res <- rt0(10, df0)
#   set.seed(123)
#   AD_res_2 <- mapreduce(1:10, ~rt0(1, df0[.x]), rbind2)
#   eq_transform(AD_res, AD_res_2, parent_of, F)
#   aeq_transform(AD_res, AD_res_2, deriv_of, F)
#
#   testthat::expect_error(rt0(9, df0))
# })
