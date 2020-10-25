testthat::context("Test utility functions")

testthat::test_that("Scale columns of a matrix by a vector", {
    mat0 <- matrix(1:4, 2, 2)
    vec0 <- c(1, 2)
    testthat::expect_identical(
        scale_columns_by_vector(mat0, vec0),
        matrix(c(1, 2, 6, 8), 2, 2)
    )
})

testthat::test_that("Scale rows of a matrix by a vector", {
    mat0 <- matrix(1:4, 2, 2)
    vec0 <- c(1, 2)
    testthat::expect_identical(
        scale_rows_by_vector(mat0, vec0),
        matrix(c(1, 4, 3, 8), 2, 2)
    )
})

testthat::test_that("Premultiplying by a diagonal matrix is equivalent to scaling the rows of a matrix by a vector", {
    mat0 <- matrix(1:4, 2, 2)
    vec0 <- c(1, 2)
    testthat::expect_identical(
        diag(vec0) %*% mat0,
        diag_v0_times_m0(mat0, vec0)
    )
})

# A function to turn a sparse matrix into dense matrix.
# This function is needed because a sparse matrix has
# 
as_matrix <- function(x) {
    x <- as.matrix(x)
    dimnames(x) <- NULL
    return(x)
}

testthat::test_that("Add a column vector to each column of a matrix", {
    mat0 <- matrix(1:6, nrow = 2)
    vec0 <- 1:2
    testthat::expect_identical(
        as_matrix(add_vector_to_matrix_column(vec0, mat0)),
        matrix(c(2, 4, 4, 6, 6, 8), nrow = 2)
    )
})

testthat::test_that("Add a column vector to each row of a matrix", {
    mat0 <- matrix(1:6, nrow = 2)
    vec0 <- 1:3
    testthat::expect_identical(
        as_matrix(add_vector_to_matrix_row(vec0, mat0)),
        matrix(c(2, 3, 5, 6, 8, 9), nrow = 2)
    )
})
