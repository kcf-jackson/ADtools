context("Test random variable simulation")

test_that("Test student-t simulation", {
	purrr::map(2:15, function(i) {
	  f <- function(df) {
	    set.seed(123)
	    rt0(n = 1, df = df)
	  }
	  f2 <- function(df) {
	    set.seed(123) 
	    rt0(n = i, df = df)
	  }
	  ctrl <- list(display = F, err_fun = abs_err, epsilon = 1e-6)
	  test_fs(list(f), list(list(df = 10 + i)), ctrl)
	  test_fs(list(f2), list(list(df = sample(20, i))), ctrl)
	})
})
