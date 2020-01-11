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

test_that("Test multivariate-t simulation", {
  purrr::map(2:10, function(i) {
    f <- function(sigma, df, delta) {
      set.seed(123 + i)
      rmvt0(n = 1, sigma = sigma, df = df, delta = delta)
    }
    f2 <- function(sigma, df, delta) {
      set.seed(123 + i)
      rmvt0(n = 2, sigma = sigma, df = df, delta = delta)
    }
    ctrl <- list(display = F, err_fun = rel_err, epsilon = 1e-4)
    param <- list(list(
      sigma = crossprod(randn(5, 5)),
      df = i,
      delta = runif(5, -1, 1)
    ))
    test_fs(list(f, f2), param, ctrl)
  })
})
