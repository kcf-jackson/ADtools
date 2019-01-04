## R package 'ADtools'

[![Travis-CI Build Status](https://travis-ci.org/<USERNAME>/<REPO>.svg?branch=master)](https://travis-ci.org/<USERNAME>/<REPO>)


## Usage
```
library(magrittr)
library(ADtools)

# Testing against finite difference
X <- randn(2, 2)
y <- randn(2, 1)
beta <- rnorm(2)

# Example 1
f2 <- function(X, y) { X %*% y }
finite_diff(f2, list(X = X, y = y))
auto_diff(f2, list(X = X, y = y)) %>% get_deriv()

# Example 2
f <- function(y, X, beta) { sum((y - X %*% beta)^2) }
finite_diff(f, list(y = y, X = X, beta = beta))
auto_diff(f, list(y = y, X = X, beta = beta)) %>% get_deriv()

# Example 2b: Suppose we are only interested in the sensitivity wrt beta
finite_diff(f, vary = list(beta = beta), fix = list(y = y, X = X))
auto_diff(f, vary = list(beta = beta), fix = list(y = y, X = X)) %>% get_deriv()
```
