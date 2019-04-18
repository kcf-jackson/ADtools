
R package 'ADtools'
-------------------

Implements the forward-mode auto-differentiation for multivariate functions using the matrix-calculus notation from Magnus and Neudecker (1988). Two key features of the package are: (i) the package incorporates various optimisaton strategies to improve performance; this includes applying memoisation to cut down object construction time, using sparse matrix representation to save derivative calculation, and creating specialised matrix operations with Rcpp to reduce computation time; (ii) the package supports differentiating random variable with respect to their parameters, targetting MCMC (and in general simulation-based) applications.

[![Travis-CI Build Status](https://travis-ci.org/kcf-jackson/ADtools.svg?branch=master)](https://travis-ci.org/kcf-jackson/ADtools)

[![Coverage status](https://codecov.io/gh/kcf-jackson/ADtools/branch/master/graph/badge.svg)](https://codecov.io/github/kcf-jackson/ADtools?branch=master)

### Installation

``` r
devtools::install_github("kcf-jackson/ADtools")
```

------------------------------------------------------------------------

### Notation

Given a function ![f: X \\mapsto Y = f(X)](https://latex.codecogs.com/png.latex?f%3A%20X%20%5Cmapsto%20Y%20%3D%20f%28X%29 "f: X \mapsto Y = f(X)"), where ![X \\in R^{m \\times n}, Y \\in R^{h \\times k}](https://latex.codecogs.com/png.latex?X%20%5Cin%20R%5E%7Bm%20%5Ctimes%20n%7D%2C%20Y%20%5Cin%20R%5E%7Bh%20%5Ctimes%20k%7D "X \in R^{m \times n}, Y \in R^{h \times k}"), the Jacobina matrix of ![f](https://latex.codecogs.com/png.latex?f "f") w.r.t. ![X](https://latex.codecogs.com/png.latex?X "X") is given by

![\\dfrac{\\partial f(X)}{\\partial X}:=\\dfrac{\\partial\\,\\text{vec}\\, f(X)}{\\partial\\, (\\text{vec}X)^T} = \\dfrac{\\partial\\,\\text{vec}\\,Y}{\\partial\\,(\\text{vec}X)^T}\\in R^{mn \\times hk}.](https://latex.codecogs.com/png.latex?%5Cdfrac%7B%5Cpartial%20f%28X%29%7D%7B%5Cpartial%20X%7D%3A%3D%5Cdfrac%7B%5Cpartial%5C%2C%5Ctext%7Bvec%7D%5C%2C%20f%28X%29%7D%7B%5Cpartial%5C%2C%20%28%5Ctext%7Bvec%7DX%29%5ET%7D%20%3D%20%5Cdfrac%7B%5Cpartial%5C%2C%5Ctext%7Bvec%7D%5C%2CY%7D%7B%5Cpartial%5C%2C%28%5Ctext%7Bvec%7DX%29%5ET%7D%5Cin%20R%5E%7Bmn%20%5Ctimes%20hk%7D. "\dfrac{\partial f(X)}{\partial X}:=\dfrac{\partial\,\text{vec}\, f(X)}{\partial\, (\text{vec}X)^T} = \dfrac{\partial\,\text{vec}\,Y}{\partial\,(\text{vec}X)^T}\in R^{mn \times hk}.")

------------------------------------------------------------------------

### Example 1. Matrix multiplication

#### Function definition

Consider ![f(X, y) = X y](https://latex.codecogs.com/png.latex?f%28X%2C%20y%29%20%3D%20X%20y "f(X, y) = X y") where ![X](https://latex.codecogs.com/png.latex?X "X") is a matrix, and ![y](https://latex.codecogs.com/png.latex?y "y") is a vector.

``` r
library(ADtools)
f <- function(X, y) X %*% y
X <- randn(2, 2)
y <- matrix(c(1, 1))
print(list(X = X, y = y, f = f(X, y)))
```

    ## $X
    ##            [,1]        [,2]
    ## [1,]  0.6008076  0.28125633
    ## [2,] -0.7162415 -0.02514947
    ## 
    ## $y
    ##      [,1]
    ## [1,]    1
    ## [2,]    1
    ## 
    ## $f
    ##           [,1]
    ## [1,]  0.882064
    ## [2,] -0.741391

#### Auto-differentiation

Since ![X](https://latex.codecogs.com/png.latex?X "X") has dimension (2, 2) and ![y](https://latex.codecogs.com/png.latex?y "y") has dimension (2, 1), the input space has dimension ![2 \\times 2 + 2 \\times 1 = 6](https://latex.codecogs.com/png.latex?2%20%5Ctimes%202%20%2B%202%20%5Ctimes%201%20%3D%206 "2 \times 2 + 2 \times 1 = 6"), and the output has dimension ![2](https://latex.codecogs.com/png.latex?2 "2"), i.e. ![f](https://latex.codecogs.com/png.latex?f "f") maps ![R^6](https://latex.codecogs.com/png.latex?R%5E6 "R^6") to ![R^2](https://latex.codecogs.com/png.latex?R%5E2 "R^2") and the Jacobian of ![f](https://latex.codecogs.com/png.latex?f "f") should be ![2 \\timesb6](https://latex.codecogs.com/png.latex?2%20%5Ctimesb6 "2 \timesb6").

``` r
# Full Jacobian matrix
f_AD <- auto_diff(f, list(X = X, y = y))
get_deriv(f_AD)   # returns a Jacobian matrix
```

    ##            d_X1 d_X2 d_X3 d_X4       d_y1        d_y2
    ## d_output_1    1    0    1    0  0.6008076  0.28125633
    ## d_output_2    0    1    0    1 -0.7162415 -0.02514947

`auto_diff` also supports computing a partial Jacobian matrix. For instance, suppose we are only interested in the derivative w.r.t. `y`, then we can run

``` r
f_AD <- auto_diff(f, vary = list(y = y), fix = list(X = X))
get_deriv(f_AD)   # returns a partial Jacobian matrix
```

    ##                  d_y1        d_y2
    ## d_output_1  0.6008076  0.28125633
    ## d_output_2 -0.7162415 -0.02514947

#### Finite-differencing

It is good practice to always check the result with finite-differencing. This can be done by calling `finite_diff` which has the same interface as `auto_diff`.

``` r
f_FD <- finite_diff(f, list(X = X, y = y))
f_FD
```

    ##            d_X1 d_X2 d_X3 d_X4       d_y1        d_y2
    ## d_output_1    1    0    1    0  0.6008076  0.28125633
    ## d_output_2    0    1    0    1 -0.7162415 -0.02514947

------------------------------------------------------------------------

### Example 2. Estimating a linear regression model

#### Simulate data from ![\\quad y\_i = X\_i \\beta + \\epsilon\_i, \\quad \\epsilon\_i \\sim N(0, 1)](https://latex.codecogs.com/png.latex?%5Cquad%20y_i%20%3D%20X_i%20%5Cbeta%20%2B%20%5Cepsilon_i%2C%20%5Cquad%20%5Cepsilon_i%20%5Csim%20N%280%2C%201%29 "\quad y_i = X_i \beta + \epsilon_i, \quad \epsilon_i \sim N(0, 1)")

``` r
set.seed(123)
n <- 1000
p <- 3
X <- randn(n, p)
beta <- randn(p, 1)
y <- X %*% beta + rnorm(n)
```

#### Inference with gradient descent

``` r
gradient_descent <- function(f, vary, fix, learning_rate = 0.01, tol = 1e-6, show = F) {
  repeat {
    df <- auto_diff(f, vary, fix)
    if (show) print(df@x)
    delta <- learning_rate * as.numeric(df@dx)
    vary <- relist(unlist(vary) - delta, vary)
    if (max(abs(delta)) < tol) break
  }
  vary
}
```

``` r
lm_loss <- function(y, X, beta) sum((y - X %*% beta)^2)

# Estimate
gradient_descent(
  f = lm_loss, vary = list(beta = rnorm(p, 1)), fix = list(y = y, X = X),  learning_rate = 1e-4
) 
```

    ## $beta
    ## [1] -0.1417494 -0.3345771 -1.4484226

``` r
# Truth
t(beta)
```

    ##            [,1]       [,2]      [,3]
    ## [1,] -0.1503075 -0.3277571 -1.448165

<!-- ### Example  2b. Fitting a 2-layer Neural Network  -->
<!-- #### Simulate data  -->
<!-- ```{r} -->
<!-- logit <- function(x) exp(x) / (1 + exp(x)) -->
<!-- X <- randn(1000, 10) -->
<!-- W1 <- randn(10, 50) -->
<!-- W2 <- randn(50, 1) -->
<!-- f1 <- f2 <- logit -->
<!-- y <- f2(f1(X %*% W1) %*% W2) -->
<!-- ``` -->
<!-- #### Inference with gradient descent -->
<!-- ```{r} -->
<!-- loss_fun <- function(y, X, W1, W2, f1, f2) { -->
<!--   Z <- f1(X %*% W1) -->
<!--   yhat <- f2(Z %*% W2) -->
<!--   sum(y - yhat)^2 -->
<!-- } -->
<!-- gradient_descent( -->
<!--   loss_fun, -->
<!--   vary = list(W1 = W1, W2 = W2), -->
<!--   fix = list(y = y, X = X, f1 = logit, f2 = logit), -->
<!--   learning_rate = 1e-4,  -->
<!--   show = T -->
<!-- ) -->
<!-- ``` -->

------------------------------------------------------------------------

### Example 3. Sensitivity analysis of MCMC algorithms

#### Simulate data from ![\\quad y\_i = X\_i \\beta + \\epsilon\_i, \\quad \\epsilon\_i \\sim N(0, 1)](https://latex.codecogs.com/png.latex?%5Cquad%20y_i%20%3D%20X_i%20%5Cbeta%20%2B%20%5Cepsilon_i%2C%20%5Cquad%20%5Cepsilon_i%20%5Csim%20N%280%2C%201%29 "\quad y_i = X_i \beta + \epsilon_i, \quad \epsilon_i \sim N(0, 1)")

``` r
set.seed(123)
n <- 30  # small data
p <- 10
X <- randn(n, p)
beta <- randn(p, 1)
y <- X %*% beta + rnorm(n)
```

#### Estimating a Bayesian linear regression model

![y \\sim N(X\\beta, \\sigma^2), \\quad \\beta \\sim N(\\mathbf{b\_0}, \\mathbf{B\_0}), \\quad \\sigma^2 \\sim IG\\left(\\dfrac{\\alpha\_0}{2}, \\dfrac{\\delta\_0}{2}\\right)](https://latex.codecogs.com/png.latex?y%20%5Csim%20N%28X%5Cbeta%2C%20%5Csigma%5E2%29%2C%20%5Cquad%20%5Cbeta%20%5Csim%20N%28%5Cmathbf%7Bb_0%7D%2C%20%5Cmathbf%7BB_0%7D%29%2C%20%5Cquad%20%5Csigma%5E2%20%5Csim%20IG%5Cleft%28%5Cdfrac%7B%5Calpha_0%7D%7B2%7D%2C%20%5Cdfrac%7B%5Cdelta_0%7D%7B2%7D%5Cright%29 "y \sim N(X\beta, \sigma^2), \quad \beta \sim N(\mathbf{b_0}, \mathbf{B_0}), \quad \sigma^2 \sim IG\left(\dfrac{\alpha_0}{2}, \dfrac{\delta_0}{2}\right)")

#### Inference using Gibbs sampler

``` r
gibbs_gaussian <- function(X, y, b_0, B_0, alpha_0, delta_0, init_sigma, num_steps = 1e4) {
  if (missing(init_sigma))
    init_sigma <- 1 / sqrt(rgamma0(1, alpha_0 / 2, scale = 2 / delta_0))

  # Initialisation
  n <- length(y)
  alpha_1 <- alpha_0 + n
  sigma_g <- init_sigma
  inv_B_0 <- solve(B_0)
  inv_B_0_times_b_0 <- inv_B_0 %*% b_0
  XTX <- crossprod(X)
  XTy <- crossprod(X, y)
  beta_res <- vector("list", num_steps)
  sigma_res <- vector("list", num_steps)

  pb <- txtProgressBar(1, num_steps, style = 3)
  for (i in 1:num_steps) {
    # Update beta
    B_g <- solve(sigma_g^(-2) * XTX + inv_B_0)
    b_g <- B_g %*% (sigma_g^(-2) * XTy + inv_B_0_times_b_0)
    beta_g <- rmvnorm0(1, b_g, B_g)

    # Update sigma
    delta_g <- delta_0 + sum((y - X %*% beta_g)^2)
    sigma_g <- 1 / sqrt(rgamma0(1, alpha_1 / 2, scale = 2 / delta_g))

    # Keep track
    beta_res[[i]] <- beta_g
    sigma_res[[i]] <- sigma_g
    setTxtProgressBar(pb, i)
  }

  list(sigma = sigma_res, beta = beta_res)
}
```

#### Auto-differentiation

``` r
gibbs_deriv <- auto_diff(
  gibbs_gaussian,
  vary = list(b_0 = numeric(p), B_0 = diag(p), alpha_0 = 4, delta_0 = 4),
  fix = list(X = X, y = y, num_steps = 5000)
)
```

#### Computing the sensitivity of the posterior mean of ![b\_0](https://latex.codecogs.com/png.latex?b_0 "b_0") w.r.t. all the prior hyperparameters

``` r
library(magrittr)
library(knitr)
library(kableExtra)

matrix_ls_to_array <- function(x) {
  structure(unlist(x), dim = c(dim(x[[1]]), length(x)), dimnames = dimnames(x[[1]]))
}

tidy_mcmc <- function(mcmc_res, var0) {
  mcmc_res[[var0]] %>% 
    purrr::map(get_deriv) %>% 
    matrix_ls_to_array()
}

tidy_table <- function(x) {
  x %>% kable() %>% kable_styling() %>% scroll_box(width = "100%")
}
```

``` r
posterior_Jacobian <- apply(tidy_mcmc(gibbs_deriv, "beta"), c(1,2), mean) 
tidy_table(posterior_Jacobian)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
d\_b\_01
</th>
<th style="text-align:right;">
d\_b\_02
</th>
<th style="text-align:right;">
d\_b\_03
</th>
<th style="text-align:right;">
d\_b\_04
</th>
<th style="text-align:right;">
d\_b\_05
</th>
<th style="text-align:right;">
d\_b\_06
</th>
<th style="text-align:right;">
d\_b\_07
</th>
<th style="text-align:right;">
d\_b\_08
</th>
<th style="text-align:right;">
d\_b\_09
</th>
<th style="text-align:right;">
d\_b\_010
</th>
<th style="text-align:right;">
d\_B\_01
</th>
<th style="text-align:right;">
d\_B\_02
</th>
<th style="text-align:right;">
d\_B\_03
</th>
<th style="text-align:right;">
d\_B\_04
</th>
<th style="text-align:right;">
d\_B\_05
</th>
<th style="text-align:right;">
d\_B\_06
</th>
<th style="text-align:right;">
d\_B\_07
</th>
<th style="text-align:right;">
d\_B\_08
</th>
<th style="text-align:right;">
d\_B\_09
</th>
<th style="text-align:right;">
d\_B\_010
</th>
<th style="text-align:right;">
d\_B\_011
</th>
<th style="text-align:right;">
d\_B\_012
</th>
<th style="text-align:right;">
d\_B\_013
</th>
<th style="text-align:right;">
d\_B\_014
</th>
<th style="text-align:right;">
d\_B\_015
</th>
<th style="text-align:right;">
d\_B\_016
</th>
<th style="text-align:right;">
d\_B\_017
</th>
<th style="text-align:right;">
d\_B\_018
</th>
<th style="text-align:right;">
d\_B\_019
</th>
<th style="text-align:right;">
d\_B\_020
</th>
<th style="text-align:right;">
d\_B\_021
</th>
<th style="text-align:right;">
d\_B\_022
</th>
<th style="text-align:right;">
d\_B\_023
</th>
<th style="text-align:right;">
d\_B\_024
</th>
<th style="text-align:right;">
d\_B\_025
</th>
<th style="text-align:right;">
d\_B\_026
</th>
<th style="text-align:right;">
d\_B\_027
</th>
<th style="text-align:right;">
d\_B\_028
</th>
<th style="text-align:right;">
d\_B\_029
</th>
<th style="text-align:right;">
d\_B\_030
</th>
<th style="text-align:right;">
d\_B\_031
</th>
<th style="text-align:right;">
d\_B\_032
</th>
<th style="text-align:right;">
d\_B\_033
</th>
<th style="text-align:right;">
d\_B\_034
</th>
<th style="text-align:right;">
d\_B\_035
</th>
<th style="text-align:right;">
d\_B\_036
</th>
<th style="text-align:right;">
d\_B\_037
</th>
<th style="text-align:right;">
d\_B\_038
</th>
<th style="text-align:right;">
d\_B\_039
</th>
<th style="text-align:right;">
d\_B\_040
</th>
<th style="text-align:right;">
d\_B\_041
</th>
<th style="text-align:right;">
d\_B\_042
</th>
<th style="text-align:right;">
d\_B\_043
</th>
<th style="text-align:right;">
d\_B\_044
</th>
<th style="text-align:right;">
d\_B\_045
</th>
<th style="text-align:right;">
d\_B\_046
</th>
<th style="text-align:right;">
d\_B\_047
</th>
<th style="text-align:right;">
d\_B\_048
</th>
<th style="text-align:right;">
d\_B\_049
</th>
<th style="text-align:right;">
d\_B\_050
</th>
<th style="text-align:right;">
d\_B\_051
</th>
<th style="text-align:right;">
d\_B\_052
</th>
<th style="text-align:right;">
d\_B\_053
</th>
<th style="text-align:right;">
d\_B\_054
</th>
<th style="text-align:right;">
d\_B\_055
</th>
<th style="text-align:right;">
d\_B\_056
</th>
<th style="text-align:right;">
d\_B\_057
</th>
<th style="text-align:right;">
d\_B\_058
</th>
<th style="text-align:right;">
d\_B\_059
</th>
<th style="text-align:right;">
d\_B\_060
</th>
<th style="text-align:right;">
d\_B\_061
</th>
<th style="text-align:right;">
d\_B\_062
</th>
<th style="text-align:right;">
d\_B\_063
</th>
<th style="text-align:right;">
d\_B\_064
</th>
<th style="text-align:right;">
d\_B\_065
</th>
<th style="text-align:right;">
d\_B\_066
</th>
<th style="text-align:right;">
d\_B\_067
</th>
<th style="text-align:right;">
d\_B\_068
</th>
<th style="text-align:right;">
d\_B\_069
</th>
<th style="text-align:right;">
d\_B\_070
</th>
<th style="text-align:right;">
d\_B\_071
</th>
<th style="text-align:right;">
d\_B\_072
</th>
<th style="text-align:right;">
d\_B\_073
</th>
<th style="text-align:right;">
d\_B\_074
</th>
<th style="text-align:right;">
d\_B\_075
</th>
<th style="text-align:right;">
d\_B\_076
</th>
<th style="text-align:right;">
d\_B\_077
</th>
<th style="text-align:right;">
d\_B\_078
</th>
<th style="text-align:right;">
d\_B\_079
</th>
<th style="text-align:right;">
d\_B\_080
</th>
<th style="text-align:right;">
d\_B\_081
</th>
<th style="text-align:right;">
d\_B\_082
</th>
<th style="text-align:right;">
d\_B\_083
</th>
<th style="text-align:right;">
d\_B\_084
</th>
<th style="text-align:right;">
d\_B\_085
</th>
<th style="text-align:right;">
d\_B\_086
</th>
<th style="text-align:right;">
d\_B\_087
</th>
<th style="text-align:right;">
d\_B\_088
</th>
<th style="text-align:right;">
d\_B\_089
</th>
<th style="text-align:right;">
d\_B\_090
</th>
<th style="text-align:right;">
d\_B\_091
</th>
<th style="text-align:right;">
d\_B\_092
</th>
<th style="text-align:right;">
d\_B\_093
</th>
<th style="text-align:right;">
d\_B\_094
</th>
<th style="text-align:right;">
d\_B\_095
</th>
<th style="text-align:right;">
d\_B\_096
</th>
<th style="text-align:right;">
d\_B\_097
</th>
<th style="text-align:right;">
d\_B\_098
</th>
<th style="text-align:right;">
d\_B\_099
</th>
<th style="text-align:right;">
d\_B\_0100
</th>
<th style="text-align:right;">
d\_alpha\_01
</th>
<th style="text-align:right;">
d\_delta\_01
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
d\_output\_1
</td>
<td style="text-align:right;">
0.0563602
</td>
<td style="text-align:right;">
0.0237924
</td>
<td style="text-align:right;">
-0.0084477
</td>
<td style="text-align:right;">
-0.0114030
</td>
<td style="text-align:right;">
-0.0002252
</td>
<td style="text-align:right;">
-0.0097792
</td>
<td style="text-align:right;">
0.0069187
</td>
<td style="text-align:right;">
-0.0247872
</td>
<td style="text-align:right;">
-0.0137287
</td>
<td style="text-align:right;">
0.0183077
</td>
<td style="text-align:right;">
-0.0435960
</td>
<td style="text-align:right;">
-0.0184568
</td>
<td style="text-align:right;">
0.0065522
</td>
<td style="text-align:right;">
0.0088451
</td>
<td style="text-align:right;">
0.0001831
</td>
<td style="text-align:right;">
0.0076033
</td>
<td style="text-align:right;">
-0.0053915
</td>
<td style="text-align:right;">
0.0192429
</td>
<td style="text-align:right;">
0.0106664
</td>
<td style="text-align:right;">
-0.0142274
</td>
<td style="text-align:right;">
-0.0418466
</td>
<td style="text-align:right;">
-0.0175721
</td>
<td style="text-align:right;">
0.0062864
</td>
<td style="text-align:right;">
0.0084711
</td>
<td style="text-align:right;">
0.0001453
</td>
<td style="text-align:right;">
0.0072940
</td>
<td style="text-align:right;">
-0.0051693
</td>
<td style="text-align:right;">
0.0184709
</td>
<td style="text-align:right;">
0.0102119
</td>
<td style="text-align:right;">
-0.0136467
</td>
<td style="text-align:right;">
-0.0543180
</td>
<td style="text-align:right;">
-0.0228709
</td>
<td style="text-align:right;">
0.0082526
</td>
<td style="text-align:right;">
0.0109957
</td>
<td style="text-align:right;">
0.0002078
</td>
<td style="text-align:right;">
0.0094097
</td>
<td style="text-align:right;">
-0.0066517
</td>
<td style="text-align:right;">
0.0238892
</td>
<td style="text-align:right;">
0.0132165
</td>
<td style="text-align:right;">
-0.0176557
</td>
<td style="text-align:right;">
-0.0710547
</td>
<td style="text-align:right;">
-0.0300214
</td>
<td style="text-align:right;">
0.0106559
</td>
<td style="text-align:right;">
0.0144757
</td>
<td style="text-align:right;">
0.0003080
</td>
<td style="text-align:right;">
0.0123342
</td>
<td style="text-align:right;">
-0.0087199
</td>
<td style="text-align:right;">
0.0312328
</td>
<td style="text-align:right;">
0.0173143
</td>
<td style="text-align:right;">
-0.0230648
</td>
<td style="text-align:right;">
-0.0154343
</td>
<td style="text-align:right;">
-0.0064639
</td>
<td style="text-align:right;">
0.0023491
</td>
<td style="text-align:right;">
0.0030912
</td>
<td style="text-align:right;">
0.0001402
</td>
<td style="text-align:right;">
0.0026747
</td>
<td style="text-align:right;">
-0.0018813
</td>
<td style="text-align:right;">
0.0067737
</td>
<td style="text-align:right;">
0.0037178
</td>
<td style="text-align:right;">
-0.0049868
</td>
<td style="text-align:right;">
0.0415205
</td>
<td style="text-align:right;">
0.0174782
</td>
<td style="text-align:right;">
-0.0062463
</td>
<td style="text-align:right;">
-0.0084214
</td>
<td style="text-align:right;">
-0.0001447
</td>
<td style="text-align:right;">
-0.0071222
</td>
<td style="text-align:right;">
0.0051149
</td>
<td style="text-align:right;">
-0.0183024
</td>
<td style="text-align:right;">
-0.0101238
</td>
<td style="text-align:right;">
0.0135122
</td>
<td style="text-align:right;">
-0.1006648
</td>
<td style="text-align:right;">
-0.0424823
</td>
<td style="text-align:right;">
0.0150458
</td>
<td style="text-align:right;">
0.0203639
</td>
<td style="text-align:right;">
0.0003471
</td>
<td style="text-align:right;">
0.0174482
</td>
<td style="text-align:right;">
-0.0122811
</td>
<td style="text-align:right;">
0.0443023
</td>
<td style="text-align:right;">
0.0245283
</td>
<td style="text-align:right;">
-0.0327341
</td>
<td style="text-align:right;">
0.0098574
</td>
<td style="text-align:right;">
0.0041591
</td>
<td style="text-align:right;">
-0.0014561
</td>
<td style="text-align:right;">
-0.0019653
</td>
<td style="text-align:right;">
-0.0000690
</td>
<td style="text-align:right;">
-0.0017237
</td>
<td style="text-align:right;">
0.0011866
</td>
<td style="text-align:right;">
-0.0042596
</td>
<td style="text-align:right;">
-0.0024178
</td>
<td style="text-align:right;">
0.0032250
</td>
<td style="text-align:right;">
0.0679873
</td>
<td style="text-align:right;">
0.0286715
</td>
<td style="text-align:right;">
-0.0101940
</td>
<td style="text-align:right;">
-0.0137367
</td>
<td style="text-align:right;">
-0.0002037
</td>
<td style="text-align:right;">
-0.0117929
</td>
<td style="text-align:right;">
0.0083669
</td>
<td style="text-align:right;">
-0.0299410
</td>
<td style="text-align:right;">
-0.0164853
</td>
<td style="text-align:right;">
0.0221319
</td>
<td style="text-align:right;">
0.1140862
</td>
<td style="text-align:right;">
0.0482443
</td>
<td style="text-align:right;">
-0.0170748
</td>
<td style="text-align:right;">
-0.0231334
</td>
<td style="text-align:right;">
-0.0004124
</td>
<td style="text-align:right;">
-0.0197665
</td>
<td style="text-align:right;">
0.0139891
</td>
<td style="text-align:right;">
-0.0501711
</td>
<td style="text-align:right;">
-0.0278446
</td>
<td style="text-align:right;">
0.0371612
</td>
<td style="text-align:right;">
-0.0018903
</td>
<td style="text-align:right;">
0.0016760
</td>
</tr>
<tr>
<td style="text-align:left;">
d\_output\_2
</td>
<td style="text-align:right;">
0.0238370
</td>
<td style="text-align:right;">
0.0894589
</td>
<td style="text-align:right;">
0.0126118
</td>
<td style="text-align:right;">
-0.0154107
</td>
<td style="text-align:right;">
0.0011086
</td>
<td style="text-align:right;">
-0.0293999
</td>
<td style="text-align:right;">
0.0190536
</td>
<td style="text-align:right;">
-0.0160957
</td>
<td style="text-align:right;">
-0.0263156
</td>
<td style="text-align:right;">
0.0234575
</td>
<td style="text-align:right;">
-0.0182316
</td>
<td style="text-align:right;">
-0.0690953
</td>
<td style="text-align:right;">
-0.0097610
</td>
<td style="text-align:right;">
0.0119210
</td>
<td style="text-align:right;">
-0.0008483
</td>
<td style="text-align:right;">
0.0227682
</td>
<td style="text-align:right;">
-0.0147671
</td>
<td style="text-align:right;">
0.0125169
</td>
<td style="text-align:right;">
0.0203839
</td>
<td style="text-align:right;">
-0.0182025
</td>
<td style="text-align:right;">
-0.0176020
</td>
<td style="text-align:right;">
-0.0661112
</td>
<td style="text-align:right;">
-0.0093969
</td>
<td style="text-align:right;">
0.0114274
</td>
<td style="text-align:right;">
-0.0008823
</td>
<td style="text-align:right;">
0.0219249
</td>
<td style="text-align:right;">
-0.0142159
</td>
<td style="text-align:right;">
0.0120301
</td>
<td style="text-align:right;">
0.0195557
</td>
<td style="text-align:right;">
-0.0174906
</td>
<td style="text-align:right;">
-0.0230276
</td>
<td style="text-align:right;">
-0.0860695
</td>
<td style="text-align:right;">
-0.0118534
</td>
<td style="text-align:right;">
0.0148336
</td>
<td style="text-align:right;">
-0.0010771
</td>
<td style="text-align:right;">
0.0283107
</td>
<td style="text-align:right;">
-0.0183489
</td>
<td style="text-align:right;">
0.0154991
</td>
<td style="text-align:right;">
0.0253298
</td>
<td style="text-align:right;">
-0.0226218
</td>
<td style="text-align:right;">
-0.0301285
</td>
<td style="text-align:right;">
-0.1128720
</td>
<td style="text-align:right;">
-0.0158729
</td>
<td style="text-align:right;">
0.0197004
</td>
<td style="text-align:right;">
-0.0013563
</td>
<td style="text-align:right;">
0.0370724
</td>
<td style="text-align:right;">
-0.0240356
</td>
<td style="text-align:right;">
0.0202812
</td>
<td style="text-align:right;">
0.0331999
</td>
<td style="text-align:right;">
-0.0295620
</td>
<td style="text-align:right;">
-0.0065227
</td>
<td style="text-align:right;">
-0.0244003
</td>
<td style="text-align:right;">
-0.0033528
</td>
<td style="text-align:right;">
0.0041492
</td>
<td style="text-align:right;">
-0.0000887
</td>
<td style="text-align:right;">
0.0080566
</td>
<td style="text-align:right;">
-0.0051912
</td>
<td style="text-align:right;">
0.0043809
</td>
<td style="text-align:right;">
0.0071351
</td>
<td style="text-align:right;">
-0.0063754
</td>
<td style="text-align:right;">
0.0174932
</td>
<td style="text-align:right;">
0.0657175
</td>
<td style="text-align:right;">
0.0092822
</td>
<td style="text-align:right;">
-0.0113731
</td>
<td style="text-align:right;">
0.0008782
</td>
<td style="text-align:right;">
-0.0214383
</td>
<td style="text-align:right;">
0.0140570
</td>
<td style="text-align:right;">
-0.0118928
</td>
<td style="text-align:right;">
-0.0193728
</td>
<td style="text-align:right;">
0.0172954
</td>
<td style="text-align:right;">
-0.0425389
</td>
<td style="text-align:right;">
-0.1596827
</td>
<td style="text-align:right;">
-0.0226798
</td>
<td style="text-align:right;">
0.0274911
</td>
<td style="text-align:right;">
-0.0021323
</td>
<td style="text-align:right;">
0.0524513
</td>
<td style="text-align:right;">
-0.0337894
</td>
<td style="text-align:right;">
0.0287913
</td>
<td style="text-align:right;">
0.0469776
</td>
<td style="text-align:right;">
-0.0419310
</td>
<td style="text-align:right;">
0.0039627
</td>
<td style="text-align:right;">
0.0155594
</td>
<td style="text-align:right;">
0.0023261
</td>
<td style="text-align:right;">
-0.0025773
</td>
<td style="text-align:right;">
0.0001285
</td>
<td style="text-align:right;">
-0.0051481
</td>
<td style="text-align:right;">
0.0032252
</td>
<td style="text-align:right;">
-0.0025116
</td>
<td style="text-align:right;">
-0.0045838
</td>
<td style="text-align:right;">
0.0041002
</td>
<td style="text-align:right;">
0.0286611
</td>
<td style="text-align:right;">
0.1077736
</td>
<td style="text-align:right;">
0.0152420
</td>
<td style="text-align:right;">
-0.0185089
</td>
<td style="text-align:right;">
0.0015043
</td>
<td style="text-align:right;">
-0.0354497
</td>
<td style="text-align:right;">
0.0230232
</td>
<td style="text-align:right;">
-0.0194488
</td>
<td style="text-align:right;">
-0.0315151
</td>
<td style="text-align:right;">
0.0283396
</td>
<td style="text-align:right;">
0.0484313
</td>
<td style="text-align:right;">
0.1812958
</td>
<td style="text-align:right;">
0.0255671
</td>
<td style="text-align:right;">
-0.0313327
</td>
<td style="text-align:right;">
0.0023575
</td>
<td style="text-align:right;">
-0.0594349
</td>
<td style="text-align:right;">
0.0385257
</td>
<td style="text-align:right;">
-0.0326159
</td>
<td style="text-align:right;">
-0.0534136
</td>
<td style="text-align:right;">
0.0477799
</td>
<td style="text-align:right;">
-0.0051469
</td>
<td style="text-align:right;">
0.0045852
</td>
</tr>
<tr>
<td style="text-align:left;">
d\_output\_3
</td>
<td style="text-align:right;">
-0.0084114
</td>
<td style="text-align:right;">
0.0126292
</td>
<td style="text-align:right;">
0.0599107
</td>
<td style="text-align:right;">
0.0018085
</td>
<td style="text-align:right;">
0.0092547
</td>
<td style="text-align:right;">
-0.0037430
</td>
<td style="text-align:right;">
-0.0123475
</td>
<td style="text-align:right;">
0.0108171
</td>
<td style="text-align:right;">
-0.0013305
</td>
<td style="text-align:right;">
0.0036047
</td>
<td style="text-align:right;">
0.0066124
</td>
<td style="text-align:right;">
-0.0097464
</td>
<td style="text-align:right;">
-0.0462855
</td>
<td style="text-align:right;">
-0.0013980
</td>
<td style="text-align:right;">
-0.0071595
</td>
<td style="text-align:right;">
0.0028852
</td>
<td style="text-align:right;">
0.0095526
</td>
<td style="text-align:right;">
-0.0083519
</td>
<td style="text-align:right;">
0.0010280
</td>
<td style="text-align:right;">
-0.0028045
</td>
<td style="text-align:right;">
0.0063383
</td>
<td style="text-align:right;">
-0.0092537
</td>
<td style="text-align:right;">
-0.0446004
</td>
<td style="text-align:right;">
-0.0013569
</td>
<td style="text-align:right;">
-0.0069243
</td>
<td style="text-align:right;">
0.0027923
</td>
<td style="text-align:right;">
0.0092018
</td>
<td style="text-align:right;">
-0.0080545
</td>
<td style="text-align:right;">
0.0009649
</td>
<td style="text-align:right;">
-0.0026972
</td>
<td style="text-align:right;">
0.0080924
</td>
<td style="text-align:right;">
-0.0121098
</td>
<td style="text-align:right;">
-0.0576397
</td>
<td style="text-align:right;">
-0.0017536
</td>
<td style="text-align:right;">
-0.0089217
</td>
<td style="text-align:right;">
0.0036022
</td>
<td style="text-align:right;">
0.0119087
</td>
<td style="text-align:right;">
-0.0104340
</td>
<td style="text-align:right;">
0.0012561
</td>
<td style="text-align:right;">
-0.0034672
</td>
<td style="text-align:right;">
0.0105729
</td>
<td style="text-align:right;">
-0.0159538
</td>
<td style="text-align:right;">
-0.0755002
</td>
<td style="text-align:right;">
-0.0021952
</td>
<td style="text-align:right;">
-0.0116488
</td>
<td style="text-align:right;">
0.0047277
</td>
<td style="text-align:right;">
0.0155451
</td>
<td style="text-align:right;">
-0.0136402
</td>
<td style="text-align:right;">
0.0016906
</td>
<td style="text-align:right;">
-0.0045330
</td>
<td style="text-align:right;">
0.0022954
</td>
<td style="text-align:right;">
-0.0034302
</td>
<td style="text-align:right;">
-0.0163782
</td>
<td style="text-align:right;">
-0.0005110
</td>
<td style="text-align:right;">
-0.0024660
</td>
<td style="text-align:right;">
0.0010268
</td>
<td style="text-align:right;">
0.0033876
</td>
<td style="text-align:right;">
-0.0029563
</td>
<td style="text-align:right;">
0.0003481
</td>
<td style="text-align:right;">
-0.0009753
</td>
<td style="text-align:right;">
-0.0062532
</td>
<td style="text-align:right;">
0.0092267
</td>
<td style="text-align:right;">
0.0441780
</td>
<td style="text-align:right;">
0.0013367
</td>
<td style="text-align:right;">
0.0068547
</td>
<td style="text-align:right;">
-0.0026708
</td>
<td style="text-align:right;">
-0.0091244
</td>
<td style="text-align:right;">
0.0079848
</td>
<td style="text-align:right;">
-0.0009623
</td>
<td style="text-align:right;">
0.0026570
</td>
<td style="text-align:right;">
0.0150615
</td>
<td style="text-align:right;">
-0.0225172
</td>
<td style="text-align:right;">
-0.1070877
</td>
<td style="text-align:right;">
-0.0032507
</td>
<td style="text-align:right;">
-0.0165956
</td>
<td style="text-align:right;">
0.0066528
</td>
<td style="text-align:right;">
0.0221574
</td>
<td style="text-align:right;">
-0.0193303
</td>
<td style="text-align:right;">
0.0023574
</td>
<td style="text-align:right;">
-0.0064605
</td>
<td style="text-align:right;">
-0.0015740
</td>
<td style="text-align:right;">
0.0021662
</td>
<td style="text-align:right;">
0.0105585
</td>
<td style="text-align:right;">
0.0003631
</td>
<td style="text-align:right;">
0.0016164
</td>
<td style="text-align:right;">
-0.0006370
</td>
<td style="text-align:right;">
-0.0022252
</td>
<td style="text-align:right;">
0.0020126
</td>
<td style="text-align:right;">
-0.0002052
</td>
<td style="text-align:right;">
0.0006340
</td>
<td style="text-align:right;">
-0.0102055
</td>
<td style="text-align:right;">
0.0151798
</td>
<td style="text-align:right;">
0.0723155
</td>
<td style="text-align:right;">
0.0022143
</td>
<td style="text-align:right;">
0.0112299
</td>
<td style="text-align:right;">
-0.0045002
</td>
<td style="text-align:right;">
-0.0149087
</td>
<td style="text-align:right;">
0.0130626
</td>
<td style="text-align:right;">
-0.0015205
</td>
<td style="text-align:right;">
0.0043639
</td>
<td style="text-align:right;">
-0.0169555
</td>
<td style="text-align:right;">
0.0256346
</td>
<td style="text-align:right;">
0.1212279
</td>
<td style="text-align:right;">
0.0036235
</td>
<td style="text-align:right;">
0.0187601
</td>
<td style="text-align:right;">
-0.0075549
</td>
<td style="text-align:right;">
-0.0249942
</td>
<td style="text-align:right;">
0.0218842
</td>
<td style="text-align:right;">
-0.0027414
</td>
<td style="text-align:right;">
0.0073877
</td>
<td style="text-align:right;">
-0.0016389
</td>
<td style="text-align:right;">
0.0014498
</td>
</tr>
<tr>
<td style="text-align:left;">
d\_output\_4
</td>
<td style="text-align:right;">
-0.0113379
</td>
<td style="text-align:right;">
-0.0154346
</td>
<td style="text-align:right;">
0.0018532
</td>
<td style="text-align:right;">
0.0531843
</td>
<td style="text-align:right;">
-0.0020445
</td>
<td style="text-align:right;">
0.0030589
</td>
<td style="text-align:right;">
-0.0086514
</td>
<td style="text-align:right;">
0.0139726
</td>
<td style="text-align:right;">
0.0140063
</td>
<td style="text-align:right;">
-0.0119918
</td>
<td style="text-align:right;">
0.0088760
</td>
<td style="text-align:right;">
0.0119499
</td>
<td style="text-align:right;">
-0.0014271
</td>
<td style="text-align:right;">
-0.0410676
</td>
<td style="text-align:right;">
0.0015750
</td>
<td style="text-align:right;">
-0.0023696
</td>
<td style="text-align:right;">
0.0066918
</td>
<td style="text-align:right;">
-0.0108091
</td>
<td style="text-align:right;">
-0.0108229
</td>
<td style="text-align:right;">
0.0092525
</td>
<td style="text-align:right;">
0.0085270
</td>
<td style="text-align:right;">
0.0116430
</td>
<td style="text-align:right;">
-0.0013978
</td>
<td style="text-align:right;">
-0.0395796
</td>
<td style="text-align:right;">
0.0015007
</td>
<td style="text-align:right;">
-0.0022922
</td>
<td style="text-align:right;">
0.0064531
</td>
<td style="text-align:right;">
-0.0104267
</td>
<td style="text-align:right;">
-0.0104424
</td>
<td style="text-align:right;">
0.0089273
</td>
<td style="text-align:right;">
0.0109093
</td>
<td style="text-align:right;">
0.0149342
</td>
<td style="text-align:right;">
-0.0016888
</td>
<td style="text-align:right;">
-0.0512409
</td>
<td style="text-align:right;">
0.0019635
</td>
<td style="text-align:right;">
-0.0029595
</td>
<td style="text-align:right;">
0.0083500
</td>
<td style="text-align:right;">
-0.0134666
</td>
<td style="text-align:right;">
-0.0135144
</td>
<td style="text-align:right;">
0.0115570
</td>
<td style="text-align:right;">
0.0142717
</td>
<td style="text-align:right;">
0.0194429
</td>
<td style="text-align:right;">
-0.0023259
</td>
<td style="text-align:right;">
-0.0669827
</td>
<td style="text-align:right;">
0.0025911
</td>
<td style="text-align:right;">
-0.0038556
</td>
<td style="text-align:right;">
0.0108948
</td>
<td style="text-align:right;">
-0.0176245
</td>
<td style="text-align:right;">
-0.0176637
</td>
<td style="text-align:right;">
0.0151233
</td>
<td style="text-align:right;">
0.0030941
</td>
<td style="text-align:right;">
0.0042361
</td>
<td style="text-align:right;">
-0.0004796
</td>
<td style="text-align:right;">
-0.0145831
</td>
<td style="text-align:right;">
0.0006209
</td>
<td style="text-align:right;">
-0.0008261
</td>
<td style="text-align:right;">
0.0023743
</td>
<td style="text-align:right;">
-0.0038306
</td>
<td style="text-align:right;">
-0.0038443
</td>
<td style="text-align:right;">
0.0032840
</td>
<td style="text-align:right;">
-0.0084267
</td>
<td style="text-align:right;">
-0.0115017
</td>
<td style="text-align:right;">
0.0013645
</td>
<td style="text-align:right;">
0.0392201
</td>
<td style="text-align:right;">
-0.0014901
</td>
<td style="text-align:right;">
0.0023595
</td>
<td style="text-align:right;">
-0.0064052
</td>
<td style="text-align:right;">
0.0103305
</td>
<td style="text-align:right;">
0.0103573
</td>
<td style="text-align:right;">
-0.0088702
</td>
<td style="text-align:right;">
0.0203068
</td>
<td style="text-align:right;">
0.0276622
</td>
<td style="text-align:right;">
-0.0033596
</td>
<td style="text-align:right;">
-0.0950374
</td>
<td style="text-align:right;">
0.0036067
</td>
<td style="text-align:right;">
-0.0055143
</td>
<td style="text-align:right;">
0.0155581
</td>
<td style="text-align:right;">
-0.0249809
</td>
<td style="text-align:right;">
-0.0250567
</td>
<td style="text-align:right;">
0.0214387
</td>
<td style="text-align:right;">
-0.0020804
</td>
<td style="text-align:right;">
-0.0027580
</td>
<td style="text-align:right;">
0.0003592
</td>
<td style="text-align:right;">
0.0093672
</td>
<td style="text-align:right;">
-0.0003717
</td>
<td style="text-align:right;">
0.0005540
</td>
<td style="text-align:right;">
-0.0015695
</td>
<td style="text-align:right;">
0.0025619
</td>
<td style="text-align:right;">
0.0024725
</td>
<td style="text-align:right;">
-0.0021009
</td>
<td style="text-align:right;">
-0.0137427
</td>
<td style="text-align:right;">
-0.0187041
</td>
<td style="text-align:right;">
0.0022449
</td>
<td style="text-align:right;">
0.0642068
</td>
<td style="text-align:right;">
-0.0024173
</td>
<td style="text-align:right;">
0.0037171
</td>
<td style="text-align:right;">
-0.0104471
</td>
<td style="text-align:right;">
0.0168764
</td>
<td style="text-align:right;">
0.0169872
</td>
<td style="text-align:right;">
-0.0144770
</td>
<td style="text-align:right;">
-0.0228891
</td>
<td style="text-align:right;">
-0.0311848
</td>
<td style="text-align:right;">
0.0037625
</td>
<td style="text-align:right;">
0.1075764
</td>
<td style="text-align:right;">
-0.0041039
</td>
<td style="text-align:right;">
0.0062193
</td>
<td style="text-align:right;">
-0.0175246
</td>
<td style="text-align:right;">
0.0282676
</td>
<td style="text-align:right;">
0.0283114
</td>
<td style="text-align:right;">
-0.0241827
</td>
<td style="text-align:right;">
-0.0015053
</td>
<td style="text-align:right;">
0.0013037
</td>
</tr>
<tr>
<td style="text-align:left;">
d\_output\_5
</td>
<td style="text-align:right;">
-0.0002281
</td>
<td style="text-align:right;">
0.0011552
</td>
<td style="text-align:right;">
0.0092272
</td>
<td style="text-align:right;">
-0.0020336
</td>
<td style="text-align:right;">
0.0411938
</td>
<td style="text-align:right;">
0.0063133
</td>
<td style="text-align:right;">
-0.0123091
</td>
<td style="text-align:right;">
-0.0084113
</td>
<td style="text-align:right;">
0.0125893
</td>
<td style="text-align:right;">
0.0138534
</td>
<td style="text-align:right;">
0.0000566
</td>
<td style="text-align:right;">
-0.0008912
</td>
<td style="text-align:right;">
-0.0071011
</td>
<td style="text-align:right;">
0.0015772
</td>
<td style="text-align:right;">
-0.0318039
</td>
<td style="text-align:right;">
-0.0048977
</td>
<td style="text-align:right;">
0.0095328
</td>
<td style="text-align:right;">
0.0064830
</td>
<td style="text-align:right;">
-0.0097241
</td>
<td style="text-align:right;">
-0.0106864
</td>
<td style="text-align:right;">
0.0000808
</td>
<td style="text-align:right;">
-0.0010292
</td>
<td style="text-align:right;">
-0.0068364
</td>
<td style="text-align:right;">
0.0015514
</td>
<td style="text-align:right;">
-0.0306168
</td>
<td style="text-align:right;">
-0.0047088
</td>
<td style="text-align:right;">
0.0091676
</td>
<td style="text-align:right;">
0.0062666
</td>
<td style="text-align:right;">
-0.0093471
</td>
<td style="text-align:right;">
-0.0103233
</td>
<td style="text-align:right;">
0.0002349
</td>
<td style="text-align:right;">
-0.0011978
</td>
<td style="text-align:right;">
-0.0090269
</td>
<td style="text-align:right;">
0.0019808
</td>
<td style="text-align:right;">
-0.0396750
</td>
<td style="text-align:right;">
-0.0060507
</td>
<td style="text-align:right;">
0.0118378
</td>
<td style="text-align:right;">
0.0081236
</td>
<td style="text-align:right;">
-0.0121105
</td>
<td style="text-align:right;">
-0.0133497
</td>
<td style="text-align:right;">
0.0003193
</td>
<td style="text-align:right;">
-0.0014313
</td>
<td style="text-align:right;">
-0.0116576
</td>
<td style="text-align:right;">
0.0024500
</td>
<td style="text-align:right;">
-0.0519704
</td>
<td style="text-align:right;">
-0.0079634
</td>
<td style="text-align:right;">
0.0155278
</td>
<td style="text-align:right;">
0.0106177
</td>
<td style="text-align:right;">
-0.0158842
</td>
<td style="text-align:right;">
-0.0174866
</td>
<td style="text-align:right;">
0.0000729
</td>
<td style="text-align:right;">
-0.0003426
</td>
<td style="text-align:right;">
-0.0025778
</td>
<td style="text-align:right;">
0.0005808
</td>
<td style="text-align:right;">
-0.0114102
</td>
<td style="text-align:right;">
-0.0017454
</td>
<td style="text-align:right;">
0.0033823
</td>
<td style="text-align:right;">
0.0023051
</td>
<td style="text-align:right;">
-0.0034387
</td>
<td style="text-align:right;">
-0.0038137
</td>
<td style="text-align:right;">
-0.0001060
</td>
<td style="text-align:right;">
0.0009740
</td>
<td style="text-align:right;">
0.0068034
</td>
<td style="text-align:right;">
-0.0015136
</td>
<td style="text-align:right;">
0.0303364
</td>
<td style="text-align:right;">
0.0045396
</td>
<td style="text-align:right;">
-0.0090522
</td>
<td style="text-align:right;">
-0.0062207
</td>
<td style="text-align:right;">
0.0092580
</td>
<td style="text-align:right;">
0.0102416
</td>
<td style="text-align:right;">
0.0003605
</td>
<td style="text-align:right;">
-0.0021406
</td>
<td style="text-align:right;">
-0.0164181
</td>
<td style="text-align:right;">
0.0036653
</td>
<td style="text-align:right;">
-0.0735208
</td>
<td style="text-align:right;">
-0.0112390
</td>
<td style="text-align:right;">
0.0218753
</td>
<td style="text-align:right;">
0.0150452
</td>
<td style="text-align:right;">
-0.0224627
</td>
<td style="text-align:right;">
-0.0247656
</td>
<td style="text-align:right;">
0.0000703
</td>
<td style="text-align:right;">
0.0002555
</td>
<td style="text-align:right;">
0.0015657
</td>
<td style="text-align:right;">
-0.0004148
</td>
<td style="text-align:right;">
0.0072621
</td>
<td style="text-align:right;">
0.0011080
</td>
<td style="text-align:right;">
-0.0021197
</td>
<td style="text-align:right;">
-0.0016153
</td>
<td style="text-align:right;">
0.0021955
</td>
<td style="text-align:right;">
0.0024568
</td>
<td style="text-align:right;">
-0.0002151
</td>
<td style="text-align:right;">
0.0014708
</td>
<td style="text-align:right;">
0.0111149
</td>
<td style="text-align:right;">
-0.0024999
</td>
<td style="text-align:right;">
0.0496180
</td>
<td style="text-align:right;">
0.0075958
</td>
<td style="text-align:right;">
-0.0148551
</td>
<td style="text-align:right;">
-0.0101515
</td>
<td style="text-align:right;">
0.0150799
</td>
<td style="text-align:right;">
0.0167054
</td>
<td style="text-align:right;">
-0.0005389
</td>
<td style="text-align:right;">
0.0022545
</td>
<td style="text-align:right;">
0.0186632
</td>
<td style="text-align:right;">
-0.0040585
</td>
<td style="text-align:right;">
0.0832905
</td>
<td style="text-align:right;">
0.0127397
</td>
<td style="text-align:right;">
-0.0248779
</td>
<td style="text-align:right;">
-0.0170094
</td>
<td style="text-align:right;">
0.0255284
</td>
<td style="text-align:right;">
0.0279119
</td>
<td style="text-align:right;">
0.0021730
</td>
<td style="text-align:right;">
-0.0019360
</td>
</tr>
<tr>
<td style="text-align:left;">
d\_output\_6
</td>
<td style="text-align:right;">
-0.0098688
</td>
<td style="text-align:right;">
-0.0294112
</td>
<td style="text-align:right;">
-0.0037266
</td>
<td style="text-align:right;">
0.0030588
</td>
<td style="text-align:right;">
0.0062910
</td>
<td style="text-align:right;">
0.0545817
</td>
<td style="text-align:right;">
-0.0142370
</td>
<td style="text-align:right;">
0.0050961
</td>
<td style="text-align:right;">
0.0106808
</td>
<td style="text-align:right;">
0.0027455
</td>
<td style="text-align:right;">
0.0074037
</td>
<td style="text-align:right;">
0.0227160
</td>
<td style="text-align:right;">
0.0029085
</td>
<td style="text-align:right;">
-0.0023674
</td>
<td style="text-align:right;">
-0.0048740
</td>
<td style="text-align:right;">
-0.0421774
</td>
<td style="text-align:right;">
0.0110376
</td>
<td style="text-align:right;">
-0.0039778
</td>
<td style="text-align:right;">
-0.0082887
</td>
<td style="text-align:right;">
-0.0020596
</td>
<td style="text-align:right;">
0.0071893
</td>
<td style="text-align:right;">
0.0214976
</td>
<td style="text-align:right;">
0.0028256
</td>
<td style="text-align:right;">
-0.0022511
</td>
<td style="text-align:right;">
-0.0046574
</td>
<td style="text-align:right;">
-0.0406024
</td>
<td style="text-align:right;">
0.0106155
</td>
<td style="text-align:right;">
-0.0038148
</td>
<td style="text-align:right;">
-0.0079311
</td>
<td style="text-align:right;">
-0.0020112
</td>
<td style="text-align:right;">
0.0095534
</td>
<td style="text-align:right;">
0.0281802
</td>
<td style="text-align:right;">
0.0033070
</td>
<td style="text-align:right;">
-0.0029339
</td>
<td style="text-align:right;">
-0.0060601
</td>
<td style="text-align:right;">
-0.0525448
</td>
<td style="text-align:right;">
0.0136975
</td>
<td style="text-align:right;">
-0.0048879
</td>
<td style="text-align:right;">
-0.0102556
</td>
<td style="text-align:right;">
-0.0026391
</td>
<td style="text-align:right;">
0.0125162
</td>
<td style="text-align:right;">
0.0371822
</td>
<td style="text-align:right;">
0.0046694
</td>
<td style="text-align:right;">
-0.0041068
</td>
<td style="text-align:right;">
-0.0079687
</td>
<td style="text-align:right;">
-0.0688471
</td>
<td style="text-align:right;">
0.0179767
</td>
<td style="text-align:right;">
-0.0064027
</td>
<td style="text-align:right;">
-0.0134932
</td>
<td style="text-align:right;">
-0.0034776
</td>
<td style="text-align:right;">
0.0027213
</td>
<td style="text-align:right;">
0.0080081
</td>
<td style="text-align:right;">
0.0009340
</td>
<td style="text-align:right;">
-0.0007835
</td>
<td style="text-align:right;">
-0.0019278
</td>
<td style="text-align:right;">
-0.0150487
</td>
<td style="text-align:right;">
0.0038938
</td>
<td style="text-align:right;">
-0.0013897
</td>
<td style="text-align:right;">
-0.0028788
</td>
<td style="text-align:right;">
-0.0007959
</td>
<td style="text-align:right;">
-0.0071772
</td>
<td style="text-align:right;">
-0.0214691
</td>
<td style="text-align:right;">
-0.0027465
</td>
<td style="text-align:right;">
0.0022568
</td>
<td style="text-align:right;">
0.0045955
</td>
<td style="text-align:right;">
0.0400256
</td>
<td style="text-align:right;">
-0.0104929
</td>
<td style="text-align:right;">
0.0037511
</td>
<td style="text-align:right;">
0.0078404
</td>
<td style="text-align:right;">
0.0020365
</td>
<td style="text-align:right;">
0.0175600
</td>
<td style="text-align:right;">
0.0524005
</td>
<td style="text-align:right;">
0.0067992
</td>
<td style="text-align:right;">
-0.0054146
</td>
<td style="text-align:right;">
-0.0111059
</td>
<td style="text-align:right;">
-0.0973958
</td>
<td style="text-align:right;">
0.0251916
</td>
<td style="text-align:right;">
-0.0091098
</td>
<td style="text-align:right;">
-0.0190344
</td>
<td style="text-align:right;">
-0.0048891
</td>
<td style="text-align:right;">
-0.0015038
</td>
<td style="text-align:right;">
-0.0050386
</td>
<td style="text-align:right;">
-0.0007668
</td>
<td style="text-align:right;">
0.0004130
</td>
<td style="text-align:right;">
0.0011777
</td>
<td style="text-align:right;">
0.0095699
</td>
<td style="text-align:right;">
-0.0023890
</td>
<td style="text-align:right;">
0.0005860
</td>
<td style="text-align:right;">
0.0018435
</td>
<td style="text-align:right;">
0.0005191
</td>
<td style="text-align:right;">
-0.0117905
</td>
<td style="text-align:right;">
-0.0353101
</td>
<td style="text-align:right;">
-0.0045251
</td>
<td style="text-align:right;">
0.0036070
</td>
<td style="text-align:right;">
0.0074441
</td>
<td style="text-align:right;">
0.0657818
</td>
<td style="text-align:right;">
-0.0171937
</td>
<td style="text-align:right;">
0.0061578
</td>
<td style="text-align:right;">
0.0126600
</td>
<td style="text-align:right;">
0.0032829
</td>
<td style="text-align:right;">
-0.0201410
</td>
<td style="text-align:right;">
-0.0597365
</td>
<td style="text-align:right;">
-0.0075670
</td>
<td style="text-align:right;">
0.0063096
</td>
<td style="text-align:right;">
0.0126189
</td>
<td style="text-align:right;">
0.1103627
</td>
<td style="text-align:right;">
-0.0287740
</td>
<td style="text-align:right;">
0.0103443
</td>
<td style="text-align:right;">
0.0217428
</td>
<td style="text-align:right;">
0.0052850
</td>
<td style="text-align:right;">
0.0047726
</td>
<td style="text-align:right;">
-0.0042470
</td>
</tr>
<tr>
<td style="text-align:left;">
d\_output\_7
</td>
<td style="text-align:right;">
0.0069767
</td>
<td style="text-align:right;">
0.0190846
</td>
<td style="text-align:right;">
-0.0123265
</td>
<td style="text-align:right;">
-0.0086759
</td>
<td style="text-align:right;">
-0.0123510
</td>
<td style="text-align:right;">
-0.0142378
</td>
<td style="text-align:right;">
0.0552586
</td>
<td style="text-align:right;">
-0.0144561
</td>
<td style="text-align:right;">
-0.0051984
</td>
<td style="text-align:right;">
-0.0034054
</td>
<td style="text-align:right;">
-0.0051502
</td>
<td style="text-align:right;">
-0.0147264
</td>
<td style="text-align:right;">
0.0095087
</td>
<td style="text-align:right;">
0.0067046
</td>
<td style="text-align:right;">
0.0095651
</td>
<td style="text-align:right;">
0.0110395
</td>
<td style="text-align:right;">
-0.0427262
</td>
<td style="text-align:right;">
0.0112147
</td>
<td style="text-align:right;">
0.0040486
</td>
<td style="text-align:right;">
0.0025785
</td>
<td style="text-align:right;">
-0.0050257
</td>
<td style="text-align:right;">
-0.0137780
</td>
<td style="text-align:right;">
0.0091628
</td>
<td style="text-align:right;">
0.0064278
</td>
<td style="text-align:right;">
0.0091804
</td>
<td style="text-align:right;">
0.0106032
</td>
<td style="text-align:right;">
-0.0411484
</td>
<td style="text-align:right;">
0.0107889
</td>
<td style="text-align:right;">
0.0038468
</td>
<td style="text-align:right;">
0.0025298
</td>
<td style="text-align:right;">
-0.0067766
</td>
<td style="text-align:right;">
-0.0182367
</td>
<td style="text-align:right;">
0.0121776
</td>
<td style="text-align:right;">
0.0083447
</td>
<td style="text-align:right;">
0.0118812
</td>
<td style="text-align:right;">
0.0137223
</td>
<td style="text-align:right;">
-0.0531991
</td>
<td style="text-align:right;">
0.0139089
</td>
<td style="text-align:right;">
0.0049729
</td>
<td style="text-align:right;">
0.0032766
</td>
<td style="text-align:right;">
-0.0088785
</td>
<td style="text-align:right;">
-0.0241698
</td>
<td style="text-align:right;">
0.0155721
</td>
<td style="text-align:right;">
0.0112048
</td>
<td style="text-align:right;">
0.0155992
</td>
<td style="text-align:right;">
0.0179838
</td>
<td style="text-align:right;">
-0.0697031
</td>
<td style="text-align:right;">
0.0182154
</td>
<td style="text-align:right;">
0.0065779
</td>
<td style="text-align:right;">
0.0043023
</td>
<td style="text-align:right;">
-0.0019306
</td>
<td style="text-align:right;">
-0.0051850
</td>
<td style="text-align:right;">
0.0034839
</td>
<td style="text-align:right;">
0.0023308
</td>
<td style="text-align:right;">
0.0036023
</td>
<td style="text-align:right;">
0.0039623
</td>
<td style="text-align:right;">
-0.0152056
</td>
<td style="text-align:right;">
0.0039637
</td>
<td style="text-align:right;">
0.0013774
</td>
<td style="text-align:right;">
0.0009719
</td>
<td style="text-align:right;">
0.0050415
</td>
<td style="text-align:right;">
0.0138386
</td>
<td style="text-align:right;">
-0.0091423
</td>
<td style="text-align:right;">
-0.0064148
</td>
<td style="text-align:right;">
-0.0090883
</td>
<td style="text-align:right;">
-0.0102508
</td>
<td style="text-align:right;">
0.0408347
</td>
<td style="text-align:right;">
-0.0106819
</td>
<td style="text-align:right;">
-0.0037955
</td>
<td style="text-align:right;">
-0.0025402
</td>
<td style="text-align:right;">
-0.0123804
</td>
<td style="text-align:right;">
-0.0339138
</td>
<td style="text-align:right;">
0.0218814
</td>
<td style="text-align:right;">
0.0154400
</td>
<td style="text-align:right;">
0.0219273
</td>
<td style="text-align:right;">
0.0253097
</td>
<td style="text-align:right;">
-0.0984186
</td>
<td style="text-align:right;">
0.0258184
</td>
<td style="text-align:right;">
0.0092276
</td>
<td style="text-align:right;">
0.0060728
</td>
<td style="text-align:right;">
0.0009806
</td>
<td style="text-align:right;">
0.0032088
</td>
<td style="text-align:right;">
-0.0020595
</td>
<td style="text-align:right;">
-0.0013939
</td>
<td style="text-align:right;">
-0.0022496
</td>
<td style="text-align:right;">
-0.0024776
</td>
<td style="text-align:right;">
0.0095947
</td>
<td style="text-align:right;">
-0.0022172
</td>
<td style="text-align:right;">
-0.0008711
</td>
<td style="text-align:right;">
-0.0006325
</td>
<td style="text-align:right;">
0.0082909
</td>
<td style="text-align:right;">
0.0228280
</td>
<td style="text-align:right;">
-0.0148554
</td>
<td style="text-align:right;">
-0.0103754
</td>
<td style="text-align:right;">
-0.0147521
</td>
<td style="text-align:right;">
-0.0171160
</td>
<td style="text-align:right;">
0.0666738
</td>
<td style="text-align:right;">
-0.0174519
</td>
<td style="text-align:right;">
-0.0060219
</td>
<td style="text-align:right;">
-0.0040897
</td>
<td style="text-align:right;">
0.0143017
</td>
<td style="text-align:right;">
0.0388497
</td>
<td style="text-align:right;">
-0.0249083
</td>
<td style="text-align:right;">
-0.0176832
</td>
<td style="text-align:right;">
-0.0248713
</td>
<td style="text-align:right;">
-0.0287369
</td>
<td style="text-align:right;">
0.1117514
</td>
<td style="text-align:right;">
-0.0292839
</td>
<td style="text-align:right;">
-0.0106583
</td>
<td style="text-align:right;">
-0.0065972
</td>
<td style="text-align:right;">
-0.0050733
</td>
<td style="text-align:right;">
0.0045040
</td>
</tr>
<tr>
<td style="text-align:left;">
d\_output\_8
</td>
<td style="text-align:right;">
-0.0248146
</td>
<td style="text-align:right;">
-0.0160314
</td>
<td style="text-align:right;">
0.0108423
</td>
<td style="text-align:right;">
0.0139584
</td>
<td style="text-align:right;">
-0.0084180
</td>
<td style="text-align:right;">
0.0050645
</td>
<td style="text-align:right;">
-0.0143925
</td>
<td style="text-align:right;">
0.0664575
</td>
<td style="text-align:right;">
0.0047204
</td>
<td style="text-align:right;">
-0.0103877
</td>
<td style="text-align:right;">
0.0191708
</td>
<td style="text-align:right;">
0.0123668
</td>
<td style="text-align:right;">
-0.0083917
</td>
<td style="text-align:right;">
-0.0108354
</td>
<td style="text-align:right;">
0.0064617
</td>
<td style="text-align:right;">
-0.0038895
</td>
<td style="text-align:right;">
0.0111680
</td>
<td style="text-align:right;">
-0.0513534
</td>
<td style="text-align:right;">
-0.0036842
</td>
<td style="text-align:right;">
0.0080640
</td>
<td style="text-align:right;">
0.0184711
</td>
<td style="text-align:right;">
0.0118332
</td>
<td style="text-align:right;">
-0.0080973
</td>
<td style="text-align:right;">
-0.0104279
</td>
<td style="text-align:right;">
0.0062698
</td>
<td style="text-align:right;">
-0.0037594
</td>
<td style="text-align:right;">
0.0107460
</td>
<td style="text-align:right;">
-0.0494791
</td>
<td style="text-align:right;">
-0.0035298
</td>
<td style="text-align:right;">
0.0077601
</td>
<td style="text-align:right;">
0.0239321
</td>
<td style="text-align:right;">
0.0154337
</td>
<td style="text-align:right;">
-0.0105312
</td>
<td style="text-align:right;">
-0.0134541
</td>
<td style="text-align:right;">
0.0081341
</td>
<td style="text-align:right;">
-0.0049073
</td>
<td style="text-align:right;">
0.0138633
</td>
<td style="text-align:right;">
-0.0640396
</td>
<td style="text-align:right;">
-0.0045405
</td>
<td style="text-align:right;">
0.0100277
</td>
<td style="text-align:right;">
0.0313272
</td>
<td style="text-align:right;">
0.0202729
</td>
<td style="text-align:right;">
-0.0136765
</td>
<td style="text-align:right;">
-0.0176765
</td>
<td style="text-align:right;">
0.0106339
</td>
<td style="text-align:right;">
-0.0064237
</td>
<td style="text-align:right;">
0.0181517
</td>
<td style="text-align:right;">
-0.0838488
</td>
<td style="text-align:right;">
-0.0059581
</td>
<td style="text-align:right;">
0.0131103
</td>
<td style="text-align:right;">
0.0068257
</td>
<td style="text-align:right;">
0.0043787
</td>
<td style="text-align:right;">
-0.0030150
</td>
<td style="text-align:right;">
-0.0038234
</td>
<td style="text-align:right;">
0.0022573
</td>
<td style="text-align:right;">
-0.0013828
</td>
<td style="text-align:right;">
0.0039662
</td>
<td style="text-align:right;">
-0.0182999
</td>
<td style="text-align:right;">
-0.0012772
</td>
<td style="text-align:right;">
0.0028442
</td>
<td style="text-align:right;">
-0.0183470
</td>
<td style="text-align:right;">
-0.0118026
</td>
<td style="text-align:right;">
0.0080379
</td>
<td style="text-align:right;">
0.0103487
</td>
<td style="text-align:right;">
-0.0062357
</td>
<td style="text-align:right;">
0.0036722
</td>
<td style="text-align:right;">
-0.0106587
</td>
<td style="text-align:right;">
0.0491436
</td>
<td style="text-align:right;">
0.0034990
</td>
<td style="text-align:right;">
-0.0076974
</td>
<td style="text-align:right;">
0.0443073
</td>
<td style="text-align:right;">
0.0285990
</td>
<td style="text-align:right;">
-0.0193337
</td>
<td style="text-align:right;">
-0.0249288
</td>
<td style="text-align:right;">
0.0150590
</td>
<td style="text-align:right;">
-0.0090161
</td>
<td style="text-align:right;">
0.0256516
</td>
<td style="text-align:right;">
-0.1186513
</td>
<td style="text-align:right;">
-0.0084315
</td>
<td style="text-align:right;">
0.0185639
</td>
<td style="text-align:right;">
-0.0043253
</td>
<td style="text-align:right;">
-0.0027742
</td>
<td style="text-align:right;">
0.0018908
</td>
<td style="text-align:right;">
0.0024447
</td>
<td style="text-align:right;">
-0.0014365
</td>
<td style="text-align:right;">
0.0008622
</td>
<td style="text-align:right;">
-0.0025240
</td>
<td style="text-align:right;">
0.0115900
</td>
<td style="text-align:right;">
0.0008386
</td>
<td style="text-align:right;">
-0.0018262
</td>
<td style="text-align:right;">
-0.0299252
</td>
<td style="text-align:right;">
-0.0192940
</td>
<td style="text-align:right;">
0.0130859
</td>
<td style="text-align:right;">
0.0168353
</td>
<td style="text-align:right;">
-0.0101906
</td>
<td style="text-align:right;">
0.0060919
</td>
<td style="text-align:right;">
-0.0173783
</td>
<td style="text-align:right;">
0.0801741
</td>
<td style="text-align:right;">
0.0056394
</td>
<td style="text-align:right;">
-0.0125468
</td>
<td style="text-align:right;">
-0.0502559
</td>
<td style="text-align:right;">
-0.0325304
</td>
<td style="text-align:right;">
0.0219212
</td>
<td style="text-align:right;">
0.0282718
</td>
<td style="text-align:right;">
-0.0170770
</td>
<td style="text-align:right;">
0.0102557
</td>
<td style="text-align:right;">
-0.0291061
</td>
<td style="text-align:right;">
0.1344764
</td>
<td style="text-align:right;">
0.0095934
</td>
<td style="text-align:right;">
-0.0211085
</td>
<td style="text-align:right;">
0.0013243
</td>
<td style="text-align:right;">
-0.0011570
</td>
</tr>
<tr>
<td style="text-align:left;">
d\_output\_9
</td>
<td style="text-align:right;">
-0.0137265
</td>
<td style="text-align:right;">
-0.0263259
</td>
<td style="text-align:right;">
-0.0012825
</td>
<td style="text-align:right;">
0.0140097
</td>
<td style="text-align:right;">
0.0125907
</td>
<td style="text-align:right;">
0.0106827
</td>
<td style="text-align:right;">
-0.0051535
</td>
<td style="text-align:right;">
0.0047603
</td>
<td style="text-align:right;">
0.0523782
</td>
<td style="text-align:right;">
-0.0142371
</td>
<td style="text-align:right;">
0.0105049
</td>
<td style="text-align:right;">
0.0203420
</td>
<td style="text-align:right;">
0.0010142
</td>
<td style="text-align:right;">
-0.0108428
</td>
<td style="text-align:right;">
-0.0097490
</td>
<td style="text-align:right;">
-0.0082652
</td>
<td style="text-align:right;">
0.0040363
</td>
<td style="text-align:right;">
-0.0037271
</td>
<td style="text-align:right;">
-0.0404625
</td>
<td style="text-align:right;">
0.0110294
</td>
<td style="text-align:right;">
0.0101294
</td>
<td style="text-align:right;">
0.0193702
</td>
<td style="text-align:right;">
0.0009710
</td>
<td style="text-align:right;">
-0.0104404
</td>
<td style="text-align:right;">
-0.0094056
</td>
<td style="text-align:right;">
-0.0079379
</td>
<td style="text-align:right;">
0.0038957
</td>
<td style="text-align:right;">
-0.0035441
</td>
<td style="text-align:right;">
-0.0389316
</td>
<td style="text-align:right;">
0.0105742
</td>
<td style="text-align:right;">
0.0132482
</td>
<td style="text-align:right;">
0.0252820
</td>
<td style="text-align:right;">
0.0010811
</td>
<td style="text-align:right;">
-0.0135023
</td>
<td style="text-align:right;">
-0.0121418
</td>
<td style="text-align:right;">
-0.0102680
</td>
<td style="text-align:right;">
0.0049602
</td>
<td style="text-align:right;">
-0.0045694
</td>
<td style="text-align:right;">
-0.0504513
</td>
<td style="text-align:right;">
0.0137103
</td>
<td style="text-align:right;">
0.0173596
</td>
<td style="text-align:right;">
0.0332739
</td>
<td style="text-align:right;">
0.0016115
</td>
<td style="text-align:right;">
-0.0177979
</td>
<td style="text-align:right;">
-0.0158696
</td>
<td style="text-align:right;">
-0.0134985
</td>
<td style="text-align:right;">
0.0064877
</td>
<td style="text-align:right;">
-0.0059953
</td>
<td style="text-align:right;">
-0.0660904
</td>
<td style="text-align:right;">
0.0179709
</td>
<td style="text-align:right;">
0.0037854
</td>
<td style="text-align:right;">
0.0072126
</td>
<td style="text-align:right;">
0.0003092
</td>
<td style="text-align:right;">
-0.0038340
</td>
<td style="text-align:right;">
-0.0035696
</td>
<td style="text-align:right;">
-0.0029617
</td>
<td style="text-align:right;">
0.0013981
</td>
<td style="text-align:right;">
-0.0013066
</td>
<td style="text-align:right;">
-0.0144291
</td>
<td style="text-align:right;">
0.0038975
</td>
<td style="text-align:right;">
-0.0100939
</td>
<td style="text-align:right;">
-0.0193270
</td>
<td style="text-align:right;">
-0.0009427
</td>
<td style="text-align:right;">
0.0103660
</td>
<td style="text-align:right;">
0.0092896
</td>
<td style="text-align:right;">
0.0077513
</td>
<td style="text-align:right;">
-0.0038087
</td>
<td style="text-align:right;">
0.0035175
</td>
<td style="text-align:right;">
0.0386649
</td>
<td style="text-align:right;">
-0.0105135
</td>
<td style="text-align:right;">
0.0244790
</td>
<td style="text-align:right;">
0.0469373
</td>
<td style="text-align:right;">
0.0023616
</td>
<td style="text-align:right;">
-0.0249967
</td>
<td style="text-align:right;">
-0.0224230
</td>
<td style="text-align:right;">
-0.0190264
</td>
<td style="text-align:right;">
0.0090884
</td>
<td style="text-align:right;">
-0.0085052
</td>
<td style="text-align:right;">
-0.0934755
</td>
<td style="text-align:right;">
0.0254202
</td>
<td style="text-align:right;">
-0.0022989
</td>
<td style="text-align:right;">
-0.0045672
</td>
<td style="text-align:right;">
-0.0002893
</td>
<td style="text-align:right;">
0.0024079
</td>
<td style="text-align:right;">
0.0022603
</td>
<td style="text-align:right;">
0.0018694
</td>
<td style="text-align:right;">
-0.0008617
</td>
<td style="text-align:right;">
0.0006742
</td>
<td style="text-align:right;">
0.0091765
</td>
<td style="text-align:right;">
-0.0024870
</td>
<td style="text-align:right;">
-0.0165100
</td>
<td style="text-align:right;">
-0.0316860
</td>
<td style="text-align:right;">
-0.0015586
</td>
<td style="text-align:right;">
0.0168787
</td>
<td style="text-align:right;">
0.0151450
</td>
<td style="text-align:right;">
0.0128629
</td>
<td style="text-align:right;">
-0.0062556
</td>
<td style="text-align:right;">
0.0057496
</td>
<td style="text-align:right;">
0.0630917
</td>
<td style="text-align:right;">
-0.0171875
</td>
<td style="text-align:right;">
-0.0278707
</td>
<td style="text-align:right;">
-0.0533980
</td>
<td style="text-align:right;">
-0.0026100
</td>
<td style="text-align:right;">
0.0284008
</td>
<td style="text-align:right;">
0.0253944
</td>
<td style="text-align:right;">
0.0215901
</td>
<td style="text-align:right;">
-0.0103951
</td>
<td style="text-align:right;">
0.0096523
</td>
<td style="text-align:right;">
0.1060262
</td>
<td style="text-align:right;">
-0.0289653
</td>
<td style="text-align:right;">
0.0026358
</td>
<td style="text-align:right;">
-0.0023481
</td>
</tr>
<tr>
<td style="text-align:left;">
d\_output\_10
</td>
<td style="text-align:right;">
0.0182002
</td>
<td style="text-align:right;">
0.0233918
</td>
<td style="text-align:right;">
0.0035771
</td>
<td style="text-align:right;">
-0.0119985
</td>
<td style="text-align:right;">
0.0139130
</td>
<td style="text-align:right;">
0.0028443
</td>
<td style="text-align:right;">
-0.0034831
</td>
<td style="text-align:right;">
-0.0103165
</td>
<td style="text-align:right;">
-0.0141297
</td>
<td style="text-align:right;">
0.0596822
</td>
<td style="text-align:right;">
-0.0143240
</td>
<td style="text-align:right;">
-0.0180968
</td>
<td style="text-align:right;">
-0.0027349
</td>
<td style="text-align:right;">
0.0093091
</td>
<td style="text-align:right;">
-0.0107211
</td>
<td style="text-align:right;">
-0.0022431
</td>
<td style="text-align:right;">
0.0026900
</td>
<td style="text-align:right;">
0.0079690
</td>
<td style="text-align:right;">
0.0108864
</td>
<td style="text-align:right;">
-0.0460892
</td>
<td style="text-align:right;">
-0.0137225
</td>
<td style="text-align:right;">
-0.0177351
</td>
<td style="text-align:right;">
-0.0026194
</td>
<td style="text-align:right;">
0.0090447
</td>
<td style="text-align:right;">
-0.0102914
</td>
<td style="text-align:right;">
-0.0021748
</td>
<td style="text-align:right;">
0.0025606
</td>
<td style="text-align:right;">
0.0076897
</td>
<td style="text-align:right;">
0.0105089
</td>
<td style="text-align:right;">
-0.0444189
</td>
<td style="text-align:right;">
-0.0175059
</td>
<td style="text-align:right;">
-0.0226885
</td>
<td style="text-align:right;">
-0.0036804
</td>
<td style="text-align:right;">
0.0115968
</td>
<td style="text-align:right;">
-0.0134119
</td>
<td style="text-align:right;">
-0.0027045
</td>
<td style="text-align:right;">
0.0033309
</td>
<td style="text-align:right;">
0.0099732
</td>
<td style="text-align:right;">
0.0136509
</td>
<td style="text-align:right;">
-0.0575110
</td>
<td style="text-align:right;">
-0.0229062
</td>
<td style="text-align:right;">
-0.0294823
</td>
<td style="text-align:right;">
-0.0045493
</td>
<td style="text-align:right;">
0.0149320
</td>
<td style="text-align:right;">
-0.0175904
</td>
<td style="text-align:right;">
-0.0035512
</td>
<td style="text-align:right;">
0.0044154
</td>
<td style="text-align:right;">
0.0130328
</td>
<td style="text-align:right;">
0.0178451
</td>
<td style="text-align:right;">
-0.0752990
</td>
<td style="text-align:right;">
-0.0049863
</td>
<td style="text-align:right;">
-0.0064606
</td>
<td style="text-align:right;">
-0.0010442
</td>
<td style="text-align:right;">
0.0033310
</td>
<td style="text-align:right;">
-0.0039744
</td>
<td style="text-align:right;">
-0.0008297
</td>
<td style="text-align:right;">
0.0009354
</td>
<td style="text-align:right;">
0.0028373
</td>
<td style="text-align:right;">
0.0038813
</td>
<td style="text-align:right;">
-0.0164752
</td>
<td style="text-align:right;">
0.0135897
</td>
<td style="text-align:right;">
0.0175449
</td>
<td style="text-align:right;">
0.0026573
</td>
<td style="text-align:right;">
-0.0089232
</td>
<td style="text-align:right;">
0.0102323
</td>
<td style="text-align:right;">
0.0018901
</td>
<td style="text-align:right;">
-0.0025431
</td>
<td style="text-align:right;">
-0.0076657
</td>
<td style="text-align:right;">
-0.0105200
</td>
<td style="text-align:right;">
0.0441452
</td>
<td style="text-align:right;">
-0.0325974
</td>
<td style="text-align:right;">
-0.0419371
</td>
<td style="text-align:right;">
-0.0062793
</td>
<td style="text-align:right;">
0.0215074
</td>
<td style="text-align:right;">
-0.0247280
</td>
<td style="text-align:right;">
-0.0049924
</td>
<td style="text-align:right;">
0.0060059
</td>
<td style="text-align:right;">
0.0184541
</td>
<td style="text-align:right;">
0.0253036
</td>
<td style="text-align:right;">
-0.1065532
</td>
<td style="text-align:right;">
0.0034233
</td>
<td style="text-align:right;">
0.0042292
</td>
<td style="text-align:right;">
0.0005361
</td>
<td style="text-align:right;">
-0.0022337
</td>
<td style="text-align:right;">
0.0024777
</td>
<td style="text-align:right;">
0.0005013
</td>
<td style="text-align:right;">
-0.0005015
</td>
<td style="text-align:right;">
-0.0020698
</td>
<td style="text-align:right;">
-0.0025085
</td>
<td style="text-align:right;">
0.0105309
</td>
<td style="text-align:right;">
0.0221218
</td>
<td style="text-align:right;">
0.0284279
</td>
<td style="text-align:right;">
0.0043034
</td>
<td style="text-align:right;">
-0.0146026
</td>
<td style="text-align:right;">
0.0166983
</td>
<td style="text-align:right;">
0.0034115
</td>
<td style="text-align:right;">
-0.0042076
</td>
<td style="text-align:right;">
-0.0124814
</td>
<td style="text-align:right;">
-0.0172478
</td>
<td style="text-align:right;">
0.0721270
</td>
<td style="text-align:right;">
0.0367173
</td>
<td style="text-align:right;">
0.0472418
</td>
<td style="text-align:right;">
0.0072331
</td>
<td style="text-align:right;">
-0.0241951
</td>
<td style="text-align:right;">
0.0280998
</td>
<td style="text-align:right;">
0.0056591
</td>
<td style="text-align:right;">
-0.0070259
</td>
<td style="text-align:right;">
-0.0208813
</td>
<td style="text-align:right;">
-0.0285606
</td>
<td style="text-align:right;">
0.1206045
</td>
<td style="text-align:right;">
0.0036639
</td>
<td style="text-align:right;">
-0.0032435
</td>
</tr>
</tbody>
</table>
