
## R package ‘ADtools’

[![Travis-CI Build
Status](https://travis-ci.org/kcf-jackson/ADtools.svg?branch=master)](https://travis-ci.org/kcf-jackson/ADtools)
[![Coverage
status](https://codecov.io/gh/kcf-jackson/ADtools/branch/master/graph/badge.svg)](https://codecov.io/github/kcf-jackson/ADtools?branch=master)

Implements the forward-mode auto-differentiation for multivariate
functions using the matrix-calculus notation from Magnus and Neudecker
(1988). Two key features of the package are: (i) the package
incorporates various optimisaton strategies to improve performance; this
includes applying memoisation to cut down object construction time,
using sparse matrix representation to save derivative calculation, and
creating specialised matrix operations with Rcpp to reduce computation
time; (ii) the package supports differentiating random variable with
respect to their parameters, targetting MCMC (and in general
simulation-based) applications.

### Installation

``` r
devtools::install_github("kcf-jackson/ADtools")
```

-----

### Notation

Given a function ![f: X \\mapsto Y =
f(X)](https://latex.codecogs.com/png.latex?f%3A%20X%20%5Cmapsto%20Y%20%3D%20f%28X%29
"f: X \\mapsto Y = f(X)"), where ![X \\in R^{m \\times n}, Y \\in R^{h
\\times
k}](https://latex.codecogs.com/png.latex?X%20%5Cin%20R%5E%7Bm%20%5Ctimes%20n%7D%2C%20Y%20%5Cin%20R%5E%7Bh%20%5Ctimes%20k%7D
"X \\in R^{m \\times n}, Y \\in R^{h \\times k}"), the Jacobina matrix
of ![f](https://latex.codecogs.com/png.latex?f "f") w.r.t.
![X](https://latex.codecogs.com/png.latex?X "X") is given by   
![\\dfrac{\\partial f(X)}{\\partial
X}:=\\dfrac{\\partial\\,\\text{vec}\\, f(X)}{\\partial\\,
(\\text{vec}X)^T} =
\\dfrac{\\partial\\,\\text{vec}\\,Y}{\\partial\\,(\\text{vec}X)^T}\\in
R^{mn \\times
hk}.](https://latex.codecogs.com/png.latex?%5Cdfrac%7B%5Cpartial%20f%28X%29%7D%7B%5Cpartial%20X%7D%3A%3D%5Cdfrac%7B%5Cpartial%5C%2C%5Ctext%7Bvec%7D%5C%2C%20f%28X%29%7D%7B%5Cpartial%5C%2C%20%28%5Ctext%7Bvec%7DX%29%5ET%7D%20%3D%20%5Cdfrac%7B%5Cpartial%5C%2C%5Ctext%7Bvec%7D%5C%2CY%7D%7B%5Cpartial%5C%2C%28%5Ctext%7Bvec%7DX%29%5ET%7D%5Cin%20R%5E%7Bmn%20%5Ctimes%20hk%7D.
"\\dfrac{\\partial f(X)}{\\partial X}:=\\dfrac{\\partial\\,\\text{vec}\\, f(X)}{\\partial\\, (\\text{vec}X)^T} = \\dfrac{\\partial\\,\\text{vec}\\,Y}{\\partial\\,(\\text{vec}X)^T}\\in R^{mn \\times hk}.")  

-----

### Example 1. Matrix multiplication

#### Function definition

Consider ![f(X, y) = X
y](https://latex.codecogs.com/png.latex?f%28X%2C%20y%29%20%3D%20X%20y
"f(X, y) = X y") where ![X](https://latex.codecogs.com/png.latex?X "X")
is a matrix, and ![y](https://latex.codecogs.com/png.latex?y "y") is a
vector.

``` r
library(ADtools)
f <- function(X, y) X %*% y
X <- randn(2, 2)
y <- matrix(c(1, 1))
print(list(X = X, y = y, f = f(X, y)))
```

    ## $X
    ##            [,1]       [,2]
    ## [1,] 0.41331040  0.7085659
    ## [2,] 0.01066195 -1.2300747
    ## 
    ## $y
    ##      [,1]
    ## [1,]    1
    ## [2,]    1
    ## 
    ## $f
    ##           [,1]
    ## [1,]  1.121876
    ## [2,] -1.219413

#### Auto-differentiation

Since ![X](https://latex.codecogs.com/png.latex?X "X") has dimension (2,
2) and ![y](https://latex.codecogs.com/png.latex?y "y") has dimension
(2, 1), the input space has dimension ![2 \\times 2 + 2 \\times 1
= 6](https://latex.codecogs.com/png.latex?2%20%5Ctimes%202%20%2B%202%20%5Ctimes%201%20%3D%206
"2 \\times 2 + 2 \\times 1 = 6"), and the output has dimension
![2](https://latex.codecogs.com/png.latex?2 "2"), i.e.
![f](https://latex.codecogs.com/png.latex?f "f") maps
![R^6](https://latex.codecogs.com/png.latex?R%5E6 "R^6") to
![R^2](https://latex.codecogs.com/png.latex?R%5E2 "R^2") and the
Jacobian of ![f](https://latex.codecogs.com/png.latex?f "f") should be
![2 \\times 6
= 12](https://latex.codecogs.com/png.latex?2%20%5Ctimes%206%20%3D%2012
"2 \\times 6 = 12").

``` r
# Full Jacobian matrix
f_AD <- auto_diff(f, at = list(X = X, y = y))
f_AD@dx   # returns a Jacobian matrix
```

    ##            d_X1 d_X2 d_X3 d_X4       d_y1       d_y2
    ## d_output_1    1    0    1    0 0.41331040  0.7085659
    ## d_output_2    0    1    0    1 0.01066195 -1.2300747

`auto_diff` also supports computing a partial Jacobian matrix. For
instance, suppose we are only interested in the derivative w.r.t. `y`,
then we can run

``` r
f_AD <- auto_diff(f, at = list(X = X, y = y), wrt = "y")
f_AD@dx   # returns a partial Jacobian matrix
```

    ##                  d_y1       d_y2
    ## d_output_1 0.41331040  0.7085659
    ## d_output_2 0.01066195 -1.2300747

#### Finite-differencing

It is good practice to always check the result with finite-differencing.
This can be done by calling `finite_diff` which has the same interface
as `auto_diff`.

``` r
f_FD <- finite_diff(f, at = list(X = X, y = y))
f_FD
```

    ##            d_X1 d_X2 d_X3 d_X4       d_y1       d_y2
    ## d_output_1    1    0    1    0 0.41331039  0.7085659
    ## d_output_2    0    1    0    1 0.01066194 -1.2300747

-----

### Example 2. Estimating a linear regression model

#### Simulate data from ![\\quad y\_i = X\_i \\beta + \\epsilon\_i, \\quad \\epsilon\_i \\sim N(0, 1)](https://latex.codecogs.com/png.latex?%5Cquad%20y_i%20%3D%20X_i%20%5Cbeta%20%2B%20%5Cepsilon_i%2C%20%5Cquad%20%5Cepsilon_i%20%5Csim%20N%280%2C%201%29 "\\quad y_i = X_i \\beta + \\epsilon_i, \\quad \\epsilon_i \\sim N(0, 1)")

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
    df <- auto_diff(f, at = append(vary, fix), wrt = names(vary))
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

-----

### Example 3. Sensitivity analysis of MCMC algorithms

#### Simulate data from ![\\quad y\_i = X\_i \\beta + \\epsilon\_i, \\quad \\epsilon\_i \\sim N(0, 1)](https://latex.codecogs.com/png.latex?%5Cquad%20y_i%20%3D%20X_i%20%5Cbeta%20%2B%20%5Cepsilon_i%2C%20%5Cquad%20%5Cepsilon_i%20%5Csim%20N%280%2C%201%29 "\\quad y_i = X_i \\beta + \\epsilon_i, \\quad \\epsilon_i \\sim N(0, 1)")

``` r
set.seed(123)
n <- 30  # small data
p <- 10
X <- randn(n, p)
beta <- randn(p, 1)
y <- X %*% beta + rnorm(n)
```

#### Estimating a Bayesian linear regression model

  
![y \\sim N(X\\beta, \\sigma^2), \\quad \\beta \\sim N(\\mathbf{b\_0},
\\mathbf{B\_0}), \\quad \\sigma^2 \\sim IG\\left(\\dfrac{\\alpha\_0}{2},
\\dfrac{\\delta\_0}{2}\\right)](https://latex.codecogs.com/png.latex?y%20%5Csim%20N%28X%5Cbeta%2C%20%5Csigma%5E2%29%2C%20%5Cquad%20%5Cbeta%20%5Csim%20N%28%5Cmathbf%7Bb_0%7D%2C%20%5Cmathbf%7BB_0%7D%29%2C%20%5Cquad%20%5Csigma%5E2%20%5Csim%20IG%5Cleft%28%5Cdfrac%7B%5Calpha_0%7D%7B2%7D%2C%20%5Cdfrac%7B%5Cdelta_0%7D%7B2%7D%5Cright%29
"y \\sim N(X\\beta, \\sigma^2), \\quad \\beta \\sim N(\\mathbf{b_0}, \\mathbf{B_0}), \\quad \\sigma^2 \\sim IG\\left(\\dfrac{\\alpha_0}{2}, \\dfrac{\\delta_0}{2}\\right)")  

#### Inference using Gibbs sampler

``` r
gibbs_gaussian <- function(X, y, b_0, B_0, alpha_0, delta_0, num_steps = 1e4) {
  # Initialisation
  init_sigma <- 1 / sqrt(rgamma0(1, alpha_0 / 2, scale = 2 / delta_0))
  
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
    beta_g <- t(rmvnorm0(1, b_g, B_g))

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
  at = list(
    b_0 = numeric(p), B_0 = diag(p), alpha_0 = 4, delta_0 = 4,
    X = X, y = y, num_steps = 5000
  ),
  wrt = c("b_0", "B_0", "alpha_0", "delta_0")
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
    purrr::map(~.x@dx) %>% 
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

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; ">

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

0.0569887

</td>

<td style="text-align:right;">

0.0240169

</td>

<td style="text-align:right;">

\-0.0085018

</td>

<td style="text-align:right;">

\-0.0114905

</td>

<td style="text-align:right;">

\-0.0002056

</td>

<td style="text-align:right;">

\-0.0099043

</td>

<td style="text-align:right;">

0.0069673

</td>

<td style="text-align:right;">

\-0.0250403

</td>

<td style="text-align:right;">

\-0.0138359

</td>

<td style="text-align:right;">

0.0184303

</td>

<td style="text-align:right;">

\-0.0441342

</td>

<td style="text-align:right;">

\-0.0186824

</td>

<td style="text-align:right;">

0.0065791

</td>

<td style="text-align:right;">

0.0089195

</td>

<td style="text-align:right;">

0.0001676

</td>

<td style="text-align:right;">

0.0077123

</td>

<td style="text-align:right;">

\-0.0054299

</td>

<td style="text-align:right;">

0.0194540

</td>

<td style="text-align:right;">

0.0107661

</td>

<td style="text-align:right;">

\-0.0143304

</td>

<td style="text-align:right;">

\-0.0422326

</td>

<td style="text-align:right;">

\-0.0177185

</td>

<td style="text-align:right;">

0.0062827

</td>

<td style="text-align:right;">

0.0085277

</td>

<td style="text-align:right;">

0.0001388

</td>

<td style="text-align:right;">

0.0073738

</td>

<td style="text-align:right;">

\-0.0051957

</td>

<td style="text-align:right;">

0.0186088

</td>

<td style="text-align:right;">

0.0102891

</td>

<td style="text-align:right;">

\-0.0137093

</td>

<td style="text-align:right;">

\-0.0548737

</td>

<td style="text-align:right;">

\-0.0230575

</td>

<td style="text-align:right;">

0.0082968

</td>

<td style="text-align:right;">

0.0110663

</td>

<td style="text-align:right;">

0.0001859

</td>

<td style="text-align:right;">

0.0095235

</td>

<td style="text-align:right;">

\-0.0066902

</td>

<td style="text-align:right;">

0.0241129

</td>

<td style="text-align:right;">

0.0132928

</td>

<td style="text-align:right;">

\-0.0177620

</td>

<td style="text-align:right;">

\-0.0718000

</td>

<td style="text-align:right;">

\-0.0302760

</td>

<td style="text-align:right;">

0.0107191

</td>

<td style="text-align:right;">

0.0145789

</td>

<td style="text-align:right;">

0.0002573

</td>

<td style="text-align:right;">

0.0124716

</td>

<td style="text-align:right;">

\-0.0087673

</td>

<td style="text-align:right;">

0.0315428

</td>

<td style="text-align:right;">

0.0174258

</td>

<td style="text-align:right;">

\-0.0232140

</td>

<td style="text-align:right;">

\-0.0156627

</td>

<td style="text-align:right;">

\-0.0065633

</td>

<td style="text-align:right;">

0.0023708

</td>

<td style="text-align:right;">

0.0031359

</td>

<td style="text-align:right;">

0.0001308

</td>

<td style="text-align:right;">

0.0027179

</td>

<td style="text-align:right;">

\-0.0019058

</td>

<td style="text-align:right;">

0.0068750

</td>

<td style="text-align:right;">

0.0037912

</td>

<td style="text-align:right;">

\-0.0050542

</td>

<td style="text-align:right;">

0.0418913

</td>

<td style="text-align:right;">

0.0176076

</td>

<td style="text-align:right;">

\-0.0062615

</td>

<td style="text-align:right;">

\-0.0084694

</td>

<td style="text-align:right;">

\-0.0001351

</td>

<td style="text-align:right;">

\-0.0071959

</td>

<td style="text-align:right;">

0.0051419

</td>

<td style="text-align:right;">

\-0.0184411

</td>

<td style="text-align:right;">

\-0.0101835

</td>

<td style="text-align:right;">

0.0135673

</td>

<td style="text-align:right;">

\-0.1016804

</td>

<td style="text-align:right;">

\-0.0428326

</td>

<td style="text-align:right;">

0.0151174

</td>

<td style="text-align:right;">

0.0204964

</td>

<td style="text-align:right;">

0.0003181

</td>

<td style="text-align:right;">

0.0176513

</td>

<td style="text-align:right;">

\-0.0123461

</td>

<td style="text-align:right;">

0.0446991

</td>

<td style="text-align:right;">

0.0246896

</td>

<td style="text-align:right;">

\-0.0329021

</td>

<td style="text-align:right;">

0.0099753

</td>

<td style="text-align:right;">

0.0042146

</td>

<td style="text-align:right;">

\-0.0014552

</td>

<td style="text-align:right;">

\-0.0019818

</td>

<td style="text-align:right;">

\-0.0000670

</td>

<td style="text-align:right;">

\-0.0017475

</td>

<td style="text-align:right;">

0.0011945

</td>

<td style="text-align:right;">

\-0.0043029

</td>

<td style="text-align:right;">

\-0.0024432

</td>

<td style="text-align:right;">

0.0032439

</td>

<td style="text-align:right;">

0.0687012

</td>

<td style="text-align:right;">

0.0289268

</td>

<td style="text-align:right;">

\-0.0102440

</td>

<td style="text-align:right;">

\-0.0138316

</td>

<td style="text-align:right;">

\-0.0001937

</td>

<td style="text-align:right;">

\-0.0119421

</td>

<td style="text-align:right;">

0.0084227

</td>

<td style="text-align:right;">

\-0.0302179

</td>

<td style="text-align:right;">

\-0.0166070

</td>

<td style="text-align:right;">

0.0222462

</td>

<td style="text-align:right;">

0.1152496

</td>

<td style="text-align:right;">

0.0486374

</td>

<td style="text-align:right;">

\-0.0171780

</td>

<td style="text-align:right;">

\-0.0232862

</td>

<td style="text-align:right;">

\-0.0003697

</td>

<td style="text-align:right;">

\-0.0199968

</td>

<td style="text-align:right;">

0.0140737

</td>

<td style="text-align:right;">

\-0.0506390

</td>

<td style="text-align:right;">

\-0.0280183

</td>

<td style="text-align:right;">

0.0373722

</td>

<td style="text-align:right;">

\-0.0018965

</td>

<td style="text-align:right;">

0.0016811

</td>

</tr>

<tr>

<td style="text-align:left;">

d\_output\_2

</td>

<td style="text-align:right;">

0.0240377

</td>

<td style="text-align:right;">

0.0904281

</td>

<td style="text-align:right;">

0.0127987

</td>

<td style="text-align:right;">

\-0.0154678

</td>

<td style="text-align:right;">

0.0011496

</td>

<td style="text-align:right;">

\-0.0297209

</td>

<td style="text-align:right;">

0.0192483

</td>

<td style="text-align:right;">

\-0.0162636

</td>

<td style="text-align:right;">

\-0.0265724

</td>

<td style="text-align:right;">

0.0235886

</td>

<td style="text-align:right;">

\-0.0183982

</td>

<td style="text-align:right;">

\-0.0700067

</td>

<td style="text-align:right;">

\-0.0099287

</td>

<td style="text-align:right;">

0.0119620

</td>

<td style="text-align:right;">

\-0.0008713

</td>

<td style="text-align:right;">

0.0230590

</td>

<td style="text-align:right;">

\-0.0149584

</td>

<td style="text-align:right;">

0.0126714

</td>

<td style="text-align:right;">

0.0206031

</td>

<td style="text-align:right;">

\-0.0183164

</td>

<td style="text-align:right;">

\-0.0177858

</td>

<td style="text-align:right;">

\-0.0670336

</td>

<td style="text-align:right;">

\-0.0095958

</td>

<td style="text-align:right;">

0.0114911

</td>

<td style="text-align:right;">

\-0.0009177

</td>

<td style="text-align:right;">

0.0222430

</td>

<td style="text-align:right;">

\-0.0144052

</td>

<td style="text-align:right;">

0.0121760

</td>

<td style="text-align:right;">

0.0198080

</td>

<td style="text-align:right;">

\-0.0176238

</td>

<td style="text-align:right;">

\-0.0232262

</td>

<td style="text-align:right;">

\-0.0869873

</td>

<td style="text-align:right;">

\-0.0120353

</td>

<td style="text-align:right;">

0.0149244

</td>

<td style="text-align:right;">

\-0.0011250

</td>

<td style="text-align:right;">

0.0286433

</td>

<td style="text-align:right;">

\-0.0185164

</td>

<td style="text-align:right;">

0.0156581

</td>

<td style="text-align:right;">

0.0255766

</td>

<td style="text-align:right;">

\-0.0227551

</td>

<td style="text-align:right;">

\-0.0303404

</td>

<td style="text-align:right;">

\-0.1139404

</td>

<td style="text-align:right;">

\-0.0160980

</td>

<td style="text-align:right;">

0.0197683

</td>

<td style="text-align:right;">

\-0.0014266

</td>

<td style="text-align:right;">

0.0374066

</td>

<td style="text-align:right;">

\-0.0242106

</td>

<td style="text-align:right;">

0.0204458

</td>

<td style="text-align:right;">

0.0334838

</td>

<td style="text-align:right;">

\-0.0296973

</td>

<td style="text-align:right;">

\-0.0066046

</td>

<td style="text-align:right;">

\-0.0247385

</td>

<td style="text-align:right;">

\-0.0034116

</td>

<td style="text-align:right;">

0.0041890

</td>

<td style="text-align:right;">

\-0.0001068

</td>

<td style="text-align:right;">

0.0081420

</td>

<td style="text-align:right;">

\-0.0052474

</td>

<td style="text-align:right;">

0.0044657

</td>

<td style="text-align:right;">

0.0072442

</td>

<td style="text-align:right;">

\-0.0064331

</td>

<td style="text-align:right;">

0.0176188

</td>

<td style="text-align:right;">

0.0663741

</td>

<td style="text-align:right;">

0.0094205

</td>

<td style="text-align:right;">

\-0.0114027

</td>

<td style="text-align:right;">

0.0009022

</td>

<td style="text-align:right;">

\-0.0216447

</td>

<td style="text-align:right;">

0.0141864

</td>

<td style="text-align:right;">

\-0.0120085

</td>

<td style="text-align:right;">

\-0.0195519

</td>

<td style="text-align:right;">

0.0173548

</td>

<td style="text-align:right;">

\-0.0428613

</td>

<td style="text-align:right;">

\-0.1613073

</td>

<td style="text-align:right;">

\-0.0230086

</td>

<td style="text-align:right;">

0.0275584

</td>

<td style="text-align:right;">

\-0.0021971

</td>

<td style="text-align:right;">

0.0529932

</td>

<td style="text-align:right;">

\-0.0341142

</td>

<td style="text-align:right;">

0.0290588

</td>

<td style="text-align:right;">

0.0474033

</td>

<td style="text-align:right;">

\-0.0421258

</td>

<td style="text-align:right;">

0.0040008

</td>

<td style="text-align:right;">

0.0157923

</td>

<td style="text-align:right;">

0.0023707

</td>

<td style="text-align:right;">

\-0.0025788

</td>

<td style="text-align:right;">

0.0001263

</td>

<td style="text-align:right;">

\-0.0052141

</td>

<td style="text-align:right;">

0.0032759

</td>

<td style="text-align:right;">

\-0.0025446

</td>

<td style="text-align:right;">

\-0.0046489

</td>

<td style="text-align:right;">

0.0041278

</td>

<td style="text-align:right;">

0.0289120

</td>

<td style="text-align:right;">

0.1089914

</td>

<td style="text-align:right;">

0.0154794

</td>

<td style="text-align:right;">

\-0.0185779

</td>

<td style="text-align:right;">

0.0015487

</td>

<td style="text-align:right;">

\-0.0358662

</td>

<td style="text-align:right;">

0.0232800

</td>

<td style="text-align:right;">

\-0.0196546

</td>

<td style="text-align:right;">

\-0.0318409

</td>

<td style="text-align:right;">

0.0285183

</td>

<td style="text-align:right;">

0.0487821

</td>

<td style="text-align:right;">

0.1830006

</td>

<td style="text-align:right;">

0.0259046

</td>

<td style="text-align:right;">

\-0.0314158

</td>

<td style="text-align:right;">

0.0024403

</td>

<td style="text-align:right;">

\-0.0600001

</td>

<td style="text-align:right;">

0.0388629

</td>

<td style="text-align:right;">

\-0.0329070

</td>

<td style="text-align:right;">

\-0.0538610

</td>

<td style="text-align:right;">

0.0480035

</td>

<td style="text-align:right;">

\-0.0051331

</td>

<td style="text-align:right;">

0.0045483

</td>

</tr>

<tr>

<td style="text-align:left;">

d\_output\_3

</td>

<td style="text-align:right;">

\-0.0085091

</td>

<td style="text-align:right;">

0.0128582

</td>

<td style="text-align:right;">

0.0606947

</td>

<td style="text-align:right;">

0.0018708

</td>

<td style="text-align:right;">

0.0093748

</td>

<td style="text-align:right;">

\-0.0038372

</td>

<td style="text-align:right;">

\-0.0124373

</td>

<td style="text-align:right;">

0.0108905

</td>

<td style="text-align:right;">

\-0.0013585

</td>

<td style="text-align:right;">

0.0035816

</td>

<td style="text-align:right;">

0.0066996

</td>

<td style="text-align:right;">

\-0.0100238

</td>

<td style="text-align:right;">

\-0.0470345

</td>

<td style="text-align:right;">

\-0.0014486

</td>

<td style="text-align:right;">

\-0.0072611

</td>

<td style="text-align:right;">

0.0030017

</td>

<td style="text-align:right;">

0.0096190

</td>

<td style="text-align:right;">

\-0.0084236

</td>

<td style="text-align:right;">

0.0010817

</td>

<td style="text-align:right;">

\-0.0028071

</td>

<td style="text-align:right;">

0.0064625

</td>

<td style="text-align:right;">

\-0.0094615

</td>

<td style="text-align:right;">

\-0.0454735

</td>

<td style="text-align:right;">

\-0.0014125

</td>

<td style="text-align:right;">

\-0.0070568

</td>

<td style="text-align:right;">

0.0028694

</td>

<td style="text-align:right;">

0.0093231

</td>

<td style="text-align:right;">

\-0.0081605

</td>

<td style="text-align:right;">

0.0009901

</td>

<td style="text-align:right;">

\-0.0026940

</td>

<td style="text-align:right;">

0.0081889

</td>

<td style="text-align:right;">

\-0.0123584

</td>

<td style="text-align:right;">

\-0.0585021

</td>

<td style="text-align:right;">

\-0.0017899

</td>

<td style="text-align:right;">

\-0.0090490

</td>

<td style="text-align:right;">

0.0037064

</td>

<td style="text-align:right;">

0.0120127

</td>

<td style="text-align:right;">

\-0.0105169

</td>

<td style="text-align:right;">

0.0013044

</td>

<td style="text-align:right;">

\-0.0034779

</td>

<td style="text-align:right;">

0.0106775

</td>

<td style="text-align:right;">

\-0.0162045

</td>

<td style="text-align:right;">

\-0.0763825

</td>

<td style="text-align:right;">

\-0.0022640

</td>

<td style="text-align:right;">

\-0.0117940

</td>

<td style="text-align:right;">

0.0048318

</td>

<td style="text-align:right;">

0.0156772

</td>

<td style="text-align:right;">

\-0.0137081

</td>

<td style="text-align:right;">

0.0017181

</td>

<td style="text-align:right;">

\-0.0045021

</td>

<td style="text-align:right;">

0.0023308

</td>

<td style="text-align:right;">

\-0.0035036

</td>

<td style="text-align:right;">

\-0.0166668

</td>

<td style="text-align:right;">

\-0.0005309

</td>

<td style="text-align:right;">

\-0.0025077

</td>

<td style="text-align:right;">

0.0010566

</td>

<td style="text-align:right;">

0.0034371

</td>

<td style="text-align:right;">

\-0.0029829

</td>

<td style="text-align:right;">

0.0003518

</td>

<td style="text-align:right;">

\-0.0009703

</td>

<td style="text-align:right;">

\-0.0063362

</td>

<td style="text-align:right;">

0.0093897

</td>

<td style="text-align:right;">

0.0447738

</td>

<td style="text-align:right;">

0.0013822

</td>

<td style="text-align:right;">

0.0069415

</td>

<td style="text-align:right;">

\-0.0027284

</td>

<td style="text-align:right;">

\-0.0091884

</td>

<td style="text-align:right;">

0.0080518

</td>

<td style="text-align:right;">

\-0.0009880

</td>

<td style="text-align:right;">

0.0026380

</td>

<td style="text-align:right;">

0.0152430

</td>

<td style="text-align:right;">

\-0.0228882

</td>

<td style="text-align:right;">

\-0.1084255

</td>

<td style="text-align:right;">

\-0.0033729

</td>

<td style="text-align:right;">

\-0.0168030

</td>

<td style="text-align:right;">

0.0068048

</td>

<td style="text-align:right;">

0.0223176

</td>

<td style="text-align:right;">

\-0.0194609

</td>

<td style="text-align:right;">

0.0023904

</td>

<td style="text-align:right;">

\-0.0063921

</td>

<td style="text-align:right;">

\-0.0015991

</td>

<td style="text-align:right;">

0.0022439

</td>

<td style="text-align:right;">

0.0107478

</td>

<td style="text-align:right;">

0.0003773

</td>

<td style="text-align:right;">

0.0016300

</td>

<td style="text-align:right;">

\-0.0006714

</td>

<td style="text-align:right;">

\-0.0022394

</td>

<td style="text-align:right;">

0.0020392

</td>

<td style="text-align:right;">

\-0.0002404

</td>

<td style="text-align:right;">

0.0006235

</td>

<td style="text-align:right;">

\-0.0103342

</td>

<td style="text-align:right;">

0.0154750

</td>

<td style="text-align:right;">

0.0733211

</td>

<td style="text-align:right;">

0.0022912

</td>

<td style="text-align:right;">

0.0113836

</td>

<td style="text-align:right;">

\-0.0046196

</td>

<td style="text-align:right;">

\-0.0150125

</td>

<td style="text-align:right;">

0.0131633

</td>

<td style="text-align:right;">

\-0.0015559

</td>

<td style="text-align:right;">

0.0043453

</td>

<td style="text-align:right;">

\-0.0171223

</td>

<td style="text-align:right;">

0.0260386

</td>

<td style="text-align:right;">

0.1225932

</td>

<td style="text-align:right;">

0.0037405

</td>

<td style="text-align:right;">

0.0189755

</td>

<td style="text-align:right;">

\-0.0077206

</td>

<td style="text-align:right;">

\-0.0251427

</td>

<td style="text-align:right;">

0.0219994

</td>

<td style="text-align:right;">

\-0.0027912

</td>

<td style="text-align:right;">

0.0073362

</td>

<td style="text-align:right;">

\-0.0016336

</td>

<td style="text-align:right;">

0.0014212

</td>

</tr>

<tr>

<td style="text-align:left;">

d\_output\_4

</td>

<td style="text-align:right;">

\-0.0115148

</td>

<td style="text-align:right;">

\-0.0155222

</td>

<td style="text-align:right;">

0.0018741

</td>

<td style="text-align:right;">

0.0538145

</td>

<td style="text-align:right;">

\-0.0020646

</td>

<td style="text-align:right;">

0.0030748

</td>

<td style="text-align:right;">

\-0.0087192

</td>

<td style="text-align:right;">

0.0141244

</td>

<td style="text-align:right;">

0.0141004

</td>

<td style="text-align:right;">

\-0.0121362

</td>

<td style="text-align:right;">

0.0090350

</td>

<td style="text-align:right;">

0.0119581

</td>

<td style="text-align:right;">

\-0.0015009

</td>

<td style="text-align:right;">

\-0.0416657

</td>

<td style="text-align:right;">

0.0016067

</td>

<td style="text-align:right;">

\-0.0023328

</td>

<td style="text-align:right;">

0.0067214

</td>

<td style="text-align:right;">

\-0.0109202

</td>

<td style="text-align:right;">

\-0.0108953

</td>

<td style="text-align:right;">

0.0093877

</td>

<td style="text-align:right;">

0.0087038

</td>

<td style="text-align:right;">

0.0117329

</td>

<td style="text-align:right;">

\-0.0014675

</td>

<td style="text-align:right;">

\-0.0402955

</td>

<td style="text-align:right;">

0.0015294

</td>

<td style="text-align:right;">

\-0.0022764

</td>

<td style="text-align:right;">

0.0065180

</td>

<td style="text-align:right;">

\-0.0105615

</td>

<td style="text-align:right;">

\-0.0105698

</td>

<td style="text-align:right;">

0.0090910

</td>

<td style="text-align:right;">

0.0111127

</td>

<td style="text-align:right;">

0.0150710

</td>

<td style="text-align:right;">

\-0.0017092

</td>

<td style="text-align:right;">

\-0.0520296

</td>

<td style="text-align:right;">

0.0019920

</td>

<td style="text-align:right;">

\-0.0029658

</td>

<td style="text-align:right;">

0.0084436

</td>

<td style="text-align:right;">

\-0.0136627

</td>

<td style="text-align:right;">

\-0.0136540

</td>

<td style="text-align:right;">

0.0117369

</td>

<td style="text-align:right;">

0.0144505

</td>

<td style="text-align:right;">

0.0195160

</td>

<td style="text-align:right;">

\-0.0023350

</td>

<td style="text-align:right;">

\-0.0675164

</td>

<td style="text-align:right;">

0.0026036

</td>

<td style="text-align:right;">

\-0.0038820

</td>

<td style="text-align:right;">

0.0109700

</td>

<td style="text-align:right;">

\-0.0177700

</td>

<td style="text-align:right;">

\-0.0177239

</td>

<td style="text-align:right;">

0.0152723

</td>

<td style="text-align:right;">

0.0031468

</td>

<td style="text-align:right;">

0.0042670

</td>

<td style="text-align:right;">

\-0.0004861

</td>

<td style="text-align:right;">

\-0.0148352

</td>

<td style="text-align:right;">

0.0006424

</td>

<td style="text-align:right;">

\-0.0008148

</td>

<td style="text-align:right;">

0.0024087

</td>

<td style="text-align:right;">

\-0.0038815

</td>

<td style="text-align:right;">

\-0.0038776

</td>

<td style="text-align:right;">

0.0033266

</td>

<td style="text-align:right;">

\-0.0085681

</td>

<td style="text-align:right;">

\-0.0115627

</td>

<td style="text-align:right;">

0.0013937

</td>

<td style="text-align:right;">

0.0396951

</td>

<td style="text-align:right;">

\-0.0015198

</td>

<td style="text-align:right;">

0.0023721

</td>

<td style="text-align:right;">

\-0.0064352

</td>

<td style="text-align:right;">

0.0104493

</td>

<td style="text-align:right;">

0.0104159

</td>

<td style="text-align:right;">

\-0.0090011

</td>

<td style="text-align:right;">

0.0206229

</td>

<td style="text-align:right;">

0.0278116

</td>

<td style="text-align:right;">

\-0.0034118

</td>

<td style="text-align:right;">

\-0.0961119

</td>

<td style="text-align:right;">

0.0036398

</td>

<td style="text-align:right;">

\-0.0055483

</td>

<td style="text-align:right;">

0.0156814

</td>

<td style="text-align:right;">

\-0.0252497

</td>

<td style="text-align:right;">

\-0.0252152

</td>

<td style="text-align:right;">

0.0216987

</td>

<td style="text-align:right;">

\-0.0021391

</td>

<td style="text-align:right;">

\-0.0027703

</td>

<td style="text-align:right;">

0.0003903

</td>

<td style="text-align:right;">

0.0095449

</td>

<td style="text-align:right;">

\-0.0003900

</td>

<td style="text-align:right;">

0.0005380

</td>

<td style="text-align:right;">

\-0.0015780

</td>

<td style="text-align:right;">

0.0026081

</td>

<td style="text-align:right;">

0.0025117

</td>

<td style="text-align:right;">

\-0.0021678

</td>

<td style="text-align:right;">

\-0.0139711

</td>

<td style="text-align:right;">

\-0.0188096

</td>

<td style="text-align:right;">

0.0022905

</td>

<td style="text-align:right;">

0.0650464

</td>

<td style="text-align:right;">

\-0.0024392

</td>

<td style="text-align:right;">

0.0037289

</td>

<td style="text-align:right;">

\-0.0105230

</td>

<td style="text-align:right;">

0.0170629

</td>

<td style="text-align:right;">

0.0171256

</td>

<td style="text-align:right;">

\-0.0146613

</td>

<td style="text-align:right;">

\-0.0232066

</td>

<td style="text-align:right;">

\-0.0313395

</td>

<td style="text-align:right;">

0.0037780

</td>

<td style="text-align:right;">

0.1086399

</td>

<td style="text-align:right;">

\-0.0041298

</td>

<td style="text-align:right;">

0.0062744

</td>

<td style="text-align:right;">

\-0.0176437

</td>

<td style="text-align:right;">

0.0285431

</td>

<td style="text-align:right;">

0.0284539

</td>

<td style="text-align:right;">

\-0.0244253

</td>

<td style="text-align:right;">

\-0.0016562

</td>

<td style="text-align:right;">

0.0014197

</td>

</tr>

<tr>

<td style="text-align:left;">

d\_output\_5

</td>

<td style="text-align:right;">

\-0.0002015

</td>

<td style="text-align:right;">

0.0011602

</td>

<td style="text-align:right;">

0.0093305

</td>

<td style="text-align:right;">

\-0.0021301

</td>

<td style="text-align:right;">

0.0416820

</td>

<td style="text-align:right;">

0.0063474

</td>

<td style="text-align:right;">

\-0.0124201

</td>

<td style="text-align:right;">

\-0.0084977

</td>

<td style="text-align:right;">

0.0127248

</td>

<td style="text-align:right;">

0.0140695

</td>

<td style="text-align:right;">

0.0000392

</td>

<td style="text-align:right;">

\-0.0009037

</td>

<td style="text-align:right;">

\-0.0072236

</td>

<td style="text-align:right;">

0.0016548

</td>

<td style="text-align:right;">

\-0.0322917

</td>

<td style="text-align:right;">

\-0.0049382

</td>

<td style="text-align:right;">

0.0096484

</td>

<td style="text-align:right;">

0.0065554

</td>

<td style="text-align:right;">

\-0.0098702

</td>

<td style="text-align:right;">

\-0.0108998

</td>

<td style="text-align:right;">

0.0000799

</td>

<td style="text-align:right;">

\-0.0010397

</td>

<td style="text-align:right;">

\-0.0070015

</td>

<td style="text-align:right;">

0.0015956

</td>

<td style="text-align:right;">

\-0.0311980

</td>

<td style="text-align:right;">

\-0.0047749

</td>

<td style="text-align:right;">

0.0093385

</td>

<td style="text-align:right;">

0.0063402

</td>

<td style="text-align:right;">

\-0.0095232

</td>

<td style="text-align:right;">

\-0.0105459

</td>

<td style="text-align:right;">

0.0002337

</td>

<td style="text-align:right;">

\-0.0011615

</td>

<td style="text-align:right;">

\-0.0091311

</td>

<td style="text-align:right;">

0.0019994

</td>

<td style="text-align:right;">

\-0.0403015

</td>

<td style="text-align:right;">

\-0.0061458

</td>

<td style="text-align:right;">

0.0119972

</td>

<td style="text-align:right;">

0.0082189

</td>

<td style="text-align:right;">

\-0.0123143

</td>

<td style="text-align:right;">

\-0.0135727

</td>

<td style="text-align:right;">

0.0002891

</td>

<td style="text-align:right;">

\-0.0014201

</td>

<td style="text-align:right;">

\-0.0117045

</td>

<td style="text-align:right;">

0.0025540

</td>

<td style="text-align:right;">

\-0.0522549

</td>

<td style="text-align:right;">

\-0.0079484

</td>

<td style="text-align:right;">

0.0155687

</td>

<td style="text-align:right;">

0.0106572

</td>

<td style="text-align:right;">

\-0.0159636

</td>

<td style="text-align:right;">

\-0.0176338

</td>

<td style="text-align:right;">

0.0000658

</td>

<td style="text-align:right;">

\-0.0003473

</td>

<td style="text-align:right;">

\-0.0025924

</td>

<td style="text-align:right;">

0.0005953

</td>

<td style="text-align:right;">

\-0.0115158

</td>

<td style="text-align:right;">

\-0.0017385

</td>

<td style="text-align:right;">

0.0033803

</td>

<td style="text-align:right;">

0.0023231

</td>

<td style="text-align:right;">

\-0.0034806

</td>

<td style="text-align:right;">

\-0.0038647

</td>

<td style="text-align:right;">

\-0.0000968

</td>

<td style="text-align:right;">

0.0009752

</td>

<td style="text-align:right;">

0.0069060

</td>

<td style="text-align:right;">

\-0.0015688

</td>

<td style="text-align:right;">

0.0307289

</td>

<td style="text-align:right;">

0.0045695

</td>

<td style="text-align:right;">

\-0.0091735

</td>

<td style="text-align:right;">

\-0.0062613

</td>

<td style="text-align:right;">

0.0093726

</td>

<td style="text-align:right;">

0.0104121

</td>

<td style="text-align:right;">

0.0003115

</td>

<td style="text-align:right;">

\-0.0021644

</td>

<td style="text-align:right;">

\-0.0166177

</td>

<td style="text-align:right;">

0.0038501

</td>

<td style="text-align:right;">

\-0.0743845

</td>

<td style="text-align:right;">

\-0.0112928

</td>

<td style="text-align:right;">

0.0220746

</td>

<td style="text-align:right;">

0.0151891

</td>

<td style="text-align:right;">

\-0.0226929

</td>

<td style="text-align:right;">

\-0.0251613

</td>

<td style="text-align:right;">

0.0000843

</td>

<td style="text-align:right;">

0.0002798

</td>

<td style="text-align:right;">

0.0016057

</td>

<td style="text-align:right;">

\-0.0004489

</td>

<td style="text-align:right;">

0.0073900

</td>

<td style="text-align:right;">

0.0011070

</td>

<td style="text-align:right;">

\-0.0021410

</td>

<td style="text-align:right;">

\-0.0016478

</td>

<td style="text-align:right;">

0.0022166

</td>

<td style="text-align:right;">

0.0025215

</td>

<td style="text-align:right;">

\-0.0001848

</td>

<td style="text-align:right;">

0.0014832

</td>

<td style="text-align:right;">

0.0112817

</td>

<td style="text-align:right;">

\-0.0026156

</td>

<td style="text-align:right;">

0.0503344

</td>

<td style="text-align:right;">

0.0076649

</td>

<td style="text-align:right;">

\-0.0150416

</td>

<td style="text-align:right;">

\-0.0102730

</td>

<td style="text-align:right;">

0.0152763

</td>

<td style="text-align:right;">

0.0170167

</td>

<td style="text-align:right;">

\-0.0004861

</td>

<td style="text-align:right;">

0.0022518

</td>

<td style="text-align:right;">

0.0188271

</td>

<td style="text-align:right;">

\-0.0042536

</td>

<td style="text-align:right;">

0.0841183

</td>

<td style="text-align:right;">

0.0127844

</td>

<td style="text-align:right;">

\-0.0250641

</td>

<td style="text-align:right;">

\-0.0171571

</td>

<td style="text-align:right;">

0.0257523

</td>

<td style="text-align:right;">

0.0282801

</td>

<td style="text-align:right;">

0.0021382

</td>

<td style="text-align:right;">

\-0.0018787

</td>

</tr>

<tr>

<td style="text-align:left;">

d\_output\_6

</td>

<td style="text-align:right;">

\-0.0098928

</td>

<td style="text-align:right;">

\-0.0297162

</td>

<td style="text-align:right;">

\-0.0037695

</td>

<td style="text-align:right;">

0.0030480

</td>

<td style="text-align:right;">

0.0063663

</td>

<td style="text-align:right;">

0.0552486

</td>

<td style="text-align:right;">

\-0.0143575

</td>

<td style="text-align:right;">

0.0052084

</td>

<td style="text-align:right;">

0.0107802

</td>

<td style="text-align:right;">

0.0029326

</td>

<td style="text-align:right;">

0.0074214

</td>

<td style="text-align:right;">

0.0230069

</td>

<td style="text-align:right;">

0.0029022

</td>

<td style="text-align:right;">

\-0.0022937

</td>

<td style="text-align:right;">

\-0.0049406

</td>

<td style="text-align:right;">

\-0.0428277

</td>

<td style="text-align:right;">

0.0111680

</td>

<td style="text-align:right;">

\-0.0040983

</td>

<td style="text-align:right;">

\-0.0083541

</td>

<td style="text-align:right;">

\-0.0022415

</td>

<td style="text-align:right;">

0.0072248

</td>

<td style="text-align:right;">

0.0217703

</td>

<td style="text-align:right;">

0.0027583

</td>

<td style="text-align:right;">

\-0.0020931

</td>

<td style="text-align:right;">

\-0.0046577

</td>

<td style="text-align:right;">

\-0.0413635

</td>

<td style="text-align:right;">

0.0107642

</td>

<td style="text-align:right;">

\-0.0039331

</td>

<td style="text-align:right;">

\-0.0079740

</td>

<td style="text-align:right;">

\-0.0022010

</td>

<td style="text-align:right;">

0.0096220

</td>

<td style="text-align:right;">

0.0285772

</td>

<td style="text-align:right;">

0.0033567

</td>

<td style="text-align:right;">

\-0.0029425

</td>

<td style="text-align:right;">

\-0.0061197

</td>

<td style="text-align:right;">

\-0.0534077

</td>

<td style="text-align:right;">

0.0138521

</td>

<td style="text-align:right;">

\-0.0050280

</td>

<td style="text-align:right;">

\-0.0103800

</td>

<td style="text-align:right;">

\-0.0028076

</td>

<td style="text-align:right;">

0.0124796

</td>

<td style="text-align:right;">

0.0373442

</td>

<td style="text-align:right;">

0.0047290

</td>

<td style="text-align:right;">

\-0.0041209

</td>

<td style="text-align:right;">

\-0.0080016

</td>

<td style="text-align:right;">

\-0.0692467

</td>

<td style="text-align:right;">

0.0179884

</td>

<td style="text-align:right;">

\-0.0065145

</td>

<td style="text-align:right;">

\-0.0135446

</td>

<td style="text-align:right;">

\-0.0036905

</td>

<td style="text-align:right;">

0.0027195

</td>

<td style="text-align:right;">

0.0080265

</td>

<td style="text-align:right;">

0.0009267

</td>

<td style="text-align:right;">

\-0.0007802

</td>

<td style="text-align:right;">

\-0.0019473

</td>

<td style="text-align:right;">

\-0.0150756

</td>

<td style="text-align:right;">

0.0038894

</td>

<td style="text-align:right;">

\-0.0014342

</td>

<td style="text-align:right;">

\-0.0029022

</td>

<td style="text-align:right;">

\-0.0008212

</td>

<td style="text-align:right;">

\-0.0071938

</td>

<td style="text-align:right;">

\-0.0216834

</td>

<td style="text-align:right;">

\-0.0027394

</td>

<td style="text-align:right;">

0.0022030

</td>

<td style="text-align:right;">

0.0046396

</td>

<td style="text-align:right;">

0.0405590

</td>

<td style="text-align:right;">

\-0.0106023

</td>

<td style="text-align:right;">

0.0038578

</td>

<td style="text-align:right;">

0.0079088

</td>

<td style="text-align:right;">

0.0022191

</td>

<td style="text-align:right;">

0.0175856

</td>

<td style="text-align:right;">

0.0529079

</td>

<td style="text-align:right;">

0.0068558

</td>

<td style="text-align:right;">

\-0.0053461

</td>

<td style="text-align:right;">

\-0.0112332

</td>

<td style="text-align:right;">

\-0.0985946

</td>

<td style="text-align:right;">

0.0253945

</td>

<td style="text-align:right;">

\-0.0093039

</td>

<td style="text-align:right;">

\-0.0191869

</td>

<td style="text-align:right;">

\-0.0052745

</td>

<td style="text-align:right;">

\-0.0015026

</td>

<td style="text-align:right;">

\-0.0050982

</td>

<td style="text-align:right;">

\-0.0007503

</td>

<td style="text-align:right;">

0.0003666

</td>

<td style="text-align:right;">

0.0011963

</td>

<td style="text-align:right;">

0.0097106

</td>

<td style="text-align:right;">

\-0.0024161

</td>

<td style="text-align:right;">

0.0006043

</td>

<td style="text-align:right;">

0.0018568

</td>

<td style="text-align:right;">

0.0005551

</td>

<td style="text-align:right;">

\-0.0118376

</td>

<td style="text-align:right;">

\-0.0357508

</td>

<td style="text-align:right;">

\-0.0045558

</td>

<td style="text-align:right;">

0.0035470

</td>

<td style="text-align:right;">

0.0075299

</td>

<td style="text-align:right;">

0.0667830

</td>

<td style="text-align:right;">

\-0.0173927

</td>

<td style="text-align:right;">

0.0063052

</td>

<td style="text-align:right;">

0.0127777

</td>

<td style="text-align:right;">

0.0035293

</td>

<td style="text-align:right;">

\-0.0201816

</td>

<td style="text-align:right;">

\-0.0602960

</td>

<td style="text-align:right;">

\-0.0076860

</td>

<td style="text-align:right;">

0.0063388

</td>

<td style="text-align:right;">

0.0127583

</td>

<td style="text-align:right;">

0.1115422

</td>

<td style="text-align:right;">

\-0.0289714

</td>

<td style="text-align:right;">

0.0105507

</td>

<td style="text-align:right;">

0.0219431

</td>

<td style="text-align:right;">

0.0056070

</td>

<td style="text-align:right;">

0.0048153

</td>

<td style="text-align:right;">

\-0.0042268

</td>

</tr>

<tr>

<td style="text-align:left;">

d\_output\_7

</td>

<td style="text-align:right;">

0.0069446

</td>

<td style="text-align:right;">

0.0192320

</td>

<td style="text-align:right;">

\-0.0124953

</td>

<td style="text-align:right;">

\-0.0086532

</td>

<td style="text-align:right;">

\-0.0124090

</td>

<td style="text-align:right;">

\-0.0143473

</td>

<td style="text-align:right;">

0.0558417

</td>

<td style="text-align:right;">

\-0.0146184

</td>

<td style="text-align:right;">

\-0.0051602

</td>

<td style="text-align:right;">

\-0.0035448

</td>

<td style="text-align:right;">

\-0.0051272

</td>

<td style="text-align:right;">

\-0.0148759

</td>

<td style="text-align:right;">

0.0097013

</td>

<td style="text-align:right;">

0.0066371

</td>

<td style="text-align:right;">

0.0096298

</td>

<td style="text-align:right;">

0.0111660

</td>

<td style="text-align:right;">

\-0.0432957

</td>

<td style="text-align:right;">

0.0113787

</td>

<td style="text-align:right;">

0.0040087

</td>

<td style="text-align:right;">

0.0027347

</td>

<td style="text-align:right;">

\-0.0050374

</td>

<td style="text-align:right;">

\-0.0139463

</td>

<td style="text-align:right;">

0.0094648

</td>

<td style="text-align:right;">

0.0063540

</td>

<td style="text-align:right;">

0.0092342

</td>

<td style="text-align:right;">

0.0107642

</td>

<td style="text-align:right;">

\-0.0418240

</td>

<td style="text-align:right;">

0.0109867

</td>

<td style="text-align:right;">

0.0038158

</td>

<td style="text-align:right;">

0.0026658

</td>

<td style="text-align:right;">

\-0.0067943

</td>

<td style="text-align:right;">

\-0.0184984

</td>

<td style="text-align:right;">

0.0123624

</td>

<td style="text-align:right;">

0.0084440

</td>

<td style="text-align:right;">

0.0120081

</td>

<td style="text-align:right;">

0.0139189

</td>

<td style="text-align:right;">

\-0.0540180

</td>

<td style="text-align:right;">

0.0141399

</td>

<td style="text-align:right;">

0.0049925

</td>

<td style="text-align:right;">

0.0034015

</td>

<td style="text-align:right;">

\-0.0087916

</td>

<td style="text-align:right;">

\-0.0242142

</td>

<td style="text-align:right;">

0.0156719

</td>

<td style="text-align:right;">

0.0111644

</td>

<td style="text-align:right;">

0.0156054

</td>

<td style="text-align:right;">

0.0179926

</td>

<td style="text-align:right;">

\-0.0700294

</td>

<td style="text-align:right;">

0.0183071

</td>

<td style="text-align:right;">

0.0065053

</td>

<td style="text-align:right;">

0.0044500

</td>

<td style="text-align:right;">

\-0.0019214

</td>

<td style="text-align:right;">

\-0.0051924

</td>

<td style="text-align:right;">

0.0035041

</td>

<td style="text-align:right;">

0.0023222

</td>

<td style="text-align:right;">

0.0036054

</td>

<td style="text-align:right;">

0.0039501

</td>

<td style="text-align:right;">

\-0.0152366

</td>

<td style="text-align:right;">

0.0040064

</td>

<td style="text-align:right;">

0.0013745

</td>

<td style="text-align:right;">

0.0009905

</td>

<td style="text-align:right;">

0.0050195

</td>

<td style="text-align:right;">

0.0139367

</td>

<td style="text-align:right;">

\-0.0093230

</td>

<td style="text-align:right;">

\-0.0063717

</td>

<td style="text-align:right;">

\-0.0091359

</td>

<td style="text-align:right;">

\-0.0103301

</td>

<td style="text-align:right;">

0.0413188

</td>

<td style="text-align:right;">

\-0.0108259

</td>

<td style="text-align:right;">

\-0.0037585

</td>

<td style="text-align:right;">

\-0.0026864

</td>

<td style="text-align:right;">

\-0.0123397

</td>

<td style="text-align:right;">

\-0.0342099

</td>

<td style="text-align:right;">

0.0222600

</td>

<td style="text-align:right;">

0.0153796

</td>

<td style="text-align:right;">

0.0220557

</td>

<td style="text-align:right;">

0.0255470

</td>

<td style="text-align:right;">

\-0.0996328

</td>

<td style="text-align:right;">

0.0261613

</td>

<td style="text-align:right;">

0.0091591

</td>

<td style="text-align:right;">

0.0063636

</td>

<td style="text-align:right;">

0.0009729

</td>

<td style="text-align:right;">

0.0032416

</td>

<td style="text-align:right;">

\-0.0021301

</td>

<td style="text-align:right;">

\-0.0013462

</td>

<td style="text-align:right;">

\-0.0022680

</td>

<td style="text-align:right;">

\-0.0025112

</td>

<td style="text-align:right;">

0.0097469

</td>

<td style="text-align:right;">

\-0.0022555

</td>

<td style="text-align:right;">

\-0.0008595

</td>

<td style="text-align:right;">

\-0.0006834

</td>

<td style="text-align:right;">

0.0082659

</td>

<td style="text-align:right;">

0.0230507

</td>

<td style="text-align:right;">

\-0.0151265

</td>

<td style="text-align:right;">

\-0.0103304

</td>

<td style="text-align:right;">

\-0.0148400

</td>

<td style="text-align:right;">

\-0.0173030

</td>

<td style="text-align:right;">

0.0675348

</td>

<td style="text-align:right;">

\-0.0176816

</td>

<td style="text-align:right;">

\-0.0059760

</td>

<td style="text-align:right;">

\-0.0042830

</td>

<td style="text-align:right;">

0.0142285

</td>

<td style="text-align:right;">

0.0391065

</td>

<td style="text-align:right;">

\-0.0251826

</td>

<td style="text-align:right;">

\-0.0176514

</td>

<td style="text-align:right;">

\-0.0249597

</td>

<td style="text-align:right;">

\-0.0289096

</td>

<td style="text-align:right;">

0.1127846

</td>

<td style="text-align:right;">

\-0.0295681

</td>

<td style="text-align:right;">

\-0.0105792

</td>

<td style="text-align:right;">

\-0.0068367

</td>

<td style="text-align:right;">

\-0.0050240

</td>

<td style="text-align:right;">

0.0044264

</td>

</tr>

<tr>

<td style="text-align:left;">

d\_output\_8

</td>

<td style="text-align:right;">

\-0.0250804

</td>

<td style="text-align:right;">

\-0.0161938

</td>

<td style="text-align:right;">

0.0109440

</td>

<td style="text-align:right;">

0.0141087

</td>

<td style="text-align:right;">

\-0.0084907

</td>

<td style="text-align:right;">

0.0051460

</td>

<td style="text-align:right;">

\-0.0146053

</td>

<td style="text-align:right;">

0.0673064

</td>

<td style="text-align:right;">

0.0047617

</td>

<td style="text-align:right;">

\-0.0103918

</td>

<td style="text-align:right;">

0.0193901

</td>

<td style="text-align:right;">

0.0123672

</td>

<td style="text-align:right;">

\-0.0085698

</td>

<td style="text-align:right;">

\-0.0107848

</td>

<td style="text-align:right;">

0.0066224

</td>

<td style="text-align:right;">

\-0.0038823

</td>

<td style="text-align:right;">

0.0111201

</td>

<td style="text-align:right;">

\-0.0522336

</td>

<td style="text-align:right;">

\-0.0036282

</td>

<td style="text-align:right;">

0.0080575

</td>

<td style="text-align:right;">

0.0187944

</td>

<td style="text-align:right;">

0.0119589

</td>

<td style="text-align:right;">

\-0.0082823

</td>

<td style="text-align:right;">

\-0.0105283

</td>

<td style="text-align:right;">

0.0064342

</td>

<td style="text-align:right;">

\-0.0038049

</td>

<td style="text-align:right;">

0.0108549

</td>

<td style="text-align:right;">

\-0.0505127

</td>

<td style="text-align:right;">

\-0.0035369

</td>

<td style="text-align:right;">

0.0078282

</td>

<td style="text-align:right;">

0.0242659

</td>

<td style="text-align:right;">

0.0156709

</td>

<td style="text-align:right;">

\-0.0106651

</td>

<td style="text-align:right;">

\-0.0137003

</td>

<td style="text-align:right;">

0.0081853

</td>

<td style="text-align:right;">

\-0.0050257

</td>

<td style="text-align:right;">

0.0142004

</td>

<td style="text-align:right;">

\-0.0650318

</td>

<td style="text-align:right;">

\-0.0046159

</td>

<td style="text-align:right;">

0.0100563

</td>

<td style="text-align:right;">

0.0314683

</td>

<td style="text-align:right;">

0.0204053

</td>

<td style="text-align:right;">

\-0.0137138

</td>

<td style="text-align:right;">

\-0.0178233

</td>

<td style="text-align:right;">

0.0105922

</td>

<td style="text-align:right;">

\-0.0065138

</td>

<td style="text-align:right;">

0.0184283

</td>

<td style="text-align:right;">

\-0.0843316

</td>

<td style="text-align:right;">

\-0.0060171

</td>

<td style="text-align:right;">

0.0130294

</td>

<td style="text-align:right;">

0.0068385

</td>

<td style="text-align:right;">

0.0043584

</td>

<td style="text-align:right;">

\-0.0030011

</td>

<td style="text-align:right;">

\-0.0038238

</td>

<td style="text-align:right;">

0.0022621

</td>

<td style="text-align:right;">

\-0.0013739

</td>

<td style="text-align:right;">

0.0038935

</td>

<td style="text-align:right;">

\-0.0183673

</td>

<td style="text-align:right;">

\-0.0012828

</td>

<td style="text-align:right;">

0.0028227

</td>

<td style="text-align:right;">

\-0.0185725

</td>

<td style="text-align:right;">

\-0.0119098

</td>

<td style="text-align:right;">

0.0081547

</td>

<td style="text-align:right;">

0.0104460

</td>

<td style="text-align:right;">

\-0.0063207

</td>

<td style="text-align:right;">

0.0037169

</td>

<td style="text-align:right;">

\-0.0108092

</td>

<td style="text-align:right;">

0.0498851

</td>

<td style="text-align:right;">

0.0035081

</td>

<td style="text-align:right;">

\-0.0076994

</td>

<td style="text-align:right;">

0.0449804

</td>

<td style="text-align:right;">

0.0289774

</td>

<td style="text-align:right;">

\-0.0196031

</td>

<td style="text-align:right;">

\-0.0252719

</td>

<td style="text-align:right;">

0.0153021

</td>

<td style="text-align:right;">

\-0.0091713

</td>

<td style="text-align:right;">

0.0260506

</td>

<td style="text-align:right;">

\-0.1207382

</td>

<td style="text-align:right;">

\-0.0085237

</td>

<td style="text-align:right;">

0.0186647

</td>

<td style="text-align:right;">

\-0.0042959

</td>

<td style="text-align:right;">

\-0.0026810

</td>

<td style="text-align:right;">

0.0019020

</td>

<td style="text-align:right;">

0.0023368

</td>

<td style="text-align:right;">

\-0.0014871

</td>

<td style="text-align:right;">

0.0008146

</td>

<td style="text-align:right;">

\-0.0023603

</td>

<td style="text-align:right;">

0.0116123

</td>

<td style="text-align:right;">

0.0007876

</td>

<td style="text-align:right;">

\-0.0018047

</td>

<td style="text-align:right;">

\-0.0303050

</td>

<td style="text-align:right;">

\-0.0194974

</td>

<td style="text-align:right;">

0.0132639

</td>

<td style="text-align:right;">

0.0170032

</td>

<td style="text-align:right;">

\-0.0103353

</td>

<td style="text-align:right;">

0.0061868

</td>

<td style="text-align:right;">

\-0.0176369

</td>

<td style="text-align:right;">

0.0814000

</td>

<td style="text-align:right;">

0.0056628

</td>

<td style="text-align:right;">

\-0.0125827

</td>

<td style="text-align:right;">

\-0.0507465

</td>

<td style="text-align:right;">

\-0.0328884

</td>

<td style="text-align:right;">

0.0220748

</td>

<td style="text-align:right;">

0.0286126

</td>

<td style="text-align:right;">

\-0.0171680

</td>

<td style="text-align:right;">

0.0104383

</td>

<td style="text-align:right;">

\-0.0296133

</td>

<td style="text-align:right;">

0.1359903

</td>

<td style="text-align:right;">

0.0097023

</td>

<td style="text-align:right;">

\-0.0211098

</td>

<td style="text-align:right;">

0.0015298

</td>

<td style="text-align:right;">

\-0.0013439

</td>

</tr>

<tr>

<td style="text-align:left;">

d\_output\_9

</td>

<td style="text-align:right;">

\-0.0138377

</td>

<td style="text-align:right;">

\-0.0266094

</td>

<td style="text-align:right;">

\-0.0013464

</td>

<td style="text-align:right;">

0.0141003

</td>

<td style="text-align:right;">

0.0127371

</td>

<td style="text-align:right;">

0.0108331

</td>

<td style="text-align:right;">

\-0.0052421

</td>

<td style="text-align:right;">

0.0048221

</td>

<td style="text-align:right;">

0.0530128

</td>

<td style="text-align:right;">

\-0.0142592

</td>

<td style="text-align:right;">

0.0106016

</td>

<td style="text-align:right;">

0.0204647

</td>

<td style="text-align:right;">

0.0009927

</td>

<td style="text-align:right;">

\-0.0108144

</td>

<td style="text-align:right;">

\-0.0098145

</td>

<td style="text-align:right;">

\-0.0083304

</td>

<td style="text-align:right;">

0.0039921

</td>

<td style="text-align:right;">

\-0.0038836

</td>

<td style="text-align:right;">

\-0.0410147

</td>

<td style="text-align:right;">

0.0110594

</td>

<td style="text-align:right;">

0.0103596

</td>

<td style="text-align:right;">

0.0196407

</td>

<td style="text-align:right;">

0.0009328

</td>

<td style="text-align:right;">

\-0.0104557

</td>

<td style="text-align:right;">

\-0.0094673

</td>

<td style="text-align:right;">

\-0.0080459

</td>

<td style="text-align:right;">

0.0037988

</td>

<td style="text-align:right;">

\-0.0038664

</td>

<td style="text-align:right;">

\-0.0396832

</td>

<td style="text-align:right;">

0.0107345

</td>

<td style="text-align:right;">

0.0134222

</td>

<td style="text-align:right;">

0.0256635

</td>

<td style="text-align:right;">

0.0011410

</td>

<td style="text-align:right;">

\-0.0136421

</td>

<td style="text-align:right;">

\-0.0123035

</td>

<td style="text-align:right;">

\-0.0104722

</td>

<td style="text-align:right;">

0.0050322

</td>

<td style="text-align:right;">

\-0.0046676

</td>

<td style="text-align:right;">

\-0.0512514

</td>

<td style="text-align:right;">

0.0138072

</td>

<td style="text-align:right;">

0.0173695

</td>

<td style="text-align:right;">

0.0334742

</td>

<td style="text-align:right;">

0.0017008

</td>

<td style="text-align:right;">

\-0.0178857

</td>

<td style="text-align:right;">

\-0.0160227

</td>

<td style="text-align:right;">

\-0.0136459

</td>

<td style="text-align:right;">

0.0066662

</td>

<td style="text-align:right;">

\-0.0059179

</td>

<td style="text-align:right;">

\-0.0664936

</td>

<td style="text-align:right;">

0.0178668

</td>

<td style="text-align:right;">

0.0037365

</td>

<td style="text-align:right;">

0.0072135

</td>

<td style="text-align:right;">

0.0003140

</td>

<td style="text-align:right;">

\-0.0037982

</td>

<td style="text-align:right;">

\-0.0036054

</td>

<td style="text-align:right;">

\-0.0029823

</td>

<td style="text-align:right;">

0.0015032

</td>

<td style="text-align:right;">

\-0.0012083

</td>

<td style="text-align:right;">

\-0.0144189

</td>

<td style="text-align:right;">

0.0038367

</td>

<td style="text-align:right;">

\-0.0102321

</td>

<td style="text-align:right;">

\-0.0195538

</td>

<td style="text-align:right;">

\-0.0009581

</td>

<td style="text-align:right;">

0.0104256

</td>

<td style="text-align:right;">

0.0093853

</td>

<td style="text-align:right;">

0.0078496

</td>

<td style="text-align:right;">

\-0.0038132

</td>

<td style="text-align:right;">

0.0036739

</td>

<td style="text-align:right;">

0.0392466

</td>

<td style="text-align:right;">

\-0.0105759

</td>

<td style="text-align:right;">

0.0247708

</td>

<td style="text-align:right;">

0.0475729

</td>

<td style="text-align:right;">

0.0024715

</td>

<td style="text-align:right;">

\-0.0252144

</td>

<td style="text-align:right;">

\-0.0227486

</td>

<td style="text-align:right;">

\-0.0193319

</td>

<td style="text-align:right;">

0.0092072

</td>

<td style="text-align:right;">

\-0.0086815

</td>

<td style="text-align:right;">

\-0.0949740

</td>

<td style="text-align:right;">

0.0255564

</td>

<td style="text-align:right;">

\-0.0021981

</td>

<td style="text-align:right;">

\-0.0043490

</td>

<td style="text-align:right;">

\-0.0002609

</td>

<td style="text-align:right;">

0.0022409

</td>

<td style="text-align:right;">

0.0021664

</td>

<td style="text-align:right;">

0.0017800

</td>

<td style="text-align:right;">

\-0.0007949

</td>

<td style="text-align:right;">

0.0006660

</td>

<td style="text-align:right;">

0.0088741

</td>

<td style="text-align:right;">

\-0.0023693

</td>

<td style="text-align:right;">

\-0.0167036

</td>

<td style="text-align:right;">

\-0.0320287

</td>

<td style="text-align:right;">

\-0.0016078

</td>

<td style="text-align:right;">

0.0169591

</td>

<td style="text-align:right;">

0.0152770

</td>

<td style="text-align:right;">

0.0130280

</td>

<td style="text-align:right;">

\-0.0062642

</td>

<td style="text-align:right;">

0.0059658

</td>

<td style="text-align:right;">

0.0639315

</td>

<td style="text-align:right;">

\-0.0172731

</td>

<td style="text-align:right;">

\-0.0280732

</td>

<td style="text-align:right;">

\-0.0539748

</td>

<td style="text-align:right;">

\-0.0027756

</td>

<td style="text-align:right;">

0.0286292

</td>

<td style="text-align:right;">

0.0256939

</td>

<td style="text-align:right;">

0.0218936

</td>

<td style="text-align:right;">

\-0.0106074

</td>

<td style="text-align:right;">

0.0097172

</td>

<td style="text-align:right;">

0.1072256

</td>

<td style="text-align:right;">

\-0.0290052

</td>

<td style="text-align:right;">

0.0027076

</td>

<td style="text-align:right;">

\-0.0023674

</td>

</tr>

<tr>

<td style="text-align:left;">

d\_output\_10

</td>

<td style="text-align:right;">

0.0183998

</td>

<td style="text-align:right;">

0.0235696

</td>

<td style="text-align:right;">

0.0035820

</td>

<td style="text-align:right;">

\-0.0121090

</td>

<td style="text-align:right;">

0.0140323

</td>

<td style="text-align:right;">

0.0028750

</td>

<td style="text-align:right;">

\-0.0034804

</td>

<td style="text-align:right;">

\-0.0103564

</td>

<td style="text-align:right;">

\-0.0142733

</td>

<td style="text-align:right;">

0.0604927

</td>

<td style="text-align:right;">

\-0.0145180

</td>

<td style="text-align:right;">

\-0.0181169

</td>

<td style="text-align:right;">

\-0.0026837

</td>

<td style="text-align:right;">

0.0092915

</td>

<td style="text-align:right;">

\-0.0108873

</td>

<td style="text-align:right;">

\-0.0023307

</td>

<td style="text-align:right;">

0.0028632

</td>

<td style="text-align:right;">

0.0081400

</td>

<td style="text-align:right;">

0.0110213

</td>

<td style="text-align:right;">

\-0.0467783

</td>

<td style="text-align:right;">

\-0.0140121

</td>

<td style="text-align:right;">

\-0.0179260

</td>

<td style="text-align:right;">

\-0.0025975

</td>

<td style="text-align:right;">

0.0089999

</td>

<td style="text-align:right;">

\-0.0105325

</td>

<td style="text-align:right;">

\-0.0022815

</td>

<td style="text-align:right;">

0.0028046

</td>

<td style="text-align:right;">

0.0079632

</td>

<td style="text-align:right;">

0.0107114

</td>

<td style="text-align:right;">

\-0.0453199

</td>

<td style="text-align:right;">

\-0.0177544

</td>

<td style="text-align:right;">

\-0.0229244

</td>

<td style="text-align:right;">

\-0.0036911

</td>

<td style="text-align:right;">

0.0116933

</td>

<td style="text-align:right;">

\-0.0135718

</td>

<td style="text-align:right;">

\-0.0027779

</td>

<td style="text-align:right;">

0.0033604

</td>

<td style="text-align:right;">

0.0100647

</td>

<td style="text-align:right;">

0.0138165

</td>

<td style="text-align:right;">

\-0.0584723

</td>

<td style="text-align:right;">

\-0.0229918

</td>

<td style="text-align:right;">

\-0.0295581

</td>

<td style="text-align:right;">

\-0.0045324

</td>

<td style="text-align:right;">

0.0150177

</td>

<td style="text-align:right;">

\-0.0176017

</td>

<td style="text-align:right;">

\-0.0035554

</td>

<td style="text-align:right;">

0.0042804

</td>

<td style="text-align:right;">

0.0129093

</td>

<td style="text-align:right;">

0.0178981

</td>

<td style="text-align:right;">

\-0.0758737

</td>

<td style="text-align:right;">

\-0.0050233

</td>

<td style="text-align:right;">

\-0.0064541

</td>

<td style="text-align:right;">

\-0.0010407

</td>

<td style="text-align:right;">

0.0033427

</td>

<td style="text-align:right;">

\-0.0040061

</td>

<td style="text-align:right;">

\-0.0008135

</td>

<td style="text-align:right;">

0.0009901

</td>

<td style="text-align:right;">

0.0029355

</td>

<td style="text-align:right;">

0.0039354

</td>

<td style="text-align:right;">

\-0.0164856

</td>

<td style="text-align:right;">

0.0137477

</td>

<td style="text-align:right;">

0.0177026

</td>

<td style="text-align:right;">

0.0026546

</td>

<td style="text-align:right;">

\-0.0089431

</td>

<td style="text-align:right;">

0.0103514

</td>

<td style="text-align:right;">

0.0019110

</td>

<td style="text-align:right;">

\-0.0025516

</td>

<td style="text-align:right;">

\-0.0076565

</td>

<td style="text-align:right;">

\-0.0106108

</td>

<td style="text-align:right;">

0.0448584

</td>

<td style="text-align:right;">

\-0.0330814

</td>

<td style="text-align:right;">

\-0.0424295

</td>

<td style="text-align:right;">

\-0.0063049

</td>

<td style="text-align:right;">

0.0217503

</td>

<td style="text-align:right;">

\-0.0250430

</td>

<td style="text-align:right;">

\-0.0050663

</td>

<td style="text-align:right;">

0.0060199

</td>

<td style="text-align:right;">

0.0185659

</td>

<td style="text-align:right;">

0.0256592

</td>

<td style="text-align:right;">

\-0.1084753

</td>

<td style="text-align:right;">

0.0033399

</td>

<td style="text-align:right;">

0.0040222

</td>

<td style="text-align:right;">

0.0004825

</td>

<td style="text-align:right;">

\-0.0021145

</td>

<td style="text-align:right;">

0.0024388

</td>

<td style="text-align:right;">

0.0005097

</td>

<td style="text-align:right;">

\-0.0005756

</td>

<td style="text-align:right;">

\-0.0020969

</td>

<td style="text-align:right;">

\-0.0024352

</td>

<td style="text-align:right;">

0.0102202

</td>

<td style="text-align:right;">

0.0223668

</td>

<td style="text-align:right;">

0.0286054

</td>

<td style="text-align:right;">

0.0042865

</td>

<td style="text-align:right;">

\-0.0146718

</td>

<td style="text-align:right;">

0.0168359

</td>

<td style="text-align:right;">

0.0034689

</td>

<td style="text-align:right;">

\-0.0042721

</td>

<td style="text-align:right;">

\-0.0125683

</td>

<td style="text-align:right;">

\-0.0174525

</td>

<td style="text-align:right;">

0.0730667

</td>

<td style="text-align:right;">

0.0369782

</td>

<td style="text-align:right;">

0.0475096

</td>

<td style="text-align:right;">

0.0072398

</td>

<td style="text-align:right;">

\-0.0244064

</td>

<td style="text-align:right;">

0.0282180

</td>

<td style="text-align:right;">

0.0056742

</td>

<td style="text-align:right;">

\-0.0068803

</td>

<td style="text-align:right;">

\-0.0207674

</td>

<td style="text-align:right;">

\-0.0287255

</td>

<td style="text-align:right;">

0.1219028

</td>

<td style="text-align:right;">

0.0038724

</td>

<td style="text-align:right;">

\-0.0033642

</td>

</tr>

</tbody>

</table>

</div>
