% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mc_rv_simulation.R
\name{rgamma0,numeric,dual,dual-method}
\alias{rgamma0,numeric,dual,dual-method}
\alias{rgamma0,numeric,dual,numeric-method}
\alias{rgamma0,numeric,numeric,dual-method}
\title{Simulate gamma random variates}
\usage{
\S4method{rgamma0}{numeric,dual,dual}(n, shape, scale, method = "inv_tf")

\S4method{rgamma0}{numeric,dual,numeric}(n, shape, scale, method = "inv_tf")

\S4method{rgamma0}{numeric,numeric,dual}(n, shape, scale, method = "inv_tf")
}
\arguments{
\item{n}{Positive integer; the number of samples.}

\item{shape}{A dual number or a scalar; the shape of the gamma distribution.}

\item{scale}{A dual number or a scalar; the scale of the gamma distribution.}

\item{method}{'base' or 'inv_tf'; 'base' refers to `stats::rgamma` while
'inv_tf' refers to inverse transform.}
}
\description{
Simulate gamma random variates
}
\note{
At least one of `shape` and `scale` should be a dual number.
}
