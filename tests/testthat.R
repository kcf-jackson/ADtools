library(testthat)
library(ADtools)

test_check("ADtools")

RcppParallel::setThreadOptions(numThreads = 2)  # CRAN advice
