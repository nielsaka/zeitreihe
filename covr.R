library(covr)
library(testthat)
# source("tests/testthat/helper-load.R")

cvr <- package_coverage()
report(cvr)
# zero_coverage(cvr)
# print(cvr, group = "functions")
# percent_coverage(cvr)
