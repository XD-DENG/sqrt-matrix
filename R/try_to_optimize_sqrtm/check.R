rm(list=ls())

source('try_to_optimize_sqrtm/revised_by_XD.R')
library(expm)
library(microbenchmark)

test_dat <- matrix(c(3,4,2,8), 2, 2)


all.equal(Re(sqrtm(test_dat) %*% sqrtm(test_dat)), test_dat)
all.equal(Re(sqrtm_revised(test_dat) %*% sqrtm_revised(test_dat)), test_dat)


microbenchmark(sqrtm(test_dat), sqrtm_revised(test_dat), times = 100)