rm(list=ls())

source('try_to_optimize_sqrtm/revised_by_XD.R')
source("Babylonian Method/Babylonian_method.R")
source("Denman Beavers Method/Denman_Beavers_iteration.R")
library(expm)

library(microbenchmark)

test_dat <- matrix(c(3,4,2,8), 2, 2)


all.equal(sqrtm(test_dat) %*% sqrtm(test_dat), test_dat)
all.equal(Babylonian_sqrtm(test_dat) %*% Babylonian_sqrtm(test_dat), test_dat)
all.equal(DB_sqrtm(test_dat) %*% DB_sqrtm(test_dat), test_dat)
all.equal(sqrtm_revised(test_dat) %*% sqrtm_revised(test_dat), test_dat)


microbenchmark(sqrtm(test_dat), 
               Babylonian_sqrtm(test_dat), 
               DB_sqrtm(test_dat),
               sqrtm_revised(test_dat), 
               times = 100)
