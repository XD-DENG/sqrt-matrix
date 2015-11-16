setwd("/Users/nus/R/sqrt-Matrix-R/")

rm(list=ls())

source('try_to_optimize_sqrtm/revised_by_XD.R')
source("Babylonian Method/Babylonian_method.R")
source("Denman Beavers Method/Denman_Beavers_iteration.R")
library(expm)

library(microbenchmark)

set.seed(3)
true_answer <- matrix(runif(10000, 10, 100), 100, 100)
test_dat <- true_answer %*% true_answer

test_dat <- rbind(cbind(1, diag(1:3)),2)

print(
  all.equal(sqrtm(test_dat) %*% sqrtm(test_dat), test_dat)
) 
# print(
#   all.equal(Babylonian_sqrtm(test_dat) %*% Babylonian_sqrtm(test_dat), test_dat)
#   )   # may encounter condition number issue
print(
  all.equal(DB_sqrtm(test_dat) %*% DB_sqrtm(test_dat), test_dat)
)
print(
  all.equal(Re(sqrtm_revised(test_dat) %*% sqrtm_revised(test_dat)), test_dat)  # here only extract the real part in case the result is made of complex numbers
)


compare_result <- microbenchmark(sqrtm(test_dat), 
               #Babylonian_sqrtm(test_dat),  # may encounter condition number issue
               DB_sqrtm(test_dat),
               sqrtm_revised(test_dat), 
               times = 100)

print(compare_result)
