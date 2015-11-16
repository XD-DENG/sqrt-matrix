# sqrt-Matrix
How to compute the square root of matrix with R language and Python.



## What have been done here
In this project, 

(1) I tried to optimize the existing R function "sqrtm" in R package "expm" with Rcpp.

(2) Introduced how to compute the square root of matrix with the existing function in Python.

(3) Implemented another two algorithms computing the square root of matrix, ***Babylonian method*** and ***Denman Beavers method*** with R language and Python separately.



## Check & Compare
I generated the code as below to compare the speed to compare the speed of these functions.

    library(microbenchmark)
    microbenchmark(sqrtm(test_dat), 
               Babylonian_sqrtm(test_dat), 
               DB_sqrtm(test_dat),
               sqrtm_revised(test_dat), 
               times = 100)


And the results are as below

![\[pic link\]](http://me.seekingqed.com/image/sqrtm_compare_results.png)

After optimzed with Rcpp, the *sqrtm* function is slower than the original *sqrtm* function which was written purely with R. This may due to the data I/O from R into C++ (my optimization was a combination of R and Rcpp).

Another finding is that the two iteration algorithms, Babylonian method and Denman Beavers Method are much faster. This may due to that *sqrtm* function considers more complicated inputs. 

(Updated on 15 November, 2015)