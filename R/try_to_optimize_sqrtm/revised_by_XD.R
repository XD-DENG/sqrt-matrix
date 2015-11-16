library(Matrix)
library(Rcpp)

sourceCpp("try_to_optimize_sqrtm/optimize_1.cpp")


sqrtm_revised <- function (x) 
{
  d <- dim(x)
  if (length(d) != 2 || d[1] != d[2]) 
    stop("'x' must be a quadratic matrix")
  n <- d[1] # XD: the dimension of the matrix
  Sch.x <- Schur(Matrix(x))   # XD: Schur Decomposition of a Matrix
  ev <- Sch.x@EValues
  
  
  # XD: here the Arg() function returns the argument conjugate for complex values
  if (getOption("verbose") && any(abs(Arg(ev) - pi) < 1e-07)) 
    message(sprintf("'x' has negative real eigenvalues; maybe ok for %s", 
                    "sqrtm()"))
  
  S <- as.matrix(Sch.x@T)
  Q <- as.matrix(Sch.x@Q)
  J.has.2 <- S[cbind(2:n, 1:(n - 1))] != 0
  k <- sum(J.has.2)
  R.index <- vector("list", n - k)
  
  
  
  # optimized part 1 - start  ---------
  
  result <- part_1_plus_2(n, k, R.index, S, as.list(ev))
  R.index <- result[[1]]
  X <- result[[2]]
  
  # optimized part 1 - end  ---------
  
  I <- diag(2)
  
  if (n - k > 1) 
    for (j in 2:(n - k)) {
      ij <- R.index[[j]]
      for (i in (j - 1):1) {
        ii <- R.index[[i]]
        sumU <- 0
        if (length(ij) == 1 & length(ii) == 1) {
          if (j - i > 1) 
            for (l in (i + 1):(j - 1)) {
              il <- R.index[[l]]
              sumU <- sumU + {
                if (length(il) == 2) 
                  X[ii, il] %*% X[il, ij]
                else X[ii, il] * X[il, ij]
              }
            }
          X[ii, ij] <- solve(X[ii, ii] + X[ij, ij], S[ii, 
                                                      ij] - sumU)
        }
        else if (length(ij) == 2 & length(ii) == 1) {
          if (j - i > 1) 
            for (l in (i + 1):(j - 1)) {
              il <- R.index[[l]]
              sumU <- sumU + {
                if (length(il) == 2) 
                  X[ii, il] %*% X[il, ij]
                else X[ii, il] * X[il, ij]
              }
            }
          X[ii, ij] <- solve(t(X[ii, ii] * I + X[ij, 
                                                 ij]), as.vector(S[ii, ij] - sumU))
        }
        else if (length(ij) == 1 & length(ii) == 2) {
          if (j - i > 1) 
            for (l in (i + 1):(j - 1)) {
              il <- R.index[[l]]
              sumU <- sumU + {
                if (length(il) == 2) 
                  X[ii, il] %*% X[il, ij]
                else X[ii, il] * X[il, ij]
              }
            }
          X[ii, ij] <- solve(X[ii, ii] + X[ij, ij] * 
                               I, S[ii, ij] - sumU)
        }
        else if (length(ij) == 2 & length(ii) == 2) {
          if (j - i > 1) 
            for (l in (i + 1):(j - 1)) {
              il <- R.index[[l]]
              sumU <- sumU + {
                if (length(il) == 2) 
                  X[ii, il] %*% X[il, ij]
                else X[ii, il] %*% t(X[il, ij])
              }
            }
          tUii <- matrix(0, 4, 4)
          tUii[1:2, 1:2] <- X[ii, ii]
          tUii[3:4, 3:4] <- X[ii, ii]
          tUjj <- matrix(0, 4, 4)
          tUjj[1:2, 1:2] <- t(X[ij, ij])[1, 1] * I
          tUjj[3:4, 3:4] <- t(X[ij, ij])[2, 2] * I
          tUjj[1:2, 3:4] <- t(X[ij, ij])[1, 2] * I
          tUjj[3:4, 1:2] <- t(X[ij, ij])[2, 1] * I
          X[ii, ij] <- solve(tUii + tUjj, as.vector(S[ii, 
                                                      ij] - sumU))
        }
      }
    }
  
  
  Q %*% X %*% solve(Q)
}