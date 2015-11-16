library(Matrix)


expm_sqrtm <- function (x) 
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

  
  # XD part-1 -------------------
  # XD: can be faster
  
  l <- 1
  i <- 1
  
  while (i < n) {
    if (S[i + 1, i] == 0) {
      R.index[[l]] <- i
    }
    else {
      R.index[[l]] <- (i:(i + 1))
      i <- i + 1
    }
    i <- i + 1
    l <- l + 1
  }

  # XD part-1 -------------------
  
  if (is.null(R.index[[n - k]])) {
    R.index[[n - k]] <- n
  }
  
  
  # XD part-2 -------------------
  # XD: can be faster
  I <- diag(2)
  X <- matrix(0, n, n)
  
  
  for (j in seq_len(n - k)) {
    ij <- R.index[[j]]
    if (length(ij) == 1) {
      X[ij, ij] <- if ((.s <- S[ij, ij]) < 0) 
        sqrt(.s + (0+0i))
      else sqrt(.s)
    }
    else {
      ev1 <- Sch.x@EValues[ij[1]]
      r1 <- Re(sqrt(ev1))
      X[ij, ij] <- r1 * I + 1/(2 * r1) * (S[ij, ij] - Re(ev1) * I)
    }
  }
  # XD part-2 -------------------

  
  # XD part - 3 begin --------------
  
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
  # XD part - 3 end --------------
  
  Q %*% X %*% solve(Q)
}