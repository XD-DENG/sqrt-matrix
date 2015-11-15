# https://en.wikipedia.org/wiki/Square_root_of_a_matrix

Babylonian_sqrtm <- function(A){
  X <- diag(dim(A)[1])
  
  error <- 1
  error_tolerance <- 1.5e-8
  
  flag <- 1
  while(error>error_tolerance){

    X_old <- X
    X=(X + A %*% solve(X))/2
    error <- max(abs(X_old - X))
    
    flag <- flag + 1
  }
  
  return(X)
}