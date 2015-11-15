DB_sqrtm <- function(A){
  
  Y <- A
  Z <- diag(dim(A)[1])
  
  error <- 1
  error_limit <- 1.5e-8
  i=1
  while(error>error_limit){
    Y_old <- Y
    Y <- (Y_old+solve(Z))/2
    Z <- (Z+solve(Y_old))/2
    error <- max(abs(Y-Y_old))
    i=i+1
  }
  
  return(Y)
  
}