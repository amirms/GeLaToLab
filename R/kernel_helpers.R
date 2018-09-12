squared.euclidean.distance.of.kernel.matrix <- function(K){
  library(matrixcalc)
  
  # stopifnot(is.positive.semi.definite(K, tol=1e-8))
  
  e <- rep(1, dim(K)[1])
  
  return( (diag(K) %*% t(e)) + e %*% t(diag(K)) -2*K )
}