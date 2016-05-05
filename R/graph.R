average.degree <- function(A) {
  
  n <- dim(A)[1]
  
  As <- make.symmetric(A)
  
  D = rowSums(A)
  
  
  return(sum(D)/n)
  
}