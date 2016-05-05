
cov = function(m1, m2) {
  #check for compatibility of matrices
  return(t(m1) %*% m2)
  
}


MLE <- function(iters){

m1 = matrix(1, nrow=3, ncol=4)
m2 = matrix(2, nrow=3, ncol=3)

colnames(m1)= c("A", "B", "C", "D")
m1[1,1] = 8
m1[1,3] = 5
m1[1,4] = 4
m1[2,2] = 3
m1[2,3] = 11
m1[3,1] = 2
m1[3,2] = 4

colnames(m2) = c("X", "Y", "Z")
m2[1,2] = 5
m2[2,1] = 1
m2[2,2] = 3
m2[3,1] = 9
m2[3,2] = 5
m2[3,3] = 6

V1 = cov(m1, m1)
V2 = cov(m2, m2)

V12 = cov(m1, m2)

iters = 0
x = svd(V12)

#repeat {
  

  E1 = V1 - x$u %*% t(x$u)
  e1 = eigen(E1)
  
  
  E2 = V2 - x$v %*% t(x$v)
  e2 = eigen(E2)
  
  
  if (iters >= 100)
    break
#}

}