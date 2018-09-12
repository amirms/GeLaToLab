## Simple kernel canonical corelation analysis
## author: alexandros karatzoglou

setGeneric("kcca",function(Kx, Ky, kpar=list(sigma = 0.1), gamma=0.1, ncomps = 10, ...) standardGeneric("kcca"))
setMethod("kcca", signature(Kx = "matrix"),
          function(Kx, Ky, kpar=list(sigma=0.1), gamma=0.1, ncomps =10, ...)
          {
            Kx <- as.matrix(Kx)
            Ky <- as.matrix(Ky)
            
            if(!(nrow(Kx)==nrow(Ky)))
              stop("Number of rows in Kx, Ky matrixes is not equal")

            n <- dim(Kx)[1]
            m <- 2
            ## Generate LH
            VK <- matrix(0,n*2,n);
            
            VK[0:n,] <- Kx
            VK[(n+1):(2*n),] <- Ky
            LH <- tcrossprod(VK, VK)
            
            for (i in 1:m)
              LH[((i-1)*n+1):(i*n),((i-1)*n+1):(i*n)] <- 0
            
            ## Generate RH
            RH <- matrix(0,n*m,n*m)
            RH[1:n,1:n] <- (Kx + diag(rep(gamma,n)))%*%Kx + diag(rep(1e-6,n))
            RH[(n+1):(2*n),(n+1):(2*n)] <- (Ky + diag(rep(gamma,n)))%*%Ky + diag(rep((1e-6),n))
            RH <- (RH+t(RH))/2
            
            ei <- .gevd(LH,RH)
            
            ret <- list()
            
            kcor <- as.double(ei$gvalues[1:ncomps])
            xcoef <- matrix(as.double(ei$gvectors[1:n,1:ncomps]),n)
            ycoef <- matrix(as.double(ei$gvectors[(n+1):(2*n),1:ncomps]),n)
            ## xvar(ret) <- rotated(xpca) %*% cca$xcoef
            ## yvar(ret) <- rotated(ypca) %*% cca$ycoef
            return(list(kcor=kcor, xcoef=xcoef, ycoef=ycoef))
          })

## gevd compute the generalized eigenvalue 
## decomposition for (a,b)
.gevd<-function(a,b=diag(nrow(a))) {
  bs<-.mfunc(b,function(x) .ginvx(sqrt(x)))
  ev<-eigen(bs%*%a%*%bs)
  return(list(gvalues=ev$values,gvectors=bs%*%ev$vectors))
}

## mfunc is a helper to compute matrix functions
.mfunc<-function(a,fn=sqrt) {
  e<-eigen(a); y<-e$vectors; v<-e$values
  return(tcrossprod(y%*%diag(fn(v)),y))
}

center.kernel <- function(K) {
  n <- nrow(K)
  return(K - 2 * diag(n) %*% K + diag(n) %*% K %*% diag(n))
}

## ginvx is a helper to compute reciprocals
.ginvx<-function(x) {ifelse(x==0,0,1/x)}
