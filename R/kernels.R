
spectrum.string.kernel = function(txts, len, lam=1) {
  require(kernlab)
  
  sk <- stringdot(type="spectrum", length=len, #lambda = lam,
                  normalized=TRUE)
  
  k <- kernelMatrix(sk, txts)
  
  dimnames(k) <- list(names(txts), names(txts))
  
  return(k)
  
}

exponential.decay.kernel  = function(txts, lam) {
  require(kernlab)
  
  sk <- stringdot(type="exponential", lambda= lam, normalized=TRUE)
  
  k <- kernelMatrix(sk, txts)
  
  k[!is.finite(k)] <- 0
  
  dimnames(k) <- list(names(txts), names(txts))
  
  return(k)
  
}

bound.range.string.kernel = function(txts, len, lam=1) {
  require(kernlab)
  
  sk <- stringdot(type="boundrange", length=len, #lambda=lam, 
                  normalized=TRUE)
  
  k <- kernelMatrix(sk, txts)
  
  dimnames(k) <- list(names(txts), names(txts))
  
  return(k)
  
}

constant.string.kernel = function(txts) {
  require(kernlab)
  
  lower_txts <- tolower(txts)
  
  sk <- stringdot(type="constant", normalized=TRUE)
  
  k <- kernelMatrix(sk, lower_txts)
  
  dimnames(k) <- list(names(txts), names(txts))
  
  return(k)
  
}


gaussian.kernel = function(bow, sigma) {
  require(kernlab)
  
  rbf = rbfdot(sigma)
  return(kernelMatrix(rbf, bow))
  
}

linear.kernel = function(bow) {
  require(kernlab)
  
  ply = polydot(degree = 1, scale = 1, offset = 0)
  
  return(kernelMatrix(ply, bow))
  
}

polynomial.kernel = function(bow, degree, scale=1, offset=0) {
  require(kernlab)
  
  ply = polydot(degree = degree, scale = scale, offset = offset)
  
  return(kernelMatrix(ply, bow))
  
}

calc.diffusion.kernel <- function(L, beta=0, is.adjacency=F){
  require(igraph)
  if(missing(L)) stop('You have to provide either a graph laplacian to calculate the diffusion kernel.')
  method="thresh"
    if(is.adjacency){
      dnames <- dimnames(L)
      L = graph.adjacency(L, mode='undirected', diag=FALSE, weighted=T)
      L = graph.laplacian(L, normalized=TRUE)
      dimnames(L) = dnames
    }
  
  if(method == "thresh"){
    eig = eigen(L)
    ## Let's discard 80% of the eigenvector with largest eigenvalues
    ncomp = round(0.8*ncol(L))+1
    V = eig$vectors[,ncol(L):ncomp]
    ## diff kernel
    R = V %*% diag(exp(-beta*eig$values[ncol(L):ncomp])) %*% t(V)
  }
  ## erstmal nur mehtod = 'thres zulassen'
  ##  else if(method == "pseudoinv")
  ##    R = pseudoinverse(L)
  ##  else if(method == "LLE"){  
  ##    if(!is.null(G)){
  ##      G = ugraph(G)
  ##      W = as(G, "matrix")
  ##      deg = rowSums(W)
  ##      Dinv = diag(1/deg)
  ##      Dinv[Dinv == Inf] = 1  
  ##      P = Dinv%*%W
  ##    }		
  ##    else{
  ##      P = L
  ##    }					
  ##    I = diag(ncol(P)) 
  ##    T = I - P
  ##    M = t(T)%*%(T)	
  ##    lam = eigen(M, symmetric=TRUE, only.values=TRUE)
  ##    K = lam$values[1]*I - M
  ##    E = I - matrix(1,ncol=ncol(I), nrow=nrow(I))/ncol(I)
  ##    R = E%*%K%*%E		
  ##  }
  dimnames(R) = dimnames(L)
  ## Normalize Kernel
  KernelDiag <- sqrt(diag(R) + 1e-10)
  R.norm <- R/(KernelDiag %*% t(KernelDiag))
  R.norm
}


Diffusionfct <- function(
  adj,
  beta=0,
  correct.neg = 1e-10
) {
  # calculates the diffusion-kernel based distance matrix for an adjacency matrix
  
  # if the adjacency matrix contains only a single gene, return a 1x1 matrix containing 0
  if (is.null(dim(adj))) return(matrix(0, 1, 1))
  
  # compute the diffusion kernal
  H <- adj - Diagonal(x=apply(adj, 1, sum))
  x <- eigen(H, symmetric=T)
  K  <- x$vectors %*% diag(exp(beta * x$values)) %*% t(x$vectors)
  Dsub <- outer(diag(K), diag(K), "+") - 2*K
  
  # correct negative distances
  if (any(Dsub < 0) && correct.neg){
    warning("negative distances set to zero")
    Dsub[Dsub < 0] <- 0
  }
  
  sqrt(Dsub)
}


####################################################################################
##
## Distance measures on graphs
##

build.graph <- function(adj, weighted=T){
  require(igraph)
  adj <- make.symmetric(adj)
  graph.adjacency(adj, mode='undirected', diag=FALSE, weighted=weighted)
}

## produces generator matrix for diffusion kernel
## if all edge weights = 1, this is the same as -graph.laplacian(g)
generator <- function(g){  
  if(is.null(E(g)$weight)) E(g)$weight = 1
  A <- get.adjacency(g,attr="weight")
  D <- diag(apply(A,1,sum))
  H <- A-D
  return(H)
} ## end function generator


## Diffusion Kernels on Graphs and Other Discrete Structures, Kondor and Lafferty, 2002
graph.diffusion <- function(g,beta=1,v=V(g),correctNeg=TRUE){
  if (class(v)!="igraph.vs") v <- my.as.igraph.vs(v-1)
  H <- generator(g)
  x <- eigen(H) ## H ~ x$vectors %*% diag(x$values) %*% t(x$vectors)
  K <- x$vectors %*% diag(exp(beta * x$values)) %*% t(x$vectors)
  D2 <- outer(diag(K),diag(K),"+") - 2*K
  if (any(D2<0) && correctNeg){
    warning("Negative 'distances' set to Zero!")
    D2[D2<0] <- 0
  }
  D <- sqrt(D2)
#   K <- K[v+1,v+1]
#   D <- D[v+1,v+1]
 #   dimnames(D) <- list(v,v)
 # dimnames(K) <- list(v,v)
  return(list(kernel=K,dist=D))
} 