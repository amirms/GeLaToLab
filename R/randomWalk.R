"
    Random walks graph kernel.
"

compare <- function(g1, g2, alpha, verbose=FALSE){
  "Compute the kernel value (similarity) between two graphs. 

Parameters
----------
g1 : networkx.Graph
First graph.
g2 : networkx.Graph
Second graph.
alpha : interger < 1
A rule of thumb for setting it is to take the largest power of 10
which is samller than 1/d^2, being d the largest degree in the 
dataset of graphs.    

Returns
-------        
k : The similarity value between g1 and g2.
"
require(pcg)
# am1 = nx.adj_matrix(g1)
# am2 = nx.adj_matrix(g2)
# x = np.zeros((len(am1),len(am2)))
x = matrix(0, dim(am1)[1],dim(am2)[1])
A = smt_filter(x,am1,am2,alpha)
# b = np.ones(dim(am1)[1]*dim(am2)[1])
newDim <- dim(am1)[1]*dim(am2)[1]
b = rep(1, newDim)
tol = 1e-6
maxit = 20
pcg(A,b,x,maxit,tol)
return(sum(x))
}


compare.cgm <- function(Ax, maxit=20){
  require(pcg)
  
  b = rep(1, dim(Ax)[1])
  
  tol = 1e-6
#   pcg.2(A,b, maxiter=maxit, tol=tol)
pcg(A,b, maxiter=maxit, tol=tol)

  return(sum(x))
}


pcg.2 <- function (A, b, M, maxiter = 1e+05, tol = 1e-06) 
{
  if (missing(M)) {
    dA <- diag(A)
    dA[which(dA == 0)] = 1e-04
    Minv = diag(1/dA, nrow = nrow(A))
  }
  else Minv = solve(M)
  x = rep(0, length(b))
  r = b - A %*% x
  z = Minv %*% (r)
  p = z
  iter = 0
  sumr2 = sum(r^2)
  while (sumr2 > tol & iter < maxiter) {
    iter = iter + 1
    Ap = A %*% p
    a = as.numeric((t(r) %*% z)/(t(p) %*% Ap))
    x = x + a * p
    r1 = r - a * Ap
    z1 = Minv %*% r1
    bet = as.numeric((t(z1) %*% r1)/(t(z) %*% r))
    p = z1 + bet * p
    z = z1
    r = r1
    sumr2 = sum(r^2)
  }
  if (iter >= maxiter) 
    print("pcg did not converge. You may increase maxiter number.")
  return(x)
}


# compute.random.walk.sim <- function(Ax, verbose=FALSE){
#   "Compute the kernel value (similarity) between two graphs. 
# 
# Parameters
# ----------
# Ax
# 
# Returns
# -------        
# k : The similarity value between g1 and g2.
# "
#   require(pcg)
#   # b = np.ones(dim(am1)[1]*dim(am2)[1])
#   newDim <- dim(am1)[1]*dim(am2)[1]
#   b = matrix(1, newDim, newDim)
#   tol = 1e-6
#   maxit = 20
#   pcg(Ax,b,x,maxit,tol)
#   return(sum(x))
# }


smt_filter <- function(x, am1, am2, alpha){
#   yy = np.dot(np.dot(am1, x), am2)
  yy = am1 %*% x %*% am2
yy = yy * alpha
vecu = x - yy
return(vecu)
}

#lex.fun <- compute_normalized_LCS
compare_rw_graphs <- function(a1, a2, alpha, lex.fun=NULL, min.nchar=5) {
  require(MASS)
  require(Matrix)
  #Make sure the row and column names of each adjacency matrix are the same
  # stopifnot(rownames(a1) == colnames(a1))
  # stopifnot(rownames(a2) == colnames(a2))
  
  vertices1 <- rownames(a1)
  vertices2 <- rownames(a2)
  
  lex.v <- lapply(vertices1, function(v1) lapply(vertices2, function(v2) lex.fun(v1, v2) ))
  lex.v <- do.call(rbind, lex.v)
  lex.v <- matrix(as.numeric(unlist(lex.v)),nrow=nrow(lex.v), dimnames=list(vertices1, vertices2))
  # rownames(lex.v) <- vertices1
  # colnames(lex.v) <- vertices2

#   lex.v <- lex.v[which(!apply(lex.v,1,FUN = function(x){all(x == 0)})),which(!apply(lex.v,2,FUN = function(x){all(x == 0)}))] 
  
  
  Vx = kronecker(lex.v, t(lex.v), make.dimnames=T)
  
#   a1 <- normalize_adjacency_matrix(a1)
#   a2 <- normalize_adjacency_matrix(a2)

  Gx = kronecker(a1, a2, make.dimnames=T)
# stopifnot(rownames(Gx) == rownames(Gx))

  columnNames <- unlist(lapply(colnames(Gx), function(name) {
    splitted <- strsplit(name, ":")[[1]]
    return(paste(splitted[2], splitted[1], sep=":"))
  }))

#NOTE rownames(Gx) is equal to rownames(Vx)
# Vx <- Vx[rownames(Gx), columnNames]
  Vx <- Vx[, columnNames]
# Gx <- Matrix(Gx, sparse=T)
  
  Ax <- Vx * Gx
  rm(Gx)
  rm(Vx)
  
  k <- dim(Ax)[1]
  
  gc()
  
  Z <- Matrix(diag(k) - alpha * Ax)
  
  inverseAx <- solve(Z)
  
  e <- rep(1/k, k)
  
 return(t(e) %*% inverseAx %*% e)
}


# normalize by outdegree
normalize_adjacency_matrix <- function(adj){
  outdegrees <- rowSums(adj)
  
  m <- adj / outdegrees
  
  m[is.nan(m)] <- 0
  
  m
}

compare_normalized <- function(self, g1, g2, alpha, verbose=False){
  "Compute the normalized kernel value between two graphs. 

A normalized version of the kernel is given by the equation: 
k_norm(g1, g2) = k(g1, g2) / sqrt(k(g1,g1) * k(g2,g2)) 

Parameters
----------
g1 : networkx.Graph
First graph.
g2 : networkx.Graph
Second graph.
alpha : interger < 1
A rule of thumb for setting it is to take the largest power of 10
which is samller than 1/d^2, being d the largest degree in the 
dataset of graphs.

Returns
-------        
k : The similarity value between g1 and g2.
"
}



compare_list <- function(graph_list, alpha, verbose=FALSE){
  "Compute the all-pairs kernel values for a list of graphs. 

This function can be used to directly compute the kernel matrix 
for a list of graphs. The direct computation of the kernel matrix is
faster than the computation of all individual pairwise kernel values.        

Parameters
----------
graph_list: list
A list of graphs (list of networkx graphs)
alpha : interger < 1
A rule of thumb for setting it is to take the largest power of 10
which is samller than 1/d^2, being d the largest degree in the 
dataset of graphs.

Return
------
K: numpy.array, shape = (len(graph_list), len(graph_list))
The similarity matrix of all graphs in graph_list.
"

n = len(graph_list)
# K = np.zeros((n,n))
for (i in range(n-1))
  for (j in range(i+1, n))
  K[i,j] = self.compare(graph_list[i], graph_list[j], alpha)

K[j,i] = K[i,j]
return(K)
}



compare_list_normalized <- function(self, graph_list, alpha, verbose=FALSE){      
  "Compute the all-pairs kernel values for a list of graphs. 
        
        A normalized version of the kernel is given by the equation: 
        k_norm(g1, g2) = k(g1, g2) / sqrt(k(g1,g1) * k(g2,g2)) 
                
        Parameters
        ----------
        graph_list: list
            A list of graphs (list of networkx graphs)
        alpha : interger < 1
            A rule of thumb for setting it is to take the largest power of 10
            which is samller than 1/d^2, being d the largest degree in the 
            dataset of graphs.     
            
        Return
        ------
        K: numpy.array, shape = (len(graph_list), len(graph_list))
        The similarity matrix of all graphs in graph_list.
        "       

k = self.compare_list(graph_list, alpha)
k_norm = np.zeros(k.shape)
for (i in range(k.shape[0]))
  for (j in range(k.shape[1]))
  k_norm[i,j] = k[i,j] / np.sqrt(k[i,i] * k[j,j])

return(k_norm)
}