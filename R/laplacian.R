#This is borrowed from implicit-bayes

laplacian <- function(adj, normalized = FALSE, ...) {
  #d <- degrees(adj)
  #   adj <- make.symmetric(cfg)
  diag(adj) <- 0 
  d <- apply(abs(adj),1,sum)
  lap <- -adj
  diag(lap) <- diag(lap) + d
  if (normalized) {
    di <- 1/sqrt(d)
    di[is.infinite(di)] <- 0
    lap <- diag(di) %*% lap %*% diag(di)
  }
  dimnames(lap) <- dimnames(adj)
  lap
  
}

pseudoinv <- function(a, tol = 1e-8, scale = FALSE) {
  udv <- svd(a)
  active <- udv$d > tol
  ainv <- udv$u[,active] %*% ((1/udv$d[active]) * t(udv$v[,active]))
  if (scale)
    ainv <- ainv / sum(diag(ainv))
  ainv
}

# RandomWalkFunctions --------------------------------------------------------


get.pseudo.inverse.random.walk  <- function(A, D){
  # Calculates pseudo inverse of the laplacian matrix (Eq. 5.1)
  # 
  # Args:
  # A: adjacency matrix of the network
  # D: diagonal matrix of the network
  #
  # Returns:
  # The pseudo inverse of the Laplacian matrix
  ones <- matrix(100, nrow = nrow(A), ncol = 1)
  some.factor <- as.numeric(t(ones) %*% ones / nrow(A))  # See equation 5.1
  print(some.factor)
  
  PseudoL <- solve(((D - A) - some.factor), ) + some.factor
  return(PseudoL)
}


get.prob.of.absorp  <- function(L, initial.node,  halt.node, end.node){
  # Calculate probability of absorption by node a
  # before absorption by node b, starting at node k (Eq. 5.2)
  # 
  # Args:
  # L: psuedo inverse of Laplacian matrix, obtained from GetPseudoInverse
  # end.node: 1st absorbing state, the desired end-point
  # halt.node: 2nd absorbing state the undesired or halt end-point
  # initial.node: initial node at which you start the random walk
  #
  # Returns: 
  # The probability of absorption by state a before absorption by state b,
  # when starting f rom state k. 
  ProbOfAbsorp <- (L[initial.node, end.node] - L[initial.node, halt.node] - 
                     L[end.node, halt.node] + L[halt.node, halt.node]) / 
    (L[end.node, end.node] + L[halt.node, halt.node] - 
       2 * L[end.node, halt.node])
  if(abs(ProbOfAbsorp) < 10^-9){
    #If ProbOfAbsorp is vanishingly small (or negative), round to zero.
    ProbOfAbsorp <- 0
  }
  return(ProbOfAbsorp)
}


get.avg.first.passage.time <- function(D, L, initial.node, end.node){
  # Calculate average first passage time (Eq. 5.4)
  # Defined as the average number of steps that a random walker,
  # starting in node i, will take to enter node k for the first time.
  #
  # Args:
  # D: diagonal matrix used to calculate pseudoinverse
  # L: psuedo inverse of Laplacian matrix, obtained from GetPseudoInverse
  # initial.node: initial starting node at which you start random walk
  # end.node: destination node at which you stop your random walk
  #
  # Returns:
  # AvgFirstPassageTime: the average number of steps that a random walker,
  # starting in node i, will take to enter node k for the first time.
  AvgFirstPassageTime <- 0  #Initialize before putting into loop
  for(j in 1:nrow(L)){
    AvgFirstPassageTime <- AvgFirstPassageTime + 
      (L[initial.node, j] - L[initial.node, end.node] - 
         L[end.node, j] + L[end.node, end.node]) * D[j, j]   
  }
  if(abs(AvgFirstPassageTime) < 10^-9){
    #Round down to zero if steps are less than 1
    AvgFirstPassageTime = 0
  }
  return(AvgFirstPassageTime)
}


get.avg.commute.time <- function(D, L, initial.node, end.node){
  # Calculate average commute time (Eq 5.6)
  # 
  # Args:
  # D: diagonal matrix used to calculate pseudoinverse
  # L: psuedo inverse of Laplacian matrix, obtained from GetPseudoInverse
  # initial.node: initial starting node
  # end.node: destination node
  # 
  # Returns:
  # AvgCommuteTime: defined as average number of steps that a random walker,
  # starting in node i, will take to enter node k for the first time and
  # go back to node i.
  AvgCommuteTime <- sum(diag(D)) * (L[initial.node, initial.node] + 
                                      L[end.node, end.node] - 2 * L[initial.node, end.node])
  if(abs(AvgCommuteTime) < 1){
    #Round down to zero if steps are less than 1
    AvgCommuteTime = 0
  }
  return(AvgCommuteTime)
}


# DiffusionModelFunctions -------------------------------------------------

solve.shifted.laplacian <- function(A, D, gamma){
  I <- diag(nrow(A))
  L <- -(A - D - gamma*I)
  L_inv_shifted <- solve(L)
  return(L_inv_shifted)
}


get.steady.state.diffusion <- function(L_inv_shifted, remove_diags=TRUE){
  length_of_matrix <- nrow(L_inv_shifted)
  influence_matrix <- matrix(0, ncol = length_of_matrix, nrow = length_of_matrix)
  for(i in 1:nrow(L_inv_shifted)){
    b <- matrix(0, ncol=1, nrow=length_of_matrix, byrow=TRUE)
    b[i] <- 1
    influence_matrix[i, ] <- as.vector(L_inv_shifted %*% b)
  }
  if (remove_diags==TRUE){
    diag(influence_matrix) <- 0
  }
  return(influence_matrix)
}

#compare graph kernels

compare.graph.kernels <- function(A){
  # ObtainMatrixInfo --------------------------------------------------------
  
  #   A <- get.adjacency(g)    # Adjacency matrix
  L <- graph.laplacian(A)    # Laplacian matrix
  D <- diag(diag(L))    # Diagonal matrix
}


# DiffusionKernel ---------------------------------------------------------
compute.diffusion.kernel <- function(A, D, gamma) { 
  #   gamma <- 1
  L_inv_shifted <- solve.shifted.laplacian(A, D, gamma)
  influence_matrix_diffusion <- get.steady.state.diffusion(L_inv_shifted, 
                                                           remove_diags=TRUE)
  
  influence_matrix_diffusion
}

# RandomWalkModel ---------------------------------------------------------

compute.avg.commute.time.kernel <- function(A,D) {
  L_inv_pseudo <- get.pseudo.inverse.random.walk(A, D)
  
  influence_matrix_randwalk <- diag(0, nrow=nrow(L_inv_pseudo))    # Initialize
  
  for(i in 1:nrow(L_inv_pseudo)){
    for(j in 1:nrow(L_inv_pseudo)){
      #     influence_matrix_randwalk[i, j] <- GetAvgFirstPassageTime(D,
      #                                                               L_inv_pseudo,
      #                                                               i,
      #                                                               j)
      influence_matrix_randwalk[i, j] <- get.avg.commute.time(D,
                                                              L_inv_pseudo,
                                                              i,
                                                              j)
    }
    print(c('All random walk times computed for node ', i))
  }
  
  dimnames(influence_matrix_randwalk) <- dimnames(A)
  influence_matrix_randwalk
}



compute.sigmoid.commute.time.kernel <- function(A,D, p) {
  L_inv_pseudo <- get.pseudo.inverse.random.walk(A, D)
  
  sigma <- sd(L_inv_pseudo)
  
  influence_matrix_randwalk <- diag(0, nrow=nrow(L_inv_pseudo))    # Initialize
  
  for(i in 1:nrow(L_inv_pseudo)){
    for(j in 1:nrow(L_inv_pseudo)){
      #     influence_matrix_randwalk[i, j] <- GetAvgFirstPassageTime(D,
      #                                                               L_inv_pseudo,
      #                                                               i,
      #                                                               j)
      
      
      sigmoid.commute.time <- exp((-p * L_inv_pseudo[i,j])/sigma)
      
      influence_matrix_randwalk[i, j] <- sigmoid.commute.time
    }
    print(c('All random walk times computed for node ', i))
  }
  
  dimnames(influence_matrix_randwalk) <- dimnames(A)
  influence_matrix_randwalk
}

compute.avg.first.passage.time.kernel <- function(A, D) {
  L_inv_pseudo <- get.pseudo.inverse.random.walk(A, D)
  
  influence_matrix_randwalk <- diag(0, nrow=nrow(L_inv_pseudo))    # Initialize
  
  for(i in 1:nrow(L_inv_pseudo)){
    for(j in 1:nrow(L_inv_pseudo)){
      influence_matrix_randwalk[i, j] <- get.avg.first.passage.time(D,
                                                                    L_inv_pseudo,
                                                                    i,
                                                                    j)
      #       influence_matrix_randwalk[i, j] <- get.avg.commute.time(D,
      #                                                               L_inv_pseudo,
      #                                                               i,
      #                                                               j)
    }
    print(c('All random walk times computed for node ', i))
  }
  
  influence_matrix_randwalk
}


compute.prob.absorption.kernel <- function(A, D) {
  L_inv_pseudo <- get.pseudo.inverse.random.walk(A, D)
  
  influence_matrix_randwalk <- diag(0, nrow=nrow(L_inv_pseudo))    # Initialize
  
  for(i in 1:nrow(L_inv_pseudo)){
    for(j in 1:nrow(L_inv_pseudo)){
      influence_matrix_randwalk[i, j] <- get.avg.first.passage.time(D,
                                                                    L_inv_pseudo,
                                                                    i,
                                                                    j)
      #       influence_matrix_randwalk[i, j] <- get.avg.commute.time(D,
      #                                                               L_inv_pseudo,
      #                                                               i,
      #                                                               j)
    }
    print(c('All random walk times computed for node ', i))
  }
  
  influence_matrix_randwalk
}

# (D - \alpha A)^{-1}D
compute.random.walk.with.restart <- function(A, D, alpha) {
  
  
  L_inv_pseudo <- get.pseudo.inverse.random.walk(A, D)
  
  L_inv_pseudo %*% D
}


#Normalization function

normalize.influence <- function(influence_matrix){
  # Normalizes influence by dividing by inverse of 'total influence'
  total_influence <- apply(influence_matrix, 1, function(x){
    1 / sum(x)
  })
  influence_matrix.norm <- matrix(NA, nrow=nrow(influence_matrix), 
                                  ncol=ncol(influence_matrix))
  for (i in 1:nrow(influence_matrix)){
    for (j in 1:ncol(influence_matrix)){
      influence_matrix.norm[i, j] <- influence_matrix[i, j] / total_influence[i]
    }
  }
  return(influence_matrix.norm)
}


# Normalize adjacency matrix by node degree.

normalize_adjacency_by_degree <- function(A, D){
  for (i in 1:nrow(A)){
    for (j in 1:ncol(A)){
      if (A[i, j]==1){
        A[i, j] = A[i, j] / diag(D)[j]
      }
    }
  }
  return(A)
}

#A corresponds to a directed graph

# compute.directed.avg.commute.time.kernel <- function(A, D) {
# 
#   #transition-probability matrix
#   P <-  solve(D) %*% A 
#     
#   eL <- eigen(t(P))    
#   eL1 <- eL$vectors[,1]
#   
#   Pi <- diag(eL1)
#   
#   Pi - (Pi%*%P + t(P) %*% Pi)/2
# }


compute.von.neumann.diffusion.kernel <- function(A, alpha) {
  K <- solve(diag(dim(A)[1]) - alpha * A)
  
  dimnames(K) <- dimnames(A)
  
  K
  
}


compute.exponential.diffusion.kernel <- function(A, alpha) {
  require(Matrix)
  K <- as.matrix(expm(alpha * A))
  
  
  dimnames(K) <- dimnames(A)
  
  K
  
}