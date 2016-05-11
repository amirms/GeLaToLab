avg_lap <- function(Laps) {
  apply(simplify2array(Laps), 1:2, mean)
}

norm_vec <- function(x) sqrt(sum(x^2))

compute_all_cuts <- function(L1, L2){
  v = geigen(L1, L2) #VERIFIED
 
  #Remove the eigenvectors associated with eigenvalues zero and infinity
  vals <- v$values
  # #     
  zerov <- min(abs(vals))
  infv <- max(vals)
  
  vec <- v$vectors
  
  vec <- vec[,-c(which(vals==zerov), which(vals==(-zerov)), which(vals==infv))]  
  
  norm_vec <- function(x) sqrt(sum(x^2))
  
  #normalize each vector - VERIFIED
  for (i in 1:dim(vec)[2])
    vec[,i] = vec[,i] / norm_vec(vec[,i])
  
  vec
}

#Input: a cut with N number of cuts
is_dominated <- function(cut, included_cuts) {
  
  for(i in 1:length(included_cuts)){
   
    stopifnot(length(cut)==length(included_cuts[[i]]$costs))
    
    if (all(unlist(included_cuts[[i]]$costs) < unlist(cut)))
      return(TRUE)
    
  }
    
  return(FALSE)
}

#return Indices
dominates <- function(cut, included_cuts) {  
  
  indices <- c()
  
  for(i in 1:length(included_cuts)){
    
    stopifnot(length(cut)==length(included_cuts[[i]]$costs))
    
    if (all(unlist(included_cuts[[i]]$costs) >= unlist(cut)))
      indices <- c(indices, i)
    
  }
  
  indices    
}


#Input: a list of costs costs

compute_pareto_front <- function(costs, I, vec) {
  #TODO all costs should have the same length
  N <- length(costs[[1]])
  
  ex = list() # The cuts to be excluded
  U = list() # The output
  iter = 0
  pick = 0
  
  included_cuts <- list()
  
  initial_costs <- lapply(costs, function(cost) cost[1])
  included_cuts[[length(included_cuts) + 1]] <- list(costs = initial_costs, U=vec[,1])
  
  for (i in 2 : N){
    
    current_costs <- lapply(costs, function(cost) cost[i])
    
    if(!is_dominated(current_costs, included_cuts)) {
      indices <- dominates(current_costs, included_cuts)
      
      if (!is.null(indices))
        included_cuts <- included_cuts[-indices]
      
      included_cuts[[length(included_cuts) + 1]] <- list(costs = current_costs, U=vec[,i])
      
    }
    
  }

  included_cuts
  
}


#computes the Pareto front for multi view learning
# Input: a list of graph laplacians for each view, i.e. Ls

#compute the normalized laplacians on the kernels

compute_generalized_pareto_multiview <- function(Ls, combine_pareto_front=T) {
  
  require(geigen)
  require(ggplot2)
  
  #TODO Stopifnot(length(Ls) >= 2)
  
  #TODO check the dimensions match  
  N = dim(Ls[[1]])[1]
  
  vecs = list()
  
  
  #TODO check to see if this makes any sense
  #ALso column bind the resulting vectors
  vec =c()
#   print("reached here")
  for(i in 1:length(Ls)){
    temp <- compute_all_cuts(Ls[[i]], avg_lap(Ls[-i]))
    vec <- cbind(vec, temp)
    vecs[[i]] <- temp #Expensive operation
  }
#   print("reached here2")
  costs <- lapply(Ls, function(L) diag(t(vec) %*% L %*% vec))
  
  #___
#   require(scatterplot3d)
#   s3d <- scatterplot3d(cost1,cost2,cost3, main="Pareto Frontier", color="blue",  pch=3)
  #___
  
  #order the costs based on the first cost
  I <- order(costs[[1]])
  costs[[1]] <- sort(costs[[1]])
  
  for (i in 2:length(costs))
    costs[[i]] <- costs[[i]][I]
  
  pareto_front <- compute_pareto_front(costs, I, vec)

  print("number of cuts")
  print(length(pareto_front))

  pareto_costs <- unlist(lapply(pareto_front, function(p) p$costs))

  #______
#   s3d$points3d(pareto_cost1, pareto_cost2, pareto_cost3, cex=1.5, pch=16, col = "dark red")
  #______

  pareto_U <- lapply(pareto_front, function(p) p$U)

  if (!combine_pareto_front)
    return(pareto_U)

  pareto_U <- matrix(unlist(pareto_U), ncol = length(pareto_U), byrow = TRUE)

  rownames(pareto_U) <- rownames(Ls[[1]])

  return(pareto_U)

}

compute.igraph.laplacian <- function(A, norm=TRUE) {
  require(igraph)
  
  dnames <- dimnames(A)  
  
  if (!isSymmetric(A))
    g <- graph.adjacency(A, mode="undirected", weighted=TRUE, diag=FALSE) 
  else
    g <- graph.adjacency(A, mode="directed", weighted=TRUE, diag=FALSE) 
  
  L <- as.matrix(graph.laplacian(g, norm, weights = NULL))
  dimnames(L) = dnames
  
  L
}


compute.laplacian <- function(A) {
  
  if (!(isSymmetric(A)))
    A <- make.symmetric(A)
  
  N <- dim(A)[2]
  
  # Compute the graph Laplacian.
  D = diag(colSums(A))
  
#   vol = sum(diag(D))
  
  D_norm = solve(sqrt(D))
  L = diag(N) - D_norm%*%A%*%D_norm
  
  dimnames(L) <- dimnames(A)  
  
  
  return(L)
}

column.normalize <- function(vec) {
  norm_vec <- function(x) sqrt(sum(x^2))
  
  #normalize each vector
  for (i in 1:dim(vec)[2])
    vec[,i] = vec[,i] / norm_vec(vec[,i])
  
  vec
}

# 
# 
# pareto_multiview <- function(L1, L2) {
#   
#   require(geigen)
#   require(ggplot2)
#   
#   t1<-theme(                              
#     plot.background = element_blank(), 
#     panel.grid.major = element_blank(), 
#     panel.grid.minor = element_blank(), 
#     panel.border = element_blank(), 
#     panel.background = element_blank(),
#     axis.line = element_line(size=.4),
#     axis.text=element_text(size=12),
#     axis.title=element_text(size=14,face="bold"),
#     legend.title=element_text(size=10, vjust=-12)
#     #     guide_colourbar.title = element_text(draw.ulim = FALSE, draw.llim = FALSE)
#     #    legend.key.height=unit(3,"line"),
#     #    legend.key.width=unit(3,"line")
#   )
#   
#   #TODO check the dimensions match
#   
#   N = dim(L1)[1]
#   
#   norm_vec <- function(x) sqrt(sum(x^2))
#   
#   v = geigen(L1,L2) #VERIFIED
#   #   v = eigen(solve(L2)%*%L1) 
#   
#   
#   #Remove the eigenvectors associated with eigenvalues zero and infinity
#   vals <- v$values
#   
#   zerov <- min(abs(vals))
#   infv <- max(vals)
#   
#   vec <- v$vectors
#   
#   vec <- vec[,-c(which(vals==zerov), which(vals==(-zerov)), which(vals==infv))]  
#   # for i = 1:N
#   # vec(:,i) = vec(:,i) / norm(vec(:,i));
#   # end
#   
#   #normalize each vector - VERIFIED
#   for (i in 1:(N-2))
#     vec[,i] = vec[,i] / norm_vec(vec[,i])
#   
#   
#   #
#   
#   #   vec <- row.normalize(vec)
#   cost1 <- diag(t(vec) %*% L1 %*% vec);
#   cost2 <- diag(t(vec)%*%L2%*%vec)
#   #[Y,I] = sort(cost1,'ascend');
#   
#   I <- order(cost1)
#   
#   cost1 <- sort(cost1)
#   
#   #   cost2 = cost2(I);
#   cost2 <- cost2[I]
#   
#   plot(cost1, cost2, col="blue")
#   
#   #   plot <- ggplot(NULL, aes(cost1, cost2)) + t1
#   #   plot <-  plot + geom_point()
#   
#   # return(list(cost1=cost1, cost2 =cost2))
#   
#   #ex = []; # The cuts to be excluded
#   #U = []; # The output
#   ex = list() # The cuts to be excluded
#   U = list() # The output
#   iter = 0
#   pick = 0
#   #   return(U)
#   while (length(U) < 1){
#     
#     # Pick the smallest cut for Graph A, excluding those in ex
#     for (i in 1:(N-2))
#       #if (ismember(I[i],ex)==false) {
#       if (!(I[i] %in% unlist(ex))) {
#         U_ind = i
#         break
#       }
#     
#     print(U_ind)
#     
#     # Compute the Pareto frontier
#     cur_cost <- cost2[U_ind]
#     start <- U_ind
#     
#     #for (i in start:(N-1))
#     for (i in start:(N-3))
#       #if ((cost2[i+1] < cur_cost) && (ismember(I(i+1),ex)==false)) {
#       if ((cost2[i+1] < cur_cost) && (!(I[i+1] %in% unlist(ex)))) {
#         #U_ind = [U_ind, i+1]
#         U_ind = cbind(U_ind, i+1)
#         cur_cost <- cost2[i+1]
#       }
#     
#     #ex = [ex, I(U_ind)']; # Exclude chosen cuts
#     print(ex)
#     print(t(I[U_ind]))
#     
#     print(U_ind)
#     
#     ex[[length(ex)+1]] =t(I[U_ind])
#     U_ind <- as.matrix(U_ind)
#     if (iter > 0) { # Skip first pass
#       #         fprintf('iter:\t%d\n', iter);
#       #for i=1:size(U_ind,2)
#       for (i in 1:dim(U_ind)[2]){
#         #   if nnz(vec(:,I(U_ind(i)))>0)>0 && nnz(vec(:,I(U_ind(i)))<0)>0
#         #U = [U, vec(:,I(U_ind(i)))];
#         U[[length(U) + 1]] = vec[,I[U_ind[i]]]
#         pick <- pick + 1
#         #fprintf('%d\t%f\t%f\n', pick, cost1(U_ind(i)), cost2(U_ind(i)));
#         points(cost1[U_ind[i]],cost2[U_ind[i]], cex=1.5, col = "dark red");
#         #text(cost1(U_ind(i)),cost2(U_ind(i)),int2str(pick),'Color',[1 0 0]);
#         #   end
#       }       
#     }
#     
#     iter <- iter + 1;
#     
#   }
#   
#   print("The number of picks:")
#   print(pick)
#   
#   print("Exlcuded cuts:")
#   print(ex)
#   # stop("as")
#   
#   # return(U_ind)
#   U <- matrix(unlist(U), ncol = length(U), byrow = TRUE)
#   
#   print(dim(U))
#   
#   
#   
#   # prname <- "junit-4.12"
#   # ggsave(filename= paste("benchmark", prname ,"pareto_frontier.png", sep="/"), plot=plot, pointsize = 15, width = 10, height = 10)
#   
#   cost1 = diag(t(U)%*%L1%*%U);
#   cost2 = diag(t(U)%*%L2%*%U);
#   
# #   plot(cost1, cost2, col="yellow")
#   
#   return(U)
#   
#   # return(list(cost1=cost1, cost2=cost2))
#   
# }
