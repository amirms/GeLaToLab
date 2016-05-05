#computes the Pareto front for multi view learning
# Input: a list of graph laplacians for each view

#compute the normalized laplacians on the kernels

generalized_pareto_multiview <- function(L1, L2, L3) {
  
  require(geigen)
  require(ggplot2)
  
  avg_lap <- function(Laps) {
    apply(simplify2array(Laps), 1:2, mean)
  }
  
  #TODO check the dimensions match

  N = dim(L1)[1]
  
  norm_vec <- function(x) sqrt(sum(x^2))
  
  Laps <- list(L1, L2, L3)
  
  compute_all_cuts <- function(L1, L2){
    v = geigen(L1, L2) #VERIFIED
    #   v = eigen(solve(L2)%*%L1) 
    #Remove the eigenvectors associated with eigenvalues zero and infinity
    vals <- v$values
# #     
    zerov <- min(abs(vals))
    infv <- max(vals)
    
    vec <- v$vectors
    
    vec <- vec[,-c(which(vals==zerov), which(vals==(-zerov)), which(vals==infv))]  
#   vec <- vec[,-c(which(vals==zerov), which(vals==(-zerov)))]  
    
    #normalize each vector - VERIFIED
    for (i in 1:dim(vec)[2])
      vec[,i] = vec[,i] / norm_vec(vec[,i])
  
    vec
  }

  is_dominated <- function(cut, included_cuts) {
    
    cost1 <- cut$cost1
    cost2 <- cut$cost2
    cost3 <- cut$cost3  
    print(cost1)
    print(cost2)
    print(cost3)
    
    for(i in 1:length(included_cuts)){
      
      cur_cost1 <- included_cuts[[i]]$cost1
      cur_cost2 <- included_cuts[[i]]$cost2
      cur_cost3 <- included_cuts[[i]]$cost3
      
      print(cur_cost1)
      print(cur_cost2)
      print(cur_cost3)
      
      print((cur_cost1 < cost1) & (cur_cost2 < cost2) & (cur_cost3 < cost3))
      
      if ((cur_cost1 < cost1) & (cur_cost2 < cost2) & (cur_cost3 < cost3))
        return(TRUE)
 
    }
    
    return(FALSE)
  }

#return Indices
dominates <- function(cut, included_cuts) {
  
  cost1 <- cut$cost1
  cost2 <- cut$cost2
  cost3 <- cut$cost3   
  
  indices <- c()
  
  for(i in 1:length(included_cuts)){
    
    cur_cost1 <- included_cuts[[i]]$cost1
    cur_cost2 <- included_cuts[[i]]$cost2
    cur_cost3 <- included_cuts[[i]]$cost3
    
    if ((cur_cost1 >= cost1) & (cur_cost2 >= cost2) & (cur_cost3 >= cost3))
      indices <- c(indices, i)
    
  }
  
  indices    
}


  compute_pareto_front <- function(cost1, cost2, cost3, I, vec) {
    #TODO cost1 and cost2 should have the same length
    N <- length(cost1)
    
    ex = list() # The cuts to be excluded
    U = list() # The output
    iter = 0
    pick = 0
    
    included_cuts <- list()
    
    included_cuts[[length(included_cuts) + 1]] <- list(cost1 = cost1[1], cost2 = cost2[1], cost3=cost3[1], U=vec[,1])
    
    for (i in 2 : N){
      
      print(paste(cost1[i], cost2[i], cost3[i]))
      print(paste(included_cuts[[1]]$cost1, included_cuts[[1]]$cost2, included_cuts[[1]]$cost3))
      
      if(!is_dominated(list(cost1=cost1[i], cost2=cost2[i], cost3=cost3[i]), included_cuts)) {
        
        indices <- dominates(list(cost1=cost1[i], cost2=cost2[i], cost3=cost3[i]), included_cuts)
        
      if (!is.null(indices))
        included_cuts <- included_cuts[-indices]
      
        included_cuts[[length(included_cuts) + 1]] <- list(cost1 = cost1[i], cost2 = cost2[i], cost3=cost3[i], U=vec[,i])
      
      }
      
    }
    
    
    included_cuts
    
    
  }

  vec1 <- compute_all_cuts(L1, avg_lap(list(L2, L3)))
  vec2 <- compute_all_cuts(L2, avg_lap(list(L1, L3)))
  vec3 <- compute_all_cuts(L3, avg_lap(list(L1, L2)))
  
  vec <- cbind(vec1, vec2, vec3)
  
  #   print(dim(L1))
  #   print(dim(vec))
  
  #   stop("as")
  
  #   vec <- row.normalize(vec)
  cost1 <- diag(t(vec) %*% L1 %*% vec);
  cost2 <- diag(t(vec)%*%L2%*%vec)
  cost3 <- diag(t(vec)%*%L3%*%vec)
  #[Y,I] = sort(cost1,'ascend');
  
  require(scatterplot3d)
  s3d <- scatterplot3d(cost1,cost2,cost3, main="Pareto Frontier", color="blue",  pch=3)
  
#   require(rgl)
#   is3d <- plot3d(cost1,cost2,cost3,  col="red", size=10)
  
  
  I <- order(cost1)
  cost1 <- sort(cost1)
  #   cost2 = cost2(I);
  cost2 <- cost2[I]
  cost3 <- cost3[I]
  
  
  pareto_front <- compute_pareto_front(cost1, cost2, cost3, I, vec)

  print("number of cuts")
  print(length(pareto_front))

  pareto_cost1 <- unlist(lapply(pareto_front, function(p) p$cost1  ))
  pareto_cost2 <- unlist(lapply(pareto_front, function(p) p$cost2  ))
  pareto_cost3 <- unlist(lapply(pareto_front, function(p) p$cost3  ))
  
  s3d$points3d(pareto_cost1, pareto_cost2, pareto_cost3, cex=1.5, pch=16, col = "dark red")

  pareto_U <- lapply(pareto_front, function(p) p$U)


  pareto_U <- matrix(unlist(pareto_U), ncol = length(pareto_U), byrow = TRUE)

rownames(pareto_U) <- rownames(L1)
  return(pareto_U)

#   
#   compute_pareto_front2 <- function(cost1, cost2, I, L1, L2, vec){
#     
#     #TODO cost1 and cost2 should have the same length
#     N <- length(cost1)
#     
#     ex = list() # The cuts to be excluded
#     U = list() # The output
#     iter = 0
#     pick = 0
#     
#     #   return(U)
#     while (length(U) < 1){
#       
#       # Pick the smallest cut for Graph A, excluding those in ex
#       for (i in 1:N)
#         #if (ismember(I[i],ex)==false) {
#         if (!(I[i] %in% unlist(ex))) {
#           U_ind = i
#           break
#         }
#       
#       print(U_ind)
#       
#       # Compute the Pareto frontier
#       cur_cost <- cost2[U_ind]
#       #     cur_cost3 <- cost3[U_ind]
#       
#       start <- U_ind
#       
#       for (i in start:(N-1))
#         #if ((cost2[i+1] < cur_cost) && (ismember(I(i+1),ex)==false)) {
#         if ((cost2[i+1] < cur_cost) && (!(I[i+1] %in% unlist(ex)))) {
#           #U_ind = [U_ind, i+1]
#           U_ind = c(U_ind, i+1)
#           cur_cost <- cost2[i+1]
#         }
# 
#       ex[[length(ex)+1]] =t(I[U_ind])
# #       U_ind <- as.matrix(U_ind)
#       if (iter > 0) { # Skip first pass
#         #         fprintf('iter:\t%d\n', iter);
#         #for i=1:size(U_ind,2)
#         for (i in 1:length(U_ind)){
#           #   if nnz(vec(:,I(U_ind(i)))>0)>0 && nnz(vec(:,I(U_ind(i)))<0)>0
#           #U = [U, vec(:,I(U_ind(i)))];
#           U[[length(U) + 1]] = vec[,I[U_ind[i]]]
#           pick <- pick + 1
#           #fprintf('%d\t%f\t%f\n', pick, cost1(U_ind(i)), cost2(U_ind(i)));
#           #           s3d$points3d(cost1[U_ind[i]],cost2[U_ind[i]], cost3[U_ind[i]], cex=1.5, col = "dark red");
#           print("printing the costs")
#           print(c(cost1[U_ind[i]],cost2[U_ind[i]]))
#           #text(cost1(U_ind(i)),cost2(U_ind(i)),int2str(pick),'Color',[1 0 0]);
#           #   end
#         }       
#       }
#       
#       iter <- iter + 1;
#       
#     }
#     
# #     return(ex)
# #     
# #     print("printing all the indices")
# #     print(U_ind)
# #     
# # #     return(U_ind)
# #     
#     return(list(U=U, U_ind=U_ind))
#     
# #     print("The number of picks:")
# #     print(pick)
#     
# #     print("Exlcuded cuts:")
# #     print(ex)
#     # stop("as")
#     
#     # return(U_ind)
#     U <- matrix(unlist(U), ncol = length(U), byrow = TRUE)
#     
#     cost1 = diag(t(U)%*%L1%*%U);
#     cost2 = diag(t(U)%*%L2%*%U);
# #     cost3 = diag(t(U)%*%L3%*%U);
#     
#     
#     return(list(U = U, cost1 = cost1, cost2=cost2, cost3=cost3))
#   }
# 
#   vec1 <- compute_all_cuts(L1, avg_lap(list(L2, L3)))
#   vec2 <- compute_all_cuts(L2, avg_lap(list(L1, L3)))
#   vec3 <- compute_all_cuts(L3, avg_lap(list(L1, L2)))
#   
#   vec <- cbind(vec1, vec2, vec3)
#   
# #   print(dim(L1))
# #   print(dim(vec))
#   
# #   stop("as")
#   
# #   vec <- row.normalize(vec)
#   cost1 <- diag(t(vec) %*% L1 %*% vec);
#   cost2 <- diag(t(vec)%*%L2%*%vec)
#   cost3 <- diag(t(vec)%*%L3%*%vec)
#   #[Y,I] = sort(cost1,'ascend');
# 
# require(scatterplot3d)
# s3d <- scatterplot3d(cost1,cost2,cost3, main="Pareto Frontier")
# 
# require(rgl)
# is3d <- plot3d(cost1,cost2,cost3,  col="red", size=10)
# 
# 
#   I <- order(cost1)
#   cost1 <- sort(cost1)
# #   cost2 = cost2(I);
#   cost2 <- cost2[I]
#   cost3 <- cost3[I]
# 
# 
# U_1<- compute_pareto_front(cost1, cost2, I, L1, L2,  vec)
# U_2 <- compute_pareto_front(cost1, cost3, I, L1, L3, vec)
# 
# 
# U_indices <- union(U_1$U_ind, U_2$U_ind)
# 
# 
# 
# s3d$points3d(cost1[U_indices],cost2[U_indices], cost3[U_indices], cex=1.5, col = "dark red")
# 
# #compute the third duo-dimension
# I <- order(cost2)
# cur_cost2 <- sort(cost2)
# 
# #   cost2 = cost2(I);
# cur_cost3 <- cost3[I]
# cur_cost1 <- cost1[I]
# 
# U_3 <- compute_pareto_front(cur_cost2, cur_cost3, I, L2, L3, vec)
# 
# t <- 0
# for(i in 1:length(U_3$U_ind)){
#   
#   index <- U_3$U_ind[i]
#   
#   if (!((cur_cost2[index] %in% cost2[U_indices]) && (cur_cost3[index] %in% cost3[U_indices])))
#       t <- t + 1
#   
# }
# 
# s3d$points3d(cur_cost1[U_3$U_ind],cur_cost2[U_3$U_ind], cur_cost3[U_3$U_ind], cex=1.5, size=10, col = "dark red")
# print("printing the number of cuts")
# 
# print(length(U_indices) + t)
# 
# # Us <- union(U_1$U, U_2$U, U_3$U)
# 
# Us <- Reduce(union,  list(U_1$U, U_2$U, U_3$U))
# 
# U <- matrix(unlist(Us), ncol = length(Us), byrow = TRUE)
# 
# rownames(U) <- rownames(L1)
# 
# 
# return(U)

# print(list(cost1[U_1],cost2[U_1], cost3[U_1]))
# 
# print(list(cost1[U_2],cost2[U_2], cost3[U_2]))


# 
# print(list(cost1[U_3],cost2[U_3], cost3[U_3]))
# 
# 
# s3d$points3d(cost1[U_3],cost2[U_3], cost3[U_3], cex=1.5, col = "dark red")
# 
# U_indices <- Reduce(union,  list(U_1, U_2, U_3))

#   plot(cost1, cost2, col="blue")







  





# prname <- "junit-4.12"
# ggsave(filename= paste("benchmark", prname ,"pareto_frontier.png", sep="/"), plot=plot, pointsize = 15, width = 10, height = 10)





return(U)

# return(list(cost1=cost1, cost2=cost2))

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
