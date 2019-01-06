#Input: a java project comprising of BoW, transaction histoy, and dependency graph
#Output: a list of classnames to be used as the intersection for all the experiments written to a file
compute_list_of_classnames <- function(prname){
#   require(proxy)
  
  setwd("~/workspace")
  
  #Load the priori decomposition
  decomposition <- read.csv(paste("benchmark", prname ,"decomposition.csv", sep="/"), sep=",",  header = TRUE)
  priori.decomp <- decomposition$x
  names(priori.decomp) <- decomposition$X
  priori.decomp <- normalizeVector(priori.decomp)
  
  #ELIMINATE packages with smaller than 4 classes
  noc <- max(priori.decomp)
  
  priori.decomp.names <- unlist(lapply(1:noc, function(x) {
    cls = priori.decomp[priori.decomp==x]
    if (length(cls) > 3)
      names(cls)
    else
      c()
  }))
  
  
  #Load the adjacency matrix
  extensions= c("java/", "org/xml/", "javax/")
  cfg <- import.bunch.matrix(paste("benchmark", prname ,"dep_graph.txt", sep="/"), exclude.ext=extensions)
  #cfg <- read.table("benchmark/jedit-5.1.0/cfg.csv", sep=",", row.names = 1, header = TRUE, check.names = FALSE)
  #   cfg <- unweight.adjacency(cfg)
  # no self-references
  diag(cfg) <- 0 
  cfg <- make.symmetric(cfg)
  d <- apply(abs(cfg),1,sum)
  
  indices <- which(d==0)
  
  while(length(indices) > 0) {
    cfg <- cfg[-indices, -indices]
    d <- apply(abs(cfg),1,sum)
    indices <- which(d==0)
  }
 
  #Load the transaction frequency
  freq <- read.table(paste("benchmark", prname , "mydata-change-freq-matrix.csv", sep="/"), sep=",", row.names = 1, header = TRUE, check.names = FALSE)  
  freq <- as.matrix(freq)
  
  #Load the bag of words
  BoW <- read.table(paste("benchmark", prname , paste("BoW", paste(prname, "BoW.csv", sep="-"), sep="/"), sep="/"), sep=",", row.names = 1, header = TRUE, check.names = FALSE)  
  
    
  classnames <- compute_intersection_names(list(rownames(cfg), rownames(BoW), rownames(freq), priori.decomp.names ))
  
  
  write.table(classnames, paste("benchmark", prname , paste("MULTIVIEW", "classnames.txt" ,sep="/") , sep="/"))
  
  
  
}

#Input: a list of string vectors
#Output: intersect of all names
compute_intersection_names <- function(names){
  intersect_names <- names[[1]]
  no_names <- length(names)

  if (no_names < 2)
    return(intersect_names)
  
  for (i in 2:no_names){
      
    intersect_names <- intersect(intersect_names, names[[i]])    
      
  }
  
  return(intersect_names)
  
}



cotraining <- function(Ks, iters, storeEachIteration=FALSE, k=10) {
    Ls <- lapply(Ks, function(K) laplacian(K, TRUE))
    
    es <- lapply(Ls, function(L) eigen(L))
    
    os <- lapply(es, function(e) order(e$values, decreasing=FALSE)[1:k] )
    
    Us <- list()
    
    for (i in 1:length(es))
      Us[[i]] <- es[[i]]$vectors[,os[[i]]]
    #   Us <- mapply(function(e,o) e$vectors[,o], es, os)
    # o1 <- 
    #   U1 <- e1$vectors[,o1]
    
    #   e2 <- eigen(L2)
    #   o2 <- order(e2$values, decreasing=FALSE)[1:k]
    #   U2 <- e2$vectors[,o2]
    
    #   result2 <- c(result2, cluster.mojosim(U2, k, priori.decomp))
    
    #   print(dim(Us))
    #   return(list(U1=U1, U2=U2))
    
    project_sum <- function(Us, K) {
      sum_UUT <- Us[[1]] %*% t(Us[[1]])
      if (length(Us) > 1)
        for (i in 2:length(Us))
          sum_UUT <- sum_UUT + (Us[[i]] %*% t(Us[[i]]))
        
        sum_UUT %*% K
    }
    
    Ss = list()
    
    make_symmetric <- function(U) {
      0.5 * (U + t(U))
    }
    
    SS_per_iteration=list()
    
    SS_per_iteration[[1]] = Ks
    
    for (iter in 1:iters) {
      print(paste("iteration", iter))
      
      for (i in 1:length(Us))
        Ss[[i]] <- make_symmetric(project_sum(Us[-i], Ks[[i]]))
      
      SS_per_iteration[[iter+1]] = Ss
      
      #     Ks <- lapply(make.symmetric(U2 %*% t(U2) %*% K1)
      
      Ls <- lapply(Ss, function(S) laplacian(S, TRUE))
      
      es <- lapply(Ls, function(L) eigen(L))
      os <- lapply(es, function(e) order(e$values, decreasing=FALSE)[1:k] )    
      for (i in 1:length(es))
        Us[[i]] <- es[[i]]$vectors[,os[[i]]]
      
    }

    if (!storeEachIteration){
    return(Ss)
    }
    
    return(SS_per_iteration)
}

plot.cotraining.per.iteration = function(Ks, iters) {
  priori.decomp <- build.dendrogam(rownames(Ks[[1]]))
  
  Rs <- cotraining(Ks, iters, TRUE)
  
  PDs <- list()
  for (i in 1:length(Rs)){
    Rss <- Rs[[i]]
    PDs[[i]]<- rep(0, length(Rss))
    for (j in 1:length(Rss)){
      kernel = Rss[[j]]
      #compute distance from kernel
      myDist <- squared.euclidean.distance.of.kernel.matrix(kernel)
      myDist <- as.dist(myDist)
      
      # pinned it to complete linkage
      clusters <- hclust(myDist, method = 'complete')
      
      clusters.tree <- ape::as.phylo(clusters)
      priori.tree <- priori.decomp$tree
      path.difference <- phangorn::path.dist(clusters.tree, priori.tree, check.labels = T)
      PDs[[i]][j] <- path.difference
    }
  }
  
  
  
  Iterations = 1:length(PDs)
  View <- unlist(lapply(seq(length(PDs[[1]])), function(index) paste("View", index, sep="")))
  #   View = c("View1","View2", "View3")
  xxxd = expand.grid(Iterations=Iterations, View=View)
  
  view1_results <- unlist(lapply(PDs, function(PD_iter) PD_iter[1]))
  view2_results <- unlist(lapply(PDs, function(PD_iter) PD_iter[2]))
  view3_results <- unlist(lapply(PDs, function(PD_iter) PD_iter[3]))
  
  xxxd$PD = c(view1_results, view2_results, view3_results)

  results <- c(view1_results, view2_results, view3_results)
  min_y_lim <- min(unlist(results)) - 0.025
  max_y_lim <- max(unlist(results)) + 0.025
  
  
  #   plot(x=1:length(result1), y=result1)
  #   plot(x=1:length(result2), y=result2)
  
  require(ggplot2)
  t1<-theme(                              
    plot.background = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    panel.background = element_blank(),
    axis.line = element_line(size=.5),
    axis.text=element_text(size=14),
    axis.title=element_text(size=14,face="bold"),
    legend.title=element_text(size=12, vjust=-12)
    #     guide_colourbar.title = element_text(draw.ulim = FALSE, draw.llim = FALSE)
    #    legend.key.height=unit(3,"line"),
    #    legend.key.width=unit(3,"line")
  )
  
  plot <- ggplot(xxxd, aes(Iterations, PD, colour = View))+ xlim(1,length(view1_results)) +t1 + ylim(min_y_lim, max_y_lim)  + geom_line(size=1.5) 
}


#Input: Similarity matrix for both views: K1, K2 (both K's are symmetric)
#       k number of clusters
#       iters: number of iterations
#Output: Assignments to k clusters

compute.cotraining <- function(Ks, k, iters, priori.decomp, prname) {
  #Compute normalized lap
#   normalized_laplacian <- function(A) {
#     D = diag(rowSums(A))
#     
#     sqrt.inv.D <- solve(sqrt(D))
#     return(sqrt.inv.D%*%A%*%sqrt.inv.D)
#   }
  
#   result1 = c(0)
#   result2 = c(0)
  
#   L1 <- compute.igraph.laplacian(K1)
#   L2 <- compute.igraph.laplacian(K2)
  
  Ls <- lapply(Ks, function(K) laplacian(K, TRUE))
  
  es <- lapply(Ls, function(L) eigen(L))
  os <- lapply(es, function(e) order(e$values, decreasing=FALSE)[1:k] )
  
  Us <- list()
  
  for (i in 1:length(es))
    Us[[i]] <- es[[i]]$vectors[,os[[i]]]
#   Us <- mapply(function(e,o) e$vectors[,o], es, os)
# o1 <- 
#   U1 <- e1$vectors[,o1]
  
#   e2 <- eigen(L2)
#   o2 <- order(e2$values, decreasing=FALSE)[1:k]
#   U2 <- e2$vectors[,o2]
  
  results <- lapply(Us, function(U) cluster.mojosim(U, k, priori.decomp))

#   result2 <- c(result2, cluster.mojosim(U2, k, priori.decomp))
  
#   print(dim(Us))
#   return(list(U1=U1, U2=U2))

  project_sum <- function(Us, K) {
    sum_UUT <- Us[[1]] %*% t(Us[[1]])
    if (length(Us) > 1)
      for (i in 2:length(Us))
        sum_UUT <- sum_UUT + (Us[[i]] %*% t(Us[[i]]))
    
    sum_UUT %*% K
  }

  Ss = list()

  make_symmetric <- function(U) {
    0.5 * (U + t(U))
  }

  for (i in 1:iters) {
    print(paste("iteration", i))
    
    for (i in 1:length(Us))
      Ss[[i]] <- make_symmetric(project_sum(Us[-i], Ks[[i]]))
    
#     Ks <- lapply(make.symmetric(U2 %*% t(U2) %*% K1)
    
      Ls <- lapply(Ss, function(S) laplacian(S, TRUE))
      
      es <- lapply(Ls, function(L) eigen(L))
      os <- lapply(es, function(e) order(e$values, decreasing=FALSE)[1:k] )    
      for (i in 1:length(es))
        Us[[i]] <- es[[i]]$vectors[,os[[i]]]
 
    for (i in 1:length(results))
      results[[i]] <- c(results[[i]], cluster.mojosim(Us[[i]], k, priori.decomp))
    
#     result1 <- c(result1, cluster.mojosim(U1, k, priori.decomp))
#     result2 <- c(result2, cluster.mojosim(U2, k, priori.decomp))
  }

#   r1 <- list(Iteration = 1:length(result1), View=rep(1, length(result1)), MoJoSim = result1)
#   r2 <- list(Iteration = 1:length(result2), View=rep(2, length(result1)), MoJoSim = result2)

  
  Iterations = 1:length(results[[1]])
  View <- unlist(lapply(seq(length(results)), function(index) paste("View", index, sep="")))
#   View = c("View1","View2", "View3")
  xxxd = expand.grid(Iterations=Iterations, View=View)
   xxxd$MoJoSim = c(results[[1]], results[[2]], results[[3]])
# xxxd$MoJoSim <- results


min_y_lim <- min(unlist(results)) - 0.025
max_y_lim <- max(unlist(results)) + 0.025


#   plot(x=1:length(result1), y=result1)
#   plot(x=1:length(result2), y=result2)
  
require(ggplot2)
t1<-theme(                              
  plot.background = element_blank(), 
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(), 
  panel.border = element_blank(), 
  panel.background = element_blank(),
  axis.line = element_line(size=.4),
  axis.text=element_text(size=12),
  axis.title=element_text(size=14,face="bold"),
  legend.title=element_text(size=10, vjust=-12)
  #     guide_colourbar.title = element_text(draw.ulim = FALSE, draw.llim = FALSE)
  #    legend.key.height=unit(3,"line"),
  #    legend.key.width=unit(3,"line")
)

plot <- ggplot(xxxd, aes(Iterations, MoJoSim, colour = View))+ xlim(1,length(results[[1]])) +t1 + ylim(min_y_lim, max_y_lim)  + geom_line() 

# ggsave(filename= "iterations.png", plot=plot, pointsize = 15, width = 10, height = 10)
print(getwd())
# ggsave(filename= paste("..", "benchmark", prname ,"iterations.png", sep="/"), plot=plot, pointsize = 15, width = 10, height = 10)

ggsave(filename= "iterations.png", plot=plot, pointsize = 15, width = 10, height = 10)


# dev.off()

  Us <- lapply(Us, function(U) row.normalize(U[, 1:k]))

#   Us[[1]] <- row.normalize(U1[, 1:k])
#   U2 <- row.normalize(U2[, 1:k])
  
  U <- cbind(Us[[1]], Us[[2]], Us[[3]])
  
  print(results[[1]])
print(results[[2]])
print(results[[3]])

  return(U)
  
}


cluster.mojosim <- function(U, k, priori.decomp) {
  U <- U[,1:k]
  U <- row.normalize(U)
  
  result <- kmeans(U, centers = k, iter.max = 300, nstart = 1200)$cluster
  
  result <- normalizeVector(result)
  
#   print(rownames(U))
#   print(names(priori.decomp))
  
#   stop("as")
  
  names(result) = names(priori.decomp)
  
#   return(compute.NMI(result, priori.decomp))
  
  mojo = compute.MoJo(result, priori.decomp)

N <- length(priori.decomp)

  return(1 - (mojo/N))
  
}