#Input: A reference named vector ref
# A named vector vec for which MoJo is to be compared relative to ref
#precondition: Both vectors are normalized and of equal length
compute.MoJo <-function(grp, ref){
  
  require(igraph)
  
  if (length(grp) != length(ref))
    stop("Incompatible Length of Vectors")
  
  k.grp <- max(grp)
  k.ref <- max(ref)
  
  v = matrix(0, nrow = k.grp, ncol=k.ref)
  for (i in 1:k.grp)
    for (j in 1:k.ref)
      v[i,j] <- length(intersect(names(grp[which(grp==i)]), names(ref[which(ref==j)])))

  
  max.v <- apply(v, 1, function(xs){max(xs)})
  
  nmoves <- length(grp) - sum(max.v)
  
  max.k <- apply(v, 1, function(xs){which(xs==max(xs))})
  
  edges <- c()
  
  for (i in 1:length(max.k)) 
    
    # i is the src.element
    for (j in 1:length(max.k[[i]]))
      
      edges <- rbind(edges, c(i, max.k[[i]][j] + k.grp))
  
  gr <- graph.data.frame(edges, directed=TRUE, vertices=NULL)
  
  V(gr)$type = c(rep(FALSE, length(max.k)), rep(TRUE, vcount(gr) - length(max.k)))
  
  mmg <- maximum.bipartite.matching(gr)$matching
  
  l.join = 0
  
  for (i in 1:k.grp)
    if (!is.na(mmg[i]))
        l.join <- l.join+1
  
  njoins = k.grp - l.join

  return(nmoves + njoins)
  
}


build.decomposition <- function(dir, pattern) {
  decomposition <- decompose.folder(dir, pattern, 1)
  
  return(decomposition$group)
  
  
}

decompose.folder <- function(dir, pattern, index) {
  result = list(
    group = NULL,
    idx = NULL
  )
  
  all <- list.files(dir, full.names = TRUE, recursive = FALSE)
  
  files <- grep(pattern, all,value = TRUE)
  
  if (length(files) > 0)
      result$group <- group.files(grep(pattern, all, value = TRUE), index)
      
  dirs <- all[file.info(all)$isdir]
  
  result$idx = index + 1
  
  i = 1
  
  print(length(dirs))
  
  while(i <= length(dirs)) {
    decomposition <- decompose.folder(dirs[i], pattern, result$idx) 
    
    result$idx <- decomposition$idx
    
    result$group <- c(result$group, decomposition$group)
    
    i <- i+1
  }
  
    return(result)
  
}

group.files <- function(files, index) {
  
  group <- rep(index, length(files))
  
  names(group) <- files
  
  print(group)
  
  return(group)
  
}

find.intersection <- function(v, g) {
  
  return(v[intersect(names(v),names(g))])
  
}

extract.authoritative.partition <- function(prname) {

  setwd("~/workspace")
    
  require(igraph)
  require(Rcpp)
    
  #Load the priori decomposition
  decomposition <- read.csv(paste("benchmark", prname ,"decomposition.csv", sep="/"), sep=",",  header = TRUE)
  priori.decomp <- decomposition$x
  names(priori.decomp) <- decomposition$X
  priori.decomp <- normalizeVector(priori.decomp)
  
  print(max(priori.decomp))
  
  
}

compute.MoJoSim <- function(clusters, classes) {
  N <- length(clusters)
  
  1 - (compute.MoJo(clusters, classes)/N)
  
}

compute.cluster.purity <- function(clusters, classes) {
  sum(apply(table(classes, clusters), 2, max)) / length(clusters)
}

compute.NMI <- function(clusters, classes) {
  require(clue)
  
  cl_agreement(as.cl_partition(clusters), as.cl_partition(classes), method="NMI")
  
}

# Rand Index
compute.RI <- function(clusters, classes) {
  require(clue)

  cl_agreement(as.cl_partition(clusters), as.cl_partition(classes), method="Rand")
  
}

# # Rand Index
compute.AdjRI <- function(clusters, classes) {
  require(clue)
  
  cl_agreement(as.cl_partition(clusters), as.cl_partition(classes), method="cRand")
  
}

compute.f1 <- function(actual, predict) {
  precision <- compute.precision(actual, predict)
  recall <- compute.recall(actual, predict)
  f1 <- (2*precision*recall)/(precision+recall)
  mean(f1)
}

compute.recall <- function(actual, predict) {
  confusion <- table(actual, predict)
  confusion <- cbind(confusion,0)
  relativeDocs <- apply(confusion[,],1, FUN = sum)
  recall <- apply(confusion[,],2,function(x) if(relativeDocs[which.max(x)]) max(x)/relativeDocs[which.max(x)] else 0)
  mean(recall[1:length(recall)-1])
}

compute.precision <- function(actual, predict) {
  confusion <- table(actual, predict)
  confusion <- cbind(confusion,0)
  precision <- apply(confusion[,],2,function(x) if(sum(x)>0)max(x)/sum(x) else 0)
  mean(precision[1:length(precision)-1])
}
