
build.partition <- function(d) {
  
  partition <- c()
  
  lista <- vector("list", length(d$lower))
  print(length(d$lower))
  
  for (i in 1:length(d$lower)) {
    lista[[i]] <- labels(d$lower[[i]])
    
  }
  
  for (i in 1:length(lista))
    for (j in 1:length(lista[[i]]))
     
    partition[lista[[i]][[j]]] <- i
  
    
return(partition)
  
}

#Input a matrix of groups to be merged
compute.join <- function(groups, threshold, prname) {
  
  print(dim(groups))
  
  agreement.matrix = matrix(0, nrow=dim(groups)[2], ncol=dim(groups)[2],
                      dimnames = list(colnames(groups), colnames(groups)))
  
   for (i in 1: dim(groups)[1]) {
    
     g <- normalizeVector(groups[i,])
     
     k <- max(g)
     
     for (j in 1:k) {
       
       subset <- g[which(g==j)]
       
       for (m in 1:length(subset))
         
         for (n in m:length(subset)) 
           if (n != m)
           {
           agreement.matrix[names(subset[m]), names(subset[n])] <-
                    agreement.matrix[names(subset[m]), names(subset[n])] +1
           agreement.matrix[names(subset[n]), names(subset[m])] <-
             agreement.matrix[names(subset[n]), names(subset[m])] +1
         }
         
     }
     
   }
 
  require(RColorBrewer)
  require(gplots)
  
  cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
  
  distCor <- function(x) as.dist(1-cor(t(x)))
  
  png(paste("benchmark", prname ,"dendrogram8.png", sep="/"), pointsize = 15, width = 3840, height = 3840)
  
  #h <- heatmap.2(mydata.Euc.dist.matrix, trace="none", scale="row", zlim=c(-3,3), reorder=TRUE,
  #          distfun=distCor, hclustfun=hclustComplete, col=rev(cols), symbreak=FALSE)
  
  #h <-heatmap.2(agreement.matrix, trace="none", scale="row", zlim=c(-3,3), reorder=TRUE,
   #             distfun=distCor, hclustfun=hclustComplete, col=rev(cols), symbreak=FALSE)
  hclustfunc <- function(x) hclust(x, method="average")
  distfunc <- function(x) as.dist((1-cor(t(x)))/2)
  d <- distfunc(agreement.matrix)
  dend <- hclustfunc(d)
  
  print(threshold)
  
  noelems = dim(agreement.matrix)[1]
  
  t <- threshold
  
  height = t / (noelems*15)
  
  cut.dend <- cut(as.dendrogram(dend), h = height)
  
  while(length(cut.dend$lower) > threshold){
    
    height = t / (noelems*15)
    
    cut.dend <- cut(as.dendrogram(dend), h = height)
    #cut.dend <-cut.dendrogram(as.dendrogram(dend), h = threshold)
    t <- t + 1
    
    print(paste("printing no of lower branches:", length(cut.dend$lower)))
    
    
  }
  
  #stop("sjh")
  
  print(paste("The dendrogram was cut at height:", height  , sep=" "))
  merged.dend = data.frame()
  #merged.dend <- merge(cut.dend$lower)
  merged.dend <- do.call(merge, cut.dend$lower)

  p <- plot(merged.dend, nodePar = list(pch = c(1,7), col = 2:1))
  
  dev.off()
  

  joined.partition <- build.partition(cut.dend)
  
  return(joined.partition)
  
}


compute.multiobj.join <- function(prname, median) {
  
  
  setwd("~/workspace")
  
  groups <- read.csv(paste("benchmark", prname ,"multiobj-groups.csv", sep="/")
                     ,check.names=FALSE, sep=",",  header = TRUE)
  
  rownames(groups) <- groups$X
  groups <- groups[,-1]
  
  #Load the results
  results <- read.csv(paste("benchmark", prname ,"multiobj-results.csv", sep="/"), sep=",",  header = TRUE)
  
  results <- results[,-1]
  
  results <- cbind(results, groups)
  
  results <- as.matrix(results)
  
  
  #results <- results[!duplicated(results[,1:3]),]
  
  results <- unique( results[ , 1:dim(results)[2]] )

  groups <- results[ , 4:dim(results)[2]]
  
  p <- compute.join(groups, median, prname)
  
  #Load the priori decomposition
  decomposition <- read.csv(paste("benchmark", prname ,"decomposition.csv", sep="/"), sep=",",  header = TRUE)
  priori.decomp <- decomposition$x
  names(priori.decomp) <- decomposition$X
  priori.decomp <- normalizeVector(priori.decomp)
  
  
  print(names(p))
  
  print(names(priori.decomp))
  
  priori.decomp <- find.intersection(priori.decomp, p)

  
  priori.decomp <- normalizeVector(priori.decomp)
  
  mojo <- compute.MoJo(p, priori.decomp)
  
  mojosim <- sapply(mojo, function(m) 1 - (m/length(p)))
  
  print(mojosim)
  
  return(p)
}