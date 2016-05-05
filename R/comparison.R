compare.adjacency.matrices <- function(prname, neighbors, size , ncpus) {
  
  setwd("~/workspace")
  
  require(igraph)
  require(Rcpp)
  require(gelato)
  
  #Load the priori decomposition
  decomposition <- read.csv(paste("benchmark", prname ,"decomposition.csv", sep="/"), sep=",",  header = TRUE)
  priori.decomp <- decomposition$x
  names(priori.decomp) <- decomposition$X
    
  #Load the adjacency matrix
  extensions= c("java/", "org/xml/", "javax/")
  cfg <- import.bunch.matrix(paste("benchmark", prname ,"dep_graph.txt", sep="/"), exclude.ext=extensions)
  
  #Build mydata
  mydata = list(cfg=as.matrix(cfg))
  mydata$cfg <- remove.diagonal(mydata$cfg)

  #Fix the priori decomposition
  priori.decomp <- priori.decomp[intersect(names(priori.decomp),rownames(mydata$cfg))]
  priori.decomp <- normalizeVector(priori.decomp)
  
  mydata$cfg <- mydata$cfg[intersect(names(priori.decomp),rownames(mydata$cfg)), intersect(names(priori.decomp),colnames(mydata$cfg))]
      
  cfg = mydata$cfg
  
  apply.hill.climb(mydata, eval.func = Cpp.MQ.evaluator, neighbors=neighbors, 
                      size = size, ncpus=ncpus, priori.decomp)
  
  
  mydata$cfg <- normalize.outdegree(cfg)
  apply.hill.climb(mydata, eval.func = Cpp.MQ.evaluator, neighbors=neighbors, 
                      size = size, ncpus=ncpus, priori.decomp)
  
  
  mydata$cfg <- unweight(cfg)
  apply.hill.climb(mydata, eval.func = Cpp.MQ.evaluator, neighbors=neighbors, 
                      size = size, ncpus=ncpus, priori.decomp)
  

}


apply.hill.climb <- function(mydata, eval.func, neighbors, size, ncpus, decomp) {
  
  
  clusters <- hill.climbing.search(mydata, rownames(mydata$cfg), eval.func, ncpus, size, neighbors)
  
  N = length(decomp)
  
  mojo <- compute.MoJo(clusters$groups, decomp)
  
  mojosim <- sapply(mojo, function(m) 1 - (m/N))
  
  print("printing mojosim")
  print(mojosim)
  
}


normalize.outdegree <- function(A) {
  D = rowSums(A)
  for (i in 1:dim(A)[1])
    A[i,] <- A[i,]/D[i]
  
  return(A)
  
}

unweight.adjacency <- function(A) {
  
  for (i in 1:dim(A)[1])
    for (j in 1:dim(A)[2]) {
      if (A[i,j] > 0)
        A[i,j] = 1
    }
  
  return(A)
  
}


compare.lexsim.matrices <- function(prname, dims, neighbors, size , ncpus) {
  
  setwd("~/workspace")
  
  require(igraph)
  require(Rcpp)
  require(gelato)
  require(lsa)
  
  #Load the priori decomposition
  decomposition <- read.csv(paste("benchmark", prname ,"decomposition.csv", sep="/"), sep=",",  header = TRUE)
  priori.decomp <- decomposition$x
  names(priori.decomp) <- decomposition$X
  
 
  lexsim <- read.table(paste("benchmark", prname , "mydata-Cos-similarity-matrix.csv", sep="/"), sep=",", row.names = 1, header = TRUE, check.names = FALSE)
    
  #Build mydata
  mydata = list(lexsim=as.matrix(lexsim))
  mydata$lexsim <- remove.diagonal(mydata$lexsim)
  
  #Fix the priori decomposition
  priori.decomp <- priori.decomp[intersect(names(priori.decomp),rownames(mydata$lexsim))]
  priori.decomp <- normalizeVector(priori.decomp)
  
  mydata$lexsim <- mydata$lexsim[intersect(names(priori.decomp),rownames(mydata$lexsim)), intersect(names(priori.decomp),colnames(mydata$lexsim))]
  
  
  print(dim(mydata$lexsim))
  
  apply.hill.climb(mydata, rownames(mydata$lexsim), eval.func = Cpp.LQ.evaluator, neighbors=neighbors, 
                   size = size, ncpus=ncpus, priori.decomp)
  
  
  

  bow <- read.table(paste("benchmark", prname , "mydata-idf-BoW-matrix.csv", sep="/"), sep=",", row.names = 1, header = TRUE, check.names = FALSE)
  space1 = lsa(t(bow), dims)
  bow <- as.textmatrix(space1)
 
  mydata$lexsim <- cosine(bow)
  mydata$lexsim <- remove.diagonal(mydata$lexsim)
  
  mydata$lexsim <- mydata$lexsim[intersect(names(priori.decomp),rownames(mydata$lexsim)), intersect(names(priori.decomp),colnames(mydata$lexsim))]
  
  apply.hill.climb(mydata, rownames(mydata$lexsim), eval.func = Cpp.LQ.evaluator, neighbors=neighbors, 
                   size = size, ncpus=ncpus, priori.decomp)
  
  
}


apply.hill.climb <- function(mydata, names, eval.func, neighbors, size, ncpus, decomp) {
  
  
  clusters <- hill.climbing.search(mydata, names, eval.func, ncpus, size, neighbors)
  
  N = length(decomp)
  
  mojo <- compute.MoJo(clusters$groups, decomp)
  
  mojosim <- sapply(mojo, function(m) 1 - (m/N))
  
  print("printing mojosim")
  print(mojosim)
  
}
