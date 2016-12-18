# setwd("C:/Users/AmirM/Documents/workspace/eclipse/org.servicifi.gelato.clustering/csv/DependencyGraph/test")
# filename <- "EditPane.java.csv"
# filename2 <- "ActionContext.java.csv"
# filename3 <- "BeanShell.java.csv"

read.graph.adjacency <- function(filename) {
  require(igraph)
  adj <- read.csv(filename, header = TRUE, sep = ",", quote = "\"",
                  dec = ".", fill = TRUE)
  
  rownames(adj) <- adj[,1]
  adj <- adj[,-1]
  
  colnames(adj) <- rownames(adj)
  
  adj <- as.matrix(adj)
  
  return(adj)
}

estimate_alpha <- function(dependencies){
#   A rule of thumb for setting it is to take the largest power of 10
#   which is samller than 1/d^2, being d the largest degree in the 
#   dataset of graphs.
  
  d <- max(unlist(lapply(dependencies, function(dep) max(rowSums(dep)))))
  
  threshold <- 1/ d^2
  
  alpha <- 1
  alpha_10 <- alpha^10
  
  while(alpha_10 >= threshold && alpha > 0.01 )
  {
    alpha <- alpha - 0.01
    
    alpha_10 <- alpha^10
  }
  
  return(alpha)
}


#lex.fun <- compute_normalized_LCS
compute.random.walk.sim <- function(prname, lex.fun= compute_normalized_LCS, alpha = 0.5, dirname="org", files=c()){
  setwd("~/workspace")
  setwd(paste("benchmark", prname , "DG", sep="/"))
  
  
  
#   setwd("C:/Users/AmirM/Documents/workspace/eclipse/org.servicifi.gelato.clustering/csv/DependencyGraph/test")
  
  dependencies = list()

  pattern <- "*.csv"
  filenames = list.files(path = dirname, pattern = pattern, all.files = FALSE,
                         full.names = TRUE, recursive = TRUE,
                         ignore.case = TRUE, include.dirs = TRUE, no.. = TRUE)
  #   srcfiles.names = c(0)
  for (i in 1:length(filenames)) {
    filename <- substr(filenames[i], 1, regexpr(".csv", filenames[i]) -1)
    
#     print(filename)
    
    if ((length(files) > 0)  && !(filename %in% files ))
      next
    
    dependencies[[filename]] = read.graph.adjacency(filenames[i])
  }
  
  # Normalize the weighted adjacency matrix by outdegree
  dependencies <- lapply(dependencies, function(dep) normalize_adjacency_matrix(dep))
  
  alpha = estimate_alpha(dependencies)
  
  r <- matrix(0, nrow=length(dependencies), ncol= length(dependencies))
    
  for(i in 1:length(dependencies))
    for (j in i: length(dependencies)){
      if (i != j) {
      
        print(i)
        print(j)
        r[i,j] <- compare_rw_graphs(dependencies[[i]], dependencies[[j]], alpha=alpha, lex.fun=lex.fun)
      }
    }
  
  r <- fill_lower_diagonal(r)
  
  dimnames(r) <- list(filenames, filenames)
    
  return (r)
  }



compute.all.sim <- function(prname){
  require(igraph)
  require(GeLaToLab)
  setwd("~/workspace")
  # Read the authoritative decomposition
  decomposition <- read.csv(paste("benchmark", prname ,"decomposition.csv", sep="/"), sep=",",  header = TRUE)
  priori.decomp <- decomposition$x
  names(priori.decomp) <- decomposition$X
  priori.decomp <- normalizeVector(priori.decomp) 
  
  #get the source code units with small packages eliminated
  priori.docs <- get_sample_docs(prname, priori.decomp, size=0.5)
  priori.decomp <- priori.decomp[priori.docs]
  
  r <- compute.random.walk.sim(prname, files=names(priori.decomp))
  
  #Fix priori decomposition 
  modules <- intersect(names(priori.decomp), rownames(r))
  
  priori.decomp <- priori.decomp[modules]
#   r <- r[modules, modules]
  r <- r[order(rownames(r)), order(colnames(r))]

  priori.decomp <- priori.decomp[order(names(priori.decomp))]
  priori.decomp <- normalizeVector(priori.decomp)
  
  if(!all(rownames(r) == names(priori.decomp)))
    stop("names don't match!")
  
  #K number of clusters
  noc <- max(priori.decomp)
  
  
  L <- normalized.symmetric.laplacian(r)
  
  partition <- spectral.clustering(L, noc)
partition <- normalizeVector(partition)

  names(partition) <- colnames(r)
  
  f1.score <- compute.f1(partition, priori.decomp)  
  adjustedRI <- compute.AdjRI(partition, priori.decomp)
  mojosim <- compute.MoJoSim(partition, priori.decomp)
  
  results <- list(mojosim = mojosim, f1.score=f1.score, adjustedRI=adjustedRI)
  
  print_clustering_results(prname, results, "DG-RW-Results.txt")
  
  return (results) 
}
