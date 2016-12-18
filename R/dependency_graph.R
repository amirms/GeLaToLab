# setwd("C:/Users/AmirM/Documents/workspace/eclipse/org.servicifi.gelato.clustering/csv/DependencyGraph/test")
# filename <- "EditPane.java.csv"
# filename2 <- "ActionContext.java.csv"
# filename3 <- "BeanShell.java.csv"

projects <- list("jdom-2.0.5", "apache-ant-1.9.3", "hadoop-0.20.2", "apache-log4j-1.2.17", "jedit-5.1.0",
                 "jfreechart-1.2.0", "jhotdraw-7.0.6", "junit-4.12","weka-3.6.11", "eclipse-jdt-core-3.8")
#
#  ,"weka-3.6.11", "eclipse-jdt-core-3.8"
# "jdom-2.0.5" Type done, 
# 

run_random_walk_clustering <- function(projects){
  library(foreach)
  library(doParallel)
  library(GeLaToLab)
  
  setwd("~/workspace")
  
  # no_cores <- detectCores() - 1
  no_cores <- 3
  cl<-makeCluster(no_cores)
  registerDoParallel(cl)
  
  foreach(i = 1:length(projects)) %dopar% {
    GeLaToLab::compute.all.sim(projects[[i]])
  }
  
  stopCluster(cl)
}

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
compute.random.walk.sim <- function(prname, lex.fun= compute_normalized_LCS, alpha = 0.5, dirname="org", files=c(), min.nchar=4){
  library(foreach)
  library(doParallel)
  library(GeLaToLab)
  require(compiler)
  
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

    #FIXME can be optimized to just read module names
    dependency = read.graph.adjacency(filenames[i])
    stopifnot(rownames(dependency) == colnames(dependency))
    dependencies[[filename]] <- dependency
  }

  #eliminate small packages based on names(dependencies)
  dependencies <- dependencies[eliminate_small_packages(names(dependencies))]
  dependencies <- dependencies[load_module_names(prname)]

  # a sample of 10%
  # dependencies <- dependencies[sample(1:length(dependencies), ceiling(0.05 * length(dependencies)))]
  
  # Normalize the weighted adjacency matrix by outdegree
  dependencies <- lapply(dependencies, function(dep) normalize_adjacency_matrix(dep))
  
  alpha = estimate_alpha(dependencies)
  
  r <- matrix(0, nrow=length(dependencies), ncol= length(dependencies))
  cmp_compare_rw_graphs <- cmpfun(compare_rw_graphs)

  for(i in 1:length(dependencies)) {
    print(i)
    a1 <- dependencies[[i]]
    if (!is.null(min.nchar)){
      a1.names.indices <- which(nchar(rownames(a1)) > min.nchar)
      a1 <- a1[a1.names.indices, a1.names.indices]
      if (is.null(dim(a1)) || any(dim(a1) == c(0,0)))
        next
    }

    for (j in i: length(dependencies)){
      if (i != j) {
        print(j)
        
        a2 <- dependencies[[j]]
        
        if (!is.null(min.nchar)){
          a2.names.indices <- which(nchar(rownames(a2)) > min.nchar)
          a2 <- a2[a2.names.indices, a2.names.indices]
          if (is.null(dim(a2)) || any(dim(a2) == c(0,0)))
            next
        }
      
        r[i,j] <- cmp_compare_rw_graphs(a1, dependencies[[j]], alpha=alpha, lex.fun=lex.fun)
      }
    }
  }
  
  r <- fill_lower_diagonal(r)
  
  dimnames(r) <- list(names(dependencies), names(dependencies))
    
  return (r)
}



compute.all.sim <- function(prname){
  require(igraph)
  require(GeLaToLab)
  setwd("~/workspace")
  # Read the authoritative decomposition
  # decomposition <- read.csv(paste("benchmark", prname ,"decomposition.csv", sep="/"), sep=",",  header = TRUE)
  # priori.decomp <- decomposition$x
  # names(priori.decomp) <- decomposition$X
  # priori.decomp <- normalizeVector(priori.decomp) 
  # 
  # #get the source code units with small packages eliminated
  # priori.docs <- get_sample_docs(prname, priori.decomp, size=0.5)
  # priori.decomp <- priori.decomp[priori.docs]
  
  kernel <- compute.random.walk.sim(prname)
  

  kernel <- kernel[order(rownames(kernel)), order(colnames(kernel))]
  if (max(kernel) > 1)
    stop("wrong similarity matrix!")
  
  #compute distance from kernel
  myDist <- 1 - kernel
  myDist <- as.dist(myDist)
  
  # pinned it to complete linkage
  clusters <- hclust(myDist, method = 'complete')
  priori.decomp <- build.dendrogam(rownames(kernel))
  
  clusters.tree <- ape::as.phylo(clusters)
  priori.tree <- priori.decomp$tree
  path.difference <- phangorn::path.dist(clusters.tree, priori.tree, check.labels = T)
  
  clusters.dend <- as.dendrogram(clusters)
  priori.dend <- priori.decomp$dend
  
  baker <- cor_bakers_gamma.dendrogram(priori.dend, clusters.dend)
  
  treeDistance = compute_tree_edit_distance_for_hc(clusters, priori.decomp$graph)
  
  mojosim.ks <- MoJo.sim.k(priori.dend, clusters.dend)
  mojosim.k <- mean(unlist(lapply(mojosim.ks, function(mj) mj[1])))
  
  # return(list(baker=baker, cophcor=cophcor, Bk=Bk, diff=path.difference, mojosim = mojosim.k))
  results <- list(baker=baker, cophcor=0, Bk=0, diff=path.difference, mojosim = mojosim.k, treeDistance = treeDistance)

  setwd("~/workspace")
  dir.create(file.path(getwd(), paste("benchmark", prname, "Results/ContextModel", sep="/")), showWarnings = FALSE)
  
  print_clustering_results(prname, results, txt.file = paste("Results/ContextModel/DG-RW-Results.txt", sep=""))
  
  return (results) 
}
