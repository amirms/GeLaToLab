build.dendrogram.from.project <- function(prname, dirname="org") {
  setwd("~/workspace")
  
  setwd(paste("benchmark", prname, sep="/"))
  
  files = list()
  
  pattern <- "*java"
  filenames = list.files(path = dirname, pattern = pattern, all.files = FALSE,
                         full.names = TRUE, recursive = TRUE,
                         ignore.case = TRUE, include.dirs = TRUE, no.. = TRUE)

  return(build.dendrogam(filenames))
}

build.dendrogam <- function(filenames){
  require(igraph)
  require(data.tree)
  require(dendextend)
  
  paths <- lapply(filenames, function(filename) {
    g <- grep("/", strsplit(filename, "")[[1]])
    lastChar <- g[length(g)]
    substr(filename, 1, lastChar-1)
  })
  paths <- unique(paths)
  
  vertices = c()
  edges = list()
  
  for (i in 1:length(paths)){
    path <- paths[[i]]
    g <- grep("/", strsplit(path, "")[[1]])

    lastSegment = NULL
    
    for (j in 1:length(g)){
    
      lastChar <- g[j]
      directory <- substr(path, 1, lastChar-1)
      
      if (!is.null(lastSegment)) {
        edges[[length(edges) + 1]] <- c(lastSegment, directory)
      }
      
      vertices <- c(vertices, directory)
      lastSegment <- directory
    }
    vertices <- c(vertices, path)
    
    if (!is.null(lastSegment)) {
      edges[[length(edges) + 1]] <- c(lastSegment, path)
    }

  }
  
  for (i in 1:length(filenames)){
    filename <- filenames[[i]]
    g <- grep("/", strsplit(filename, "")[[1]])
    
    lastChar <- g[length(g)]
    directory <- substr(filename, 1, lastChar-1)
      
    vertices <- c(vertices, filename)
    
    edges[[length(edges) + 1]] <- c(directory, filename)
    
  }
  
  vertices <- unique(vertices)
  
  el <- matrix( unlist(edges), nc = 2, byrow = TRUE)
  g <- graph_from_edgelist(el)
  df <- as_data_frame(g)
  df <- df[!duplicated(df), ]
  
  tree <- convert.edge.list.to.tree(df, filenames, vertices)
  
  dt <- FromDataFrameNetwork(df)
  # class(dt) <- "hclust"
  
  # tip.labels <- filenames
  # all.labels <- vertices
  
  dend <- as.dendrogram(dt)
  
  return(list(dend=dend, tree=tree))
}

convert.edge.list.to.tree <- function(edges, tip.labels, all.labels) {

  
  node.labels <- setdiff(all.labels, tip.labels)
  
  normalized_edges = matrix(0, ncol = 2, nrow = dim(edges)[1])
  
  n<-c((length(tip.labels)+1):length(all.labels),1:length(tip.labels))
  
  for (i in 1:dim(edges)[1])
    for (j in 1:dim(edges)[2]) {
      normalized_edges[i,j] <- n[which(all.labels %in% edges[i,j])]
      
    }
  
  
  #Now just create the phylo object
  tree<-list(
    edge=normalized_edges,
    tip.label= tip.labels,
    Nnode=length(node.labels),
    # edge.length=rep(1, dim(normalized_edges)[1]),
    node.label = node.labels)
  
  class(tree)<-"phylo"
  tree <- ape::reorder.phylo(tree)
  
  return(tree)
  
}


sim2dist <- function(S, type = "diff"){
  d <- switch(type,
        diff = 1-S,
        sqrt = sqrt(1-S),
       log = -log(S),
       inv = (1/(S+1))
       )
  d[is.infinite(d)] <- 0 
  as.dist(d)
  
}
dist.types <- c("diff", "sqrt", "inv")

#so far, (inv, 1), (diff, 2), (sqrt, 3) => diff = 1-S is the best

# compute.cophenetic.correlation <- function(dend1, dend2) {
#   cor(cophenetic(dend1), cophenetic(dend2))
# }

compute_hierarchical_clustering <- function(semantic, myBoF){
  require(dendextend)
  
  #SVD to compute to USU^T
  USUt <- svd(semantic)
  S <- USUt$u %*% diag(sqrt(USUt$d))
  
  #Compute cosine similarity
  Phi_d <- apply_tf_idf(myBoF) %*% S
  
  dimnames(Phi_d) <- dimnames(myBoF)
  Phi_d <- Phi_d[order(rownames(Phi_d)),]
  
  print("Dimensions of Phi_d before Cleansing")
  print(dim(Phi_d))
  
  #Remove empty rows
  Phi_d <- Phi_d[ apply(Phi_d!=0, 1, any), , drop=FALSE] 
  # #Remove duplicated rows
  # Phi_d <- Phi_d[!duplicated(Phi_d),]
  
  
  print("Dimensions of Phi_d after Cleansing")
  print(dim(Phi_d))
  
  kernel <- compute_cosine_kernel(Phi_d)
  # kernel <- kernel[order(rownames(kernel)), order(colnames(kernel))]
  
  if (max(kernel) > 1)
    stop("wrong similarity matrix!")

  #compute distance from kernel
  myDist <- 1 - kernel
  myDist <- as.dist(myDist)
  
  # pinned it to complete linkage
  clusters <- hclust(myDist, method = 'complete')
  priori.decomp <- build.dendrogam(rownames(Phi_d))
  
  clusters.tree <- ape::as.phylo(clusters)
  priori.tree <- priori.decomp$tree
  path.difference <- phangorn::treedist(clusters.tree, priori.tree, check.labels = T)[2]
  
  clusters.dend <- as.dendrogram(clusters)
  priori.dend <- priori.decomp$dend
  
  baker <- cor_bakers_gamma.dendrogram(priori.dend, clusters.dend)
  cophcor <- dendextend::cor_cophenetic(priori.dend, clusters.dend)
  Bks <- Bk2(priori.dend, clusters.dend)
  Bk <- mean(unlist(lapply(Bks, function(b) b[1])))
  
  return(list(baker=baker, cophcor=cophcor, Bk=Bk, diff=path.difference))
}

Bk2 <- function(tree1, tree2, include_EV = TRUE, warn = dendextend_options("warn")) 
{
  require(dendextend)

  if (warn) {
    tree1_labels <- labels(tree1)
    tree2_labels <- labels(tree2)
    length_tree1_labels <- length(tree1_labels)
    length_tree2_labels <- length(tree2_labels)
    if (length_tree1_labels != length_tree2_labels) 
      stop("The two clusters don't have the same number of items!")
    if (!all(sort(tree1_labels) == sort(tree2_labels))) 
      stop("Your trees are having leaves with different names - please correct it in order to use this function")
  }
    #find the minimum of the two possible height
    tree1_heights_per_k <- heights_per_k.dendrogram(tree1)
    tree2_heights_per_k <- heights_per_k.dendrogram(tree2)
    
    if (length(tree1_heights_per_k) > length(tree2_heights_per_k))
        dend_heights_per_k <- tree2_heights_per_k
    else 
        dend_heights_per_k <- tree1_heights_per_k
    
    # k <- everything except 1 and nleaves
    ks_to_use <- names(dend_heights_per_k[which(!(names(dend_heights_per_k) %in% c(1,nleaves(tree1))))])

  cutree_tree1 <- lapply(ks_to_use, function(k) cutree.k.dendrogram(tree1, k = as.integer(k)))
  cutree_tree2 <- lapply(ks_to_use, function(k) dendextend::cutree(tree2, k = as.integer(k)))
  
  cutree_tree1 <- lapply(cutree_tree1, function(t) t[order(t)])
  cutree_tree1 <- lapply(cutree_tree1, function(t) normalizeVector(t))
  
  if (length(ks_to_use) == 1) {
    cutree_tree1 <- as.matrix(cutree_tree1)
    cutree_tree2 <- as.matrix(cutree_tree2)
  }
  n_ks <- length(cutree_tree1)
  Bk_for_each_k <- function(i_k) {
    FM_index(cutree_tree1[[i_k]], cutree_tree2[[i_k]], assume_sorted_vectors = FALSE, 
             include_EV = include_EV, warn = warn)
  }
  the_Bks <- lapply(seq_len(n_ks), Bk_for_each_k)
  names(the_Bks) <- ks_to_use
  return(the_Bks)
}

labels <- function(dend) {
  l <- dend %>% get_nodes_attr("label")
  l[which(!is.na(l))]
}

# methods <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
#Fight is between ward.D2, complete, and average
methods <- c("ward.D2", "complete", "average")

# So far 'complete'

test_best_method <- function(myDist, method="complete"){
  clusters <- hclust(myDist, method = method)
  priori.decomp <- build.dendrogam(rownames(Phi_d))
  
  clusters.tree <- ape::as.phylo(clusters)
  priori.tree <- priori.decomp$tree
  diff <- phangorn::treedist(clusters.tree, priori.tree, check.labels = T)
  
  clusters.dend <- as.dendrogram(clusters)
  priori.dend <- priori.decomp$dend
  
  baker <- cor_bakers_gamma.dendrogram(priori.dend, clusters.dend)
  cophcor <- dendextend::cor_cophenetic(priori.dend, clusters.dend)
  Bks <- Bk2(priori.dend, clusters.dend)
  Bk <- mean(unlist(lapply(Bks, function(b) b[1])))
  
  return(list(baker=baker, cophcor=cophcor, Bk=Bk, diff=diff))
}

# r2 <- lapply(methods, function(m) test_best_method(myDist, m))
# baker.order <- order(unlist(lapply(r, function(x) x$baker)))
# cophcor.order <- order(unlist(lapply(r, function(x) x$cophcor)))
# unlist(lapply(r, function(x) x$baker))