recover.DAG = function(K, rootNodeIdx, isTree=TRUE) {
  require(igraph)
  #stopifnot(isSymmetric(K), "K is not symmetric");
  upper_sim <- K[upper.tri(K, diag = FALSE)]
  #K[lower.tri(K, diag = T)] <- 0
  sorted_sim <- sort(K, decreasing = T)
  
  #create an empty graph to check to cycles
  g <- make_empty_graph(n = dim(K)[1], directed = TRUE)
  #the corresponding adjacency matrix
  adj <- matrix(0, nrow=dim(K)[1], ncol= dim(K)[2])
  
  i =1;
  while(i <= length(sorted_sim)) {
    current_max_sim = sorted_sim[i];
    indices <- which(K == current_max_sim, arr.ind = TRUE);
    
    for (j in 1:dim(indices)[1]){
      i<- i+1;
      
      rowNum <- indices[j,1]
      colNum <- indices[j,2]
      
      print(paste("connecting row: ", rowNum, " to ", colNum))
      
      #Is this the root Node, then no outgoing edges
      if (rowNum == rootNodeIdx) {
        next;
      }
      
      # THINK if isTree == destination is a non-int erface
      if(isTree){
        # Does node[rowNum] already have a parent
        if (max(adj[rowNum, ]) > 0) {
          next;
        } 
      }
      
      #does adding (rowNum, colNUm) to g, create a cycle?
      g_clone = g
      g_clone <- add_edges(g_clone, c(rowNum, colNum));
      #if (is_dag(g_clone)) {
        #add the edge to g and cc
        g <- g_clone
        adj[rowNum, colNum] <- 1
      #}

    }
    
  }
  
  adj
  
}

recover.Tree = function(K, rootNodeIdx) {
  require(igraph)
  #stopifnot(isSymmetric(K), "K is not symmetric");
  upper_sim <- K[upper.tri(K, diag = FALSE)]
  #K[lower.tri(K, diag = T)] <- 0
  sorted_sim <- sort(K, decreasing = T)
  
  #create an empty graph to check to cycles
  g <- make_empty_graph(n = dim(K)[1], directed = FALSE)
  
  i =1;
  while(i <= length(sorted_sim)) {
    current_max_sim = sorted_sim[i];
    indices <- which(K == current_max_sim, arr.ind = TRUE);
    
    for (j in 1:dim(indices)[1]){
      i<- i+1;
      
      rowNum <- indices[j,1]
      colNum <- indices[j,2]
      
      print(paste("connecting row: ", rowNum, " to ", colNum))
  
      
      #does adding (rowNum, colNUm) to g, create a cycle?
      g_clone = g
      g_clone <- add_edges(g_clone, c(rowNum, colNum));
      if (is_dag(g_clone)) {
      #add the edge to g and cc
      g <- g_clone
      }
      
    }
    
  }
  
  g
  
} 
