# 'clusters' is an instanceof 'hclust'
convert_hclust_to_adjacency_tree <- function(clusters) {
  require(igraph)
  require(ggraph)
  
  g <- den_to_igraph(clusters)
  
  names <- V(g)$name
  labels <- V(g)$label
  label_indices <- which(labels != "")
  names[label_indices] <- labels[label_indices]
  ay <- set.vertex.attribute(g, "name", value=names)
  
  adj <- as.matrix(as_adjacency_matrix(ay))

  #root is the one without incoming edges
  list(adj= as.matrix(as_adjacency_matrix(ay)), root= names[which(colSums(adj) == 0)])
}

convert_graph_to_adjacency_tree <- function(graph) {
  require(igraph)
  
  adj = as.matrix(as_adjacency_matrix(graph))
  names <- colnames(adj)
  adj[which(adj > 0)] <- 1

  list(adj= adj, root= names[which(colSums(adj) == 0)])
}

get_children <- function(node, m) {
  if (any(m[node, ]>0))
    return(names(which(m[node,]>0)))
  
  return(vector(mode="numeric", length=0))
}

compute_tree_edit_distance_for_hc <- function(clusters, graph){
  tree1 = convert_hclust_to_adjacency_tree(clusters)
  tree2 = convert_graph_to_adjacency_tree(graph)
  
  tree1$adj <- tree1$adj[order(rownames(tree1$adj)), order(colnames(tree1$adj))]
  tree2$adj <- tree2$adj[order(rownames(tree2$adj)), order(colnames(tree2$adj))]
  
  compute_tree_edit_distance(tree1, tree2)
}

#Input: tree = list(adjacency, root)
compute_tree_edit_distance <- function(tree1, tree2) {
  get_children1 <- function(node) {get_children(node, tree1$adj)}
  get_children2 <- function(node) {get_children(node, tree2$adj)}
  
  insert_cost <- function(x) {1}
  remove_cost <- function(x) {1}
  update_cost <- function(x, y) {
    if (x == y) {
       return(0)
    }
    
    # if (grepl('.java$', x) || grepl('.java$', y)) {
    #    return(Inf)
    # }

    if((length(get_children1(x))==0) || (length(get_children2(y))==0) ){
      return(Inf)
    }
    
    return(0)
  }
  
  tree_distance(tree1$root, tree2$root, get_children1, get_children2, insert_cost = insert_cost, remove_cost = remove_cost, update_cost = update_cost)
}