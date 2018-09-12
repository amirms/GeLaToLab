compute_path_difference <- function(clusters, graph, labels){
  require(igraph)
    tree1 = convert_hclust_to_adjacency_tree(clusters)

d1 <- igraph::distances(graph_from_adjacency_matrix(tree1$adj))
d2 <-igraph::distances(graph)

d1 <- d1[labels, labels]
d2 <- d2[labels, labels]

sqrt(sum((d1-d2)^2)/2)

}