# setwd("C:/Users/AmirM/Documents/workspace/eclipse/org.servicifi.gelato.clustering/csv/DependencyGraph/test")
# filename <- "EditPane.java.csv"
# filename2 <- "ActionContext.java.csv"
# filename3 <- "BeanShell.java.csv"

PATH = "C:/Users/AmirM/Documents/workspace/eclipse/org.servicifi.gelato.clustering/csv"
PRNAME = "test-project"

visualize.semantic.network <- function(path = PATH, prname=PRNAME) {
  require(igraph)
  setwd(paste(path, "SemanticNetwork", sep="/"))
  
  filename <- paste(paste(prname, "SN", sep="-"), "csv", sep=".")
  
  adj <- read.csv(filename, header = TRUE, sep = ",", quote = "\"",
                  dec = ".", fill = TRUE)
  
  rownames(adj) <- adj[,1]
  adj <- adj[,-1]
  
  colnames(adj) <- rownames(adj)
  
  adj <- as.matrix(adj)

  g <- graph_from_adjacency_matrix(adj, mode = "directed", weighted = T, diag = TRUE)
  
  plot(g)
  # 
  # return(g)
}

#filename1 = "org/employee/Employee.java.csv"
#folder1 = ""

visualize.semantic.dependency.graph <- function(path = PATH, prname=PRNAME, filename) {
  require(igraph)
  
  setwd(paste(path, "DependencyGraph", prname, sep="/"))
  
  adj <- read.csv(filename, header = TRUE, sep = ",", quote = "\"",
                  dec = ".", fill = TRUE)
  
  rownames(adj) <- adj[,1]
  adj <- adj[,-1]
  
  colnames(adj) <- rownames(adj)
  
  adj <- as.matrix(adj)
  
  g <- graph_from_adjacency_matrix(adj, mode = "directed", weighted = T, diag = TRUE)
  
  plot(g)
  # 
  # return(g)
}