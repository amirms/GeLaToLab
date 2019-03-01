plot_dependency_graph <- function() {
library(igraph)
library(readtext)

g <- igraph::make_empty_graph();

#read files
directory = "C:\\workspace\\org.servicifi.gelato.dependency\\myData\\test-project"
files = c("test.ee.txt", "train.ee.txt", "valid.ee.txt");

d = lapply(files, function(f) readtext(paste(directory, f, sep="\\")))
z <- unlist(lapply(d, function(x) strsplit(x$text, '\n') [[1]]))
z <- gsub("#", "$", z)

all <- read.table(text = z,  sep=";",
                  col.names=c("head", "tail", "relation"))

# all <- as.matrix(rbind(d[[1]], d[[2]], d[[3]]))

all <- as.matrix(all)
print(all)

identifiers <- unique(c(all[,1], all[,2]))

g <- make_empty_graph() %>%
  add_vertices(length(identifiers))   %>%
  set_vertex_attr("label", value = identifiers)

for (i in 1:dim(all)[1]){
  headIndex = which(identifiers == all[i,1])
  tailIndex = which(identifiers == all[i,2])
  
  g <- add_edges(g, c(headIndex, tailIndex), label=all[i,3])
}

tkplot(g, vertex.size=10,       vertex.color="green")
}