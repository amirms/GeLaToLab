
visualize <- function(m, grp) {
  require(igraph)
  
  ig <- graph.adjacency(m, mode="undirected", weighted=TRUE)
  
  grp.names = lapply(1:max(grp), function(x) which(grp==x) )
  
  #plot(ig, mark.groups=list(c("DR", "2C", "gt"), c("sv", "ME", "I2"), c("zU", "W3", "uz", "Yj")))
plot(ig, mark.groups=grp.names, layout=layout.auto)  
#coords <- layout.fruchterman.reingold(ig, dim=3)
#rglplot(ig, layout=coords)

#colbar <- rainbow(max(grp)+1)
#V(ig)$color <- colbar[grp+1]

#plot(ig, layout=layout.spring)

}
