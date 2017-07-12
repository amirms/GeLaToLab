## Inspired by https://github.com/timtadh/zhang-shasha/blob/master/zss/compare.py

tree_edit_distance_test1 <- function(){
  #TEST DATA REMOVE LATER
  m <- matrix(0, nrow=6,ncol=6)
  rownames(m) <- c('a', 'b', 'c', 'd', 'e', 'f')
  colnames(m) <- rownames(m)
  m['f', 'd'] <- 1
  m['f', 'e'] <- 1
  m['d', 'a'] <- 1
  m['d', 'c'] <- 1
  m['c', 'b'] <- 1
  
  children <- function(node) {
    if (any(m[node, ]>0))
      return(names(which(m[node,]>0)))
    return(vector(mode="numeric", length=0))
  }
  
  m2 <- matrix(0, nrow=6,ncol=6)
  rownames(m2) <- c('a', 'b', 'c', 'd', 'e', 'f')
  colnames(m2) <- rownames(m2)
  m2['f', 'c'] <- 1
  m2['f', 'e'] <- 1
  m2['c', 'd'] <- 1
  m2['d', 'a'] <- 1
  m2['d', 'b'] <- 1
  
  children2 <- function(node) {
    if (any(m2[node, ]>0))
      return(names(which(m2[node,]>0)))
    return(vector(mode="numeric", length=0))
  }
  
  insert_cost <- function(x) {10*strdist('', x)}
  remove_cost <- function(x) {10*strdist(x,'')}
  update_cost <- function(x, y) {return(weird_dist(x,y))}

  
  weird_dist <- function(A, B){
    return(10*strdist(A, B))
  }
}

strdist <- function(a, b){
  if (a == b)
    return(0)
  
  return(1)
}

tree_edit_distance_test2 <- function(){
  m <- matrix(0, nrow=6,ncol=6)
  rownames(m) <- c('a', 'h', 'c', 'l', 'e', 'f')
  colnames(m) <- rownames(m)
  m['f', 'a'] <- 1
  m['a', 'h'] <- 1
  m['a', 'c'] <- 1
  m['c', 'l'] <- 1
  m['f', 'e'] <- 1
  
  children <- function(node) {
    if (any(m[node, ]>0))
      return(names(which(m[node,]>0)))
    return(vector(mode="numeric", length=0))
  }
  
  m2 <- matrix(0, nrow=6,ncol=6)
  rownames(m2) <- c('a', 'b', 'c', 'd', 'e', 'f')
  colnames(m2) <- rownames(m2)
  m2['f', 'a'] <- 1
  m2['a', 'd'] <- 1
  m2['a', 'c'] <- 1
  m2['c', 'b'] <- 1
  m2['f', 'e'] <- 1
  
  children2 <- function(node) {
    if (any(m2[node, ]>0))
      return(sort(names(which(m2[node,]>0))))
    return(vector(mode="numeric", length=0))
  }
  
  insert_cost <- function(x) {strdist('', x)}
  remove_cost <- function(x) {strdist(x,'')}
  update_cost <- function(x, y) {return(strdist(x,y))}
}