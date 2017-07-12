## Inspired by https://github.com/timtadh/zhang-shasha/blob/master/zss/compare.py

#TEST DATA REMOVE LATER
# m <- matrix(0, nrow=3,ncol=3)
# rownames(m) <- c('a', 'b', 'c')
# colnames(m) <- rownames(m)
# m['a', 'b'] <- 1
# m['a', 'c'] <- 1
# 
# children <- function(node) {
#   if (any(m[node, ]>0))
#     return(names(which(m[node,]>0)))
#  return(vector(mode="numeric", length=0))
# }
# 
# m2 <- matrix(0, nrow=4,ncol=4)
# rownames(m2) <- c('a', 'b', 'c', 'd')
# colnames(m2) <- rownames(m2)
# m2['a', 'b'] <- 1
# m2['b', 'c'] <- 1
# m2['a', 'd'] <- 1
# # m2['a', 'e'] <- 1
# 
# children2 <- function(node) {
#   if (any(m2[node, ]>0))
#     return(names(which(m2[node,]>0)))
#   return(vector(mode="numeric", length=0))
# }
# 
# insert_cost <- function(x) {1}
# remove_cost <- function(x) {1}
# update_cost <- function(x, y) {1}

tree_distance = function(A, B, get_children_A, get_children_B, insert_cost, remove_cost, update_cost){
    "Computes the exact tree edit distance between trees A and B with a
  richer API than :py:func:`zss.simple_distance`.
  Use this function if either of these things are true:
  * The cost to insert a node is **not** equivalent to the cost of changing
  an empty node to have the new node's label
  * The cost to remove a node is **not** equivalent to the cost of changing
  it to a node with an empty label
  Otherwise, use :py:func:`zss.simple_distance`.
  :param A: The root of a tree.
  :param B: The root of a tree.
  :param get_children:
    A function ``get_children(node) == [node children]``.  Defaults to
  :py:func:`zss.Node.get_children`.
  :param insert_cost:
    A function ``insert_cost(node) == cost to insert node >= 0``.
  :param remove_cost:
    A function ``remove_cost(node) == cost to remove node >= 0``.
  :param update_cost:
    A function ``update_cost(a, b) == cost to change a into b >= 0``.
  :return: An integer distance [0, inf+)
  "
  A = AnnotatedTree(A, get_children_A)
  B = AnnotatedTree(B, get_children_B)
  treedists = matrix(0, nrow = length(A$nodes), ncol = length(B$nodes))

  treedist = function(i, j, treedists) {
    Al = A$lmds
    Bl = B$lmds
    An = A$nodes
    Bn = B$nodes
    
    m = i - Al[[i]] + 2
    n = j - Bl[[j]] + 2
    fd = matrix(0, nrow= m,ncol=n)
    
    ioff = Al[[i]] - 2
    joff = Bl[[j]] - 2
    
    for (x in 2:m){ # delta(l(i1)..i, Theta) = delta(l(1i)..1-1, Theta) + gamma(v -> lambda)
      # fd[x][0] = fd[x-1][0] + remove_cost(An[x+ioff])
      fd[x,1] = fd[x-1,1] + remove_cost(An[[x+ioff]])
    }
    for (y in 2:n){ # delta(Theta, l(j1)..j) = delta(Theta, l(j1)..j-1) + gamma(lambda -> w)
      # fd[0][y] = fd[0][y-1] + insert_cost(Bn[y+joff])
      fd[1,y] = fd[1,y-1] + insert_cost(Bn[[y+joff]])
    }
    
    for (x in 2:m){ ## the plus one is for the xrange impl
      for (y in 2:n){
        # only need to check if x is an ancestor of i
        # and y is an ancestor of j
        
        # print(paste("x+ioff: ", x+ioff))
        # print(paste("y+joff: ", y+joff))
        # 
        # print(paste("Al[[i]]: ", Al[[i]]))
        # print(paste("Al[[x+ioff]]: ", Al[[x+ioff]]))
        # print(paste("Bl[[j]]: ", Bl[[j]]))
        # print(paste("Bl[[y+joff]]: ", Bl[[y+joff]]))
        # print(paste("Is it equal: ", ((Al[[i]] == Al[[x+ioff]]) && (Bl[[j]] == Bl[[y+joff]]))))
        # 
        # print(paste("An[[x+ioff]]: ", An[[x+ioff]]))
        # print(paste("Bn[[y+joff]]: ", Bn[[y+joff]]))
        if ((Al[[i]] == Al[[x+ioff]]) && (Bl[[j]] == Bl[[y+joff]])) {
          #                   +-
          #                   | delta(l(i1)..i-1, l(j1)..j) + gamma(v -> lambda)
          # delta(F1 , F2 ) = min-+ delta(l(i1)..i , l(j1)..j-1) + gamma(lambda -> w)
          #                   | delta(l(i1)..i-1, l(j1)..j-1) + gamma(v -> w)
          #                   +-
          # fd[x][y] = min(
          # fd[x-1][y] + remove_cost(An[x+ioff]),
          # fd[x][y-1] + insert_cost(Bn[y+joff]),
          # fd[x-1][y-1] + update_cost(An[x+ioff], Bn[y+joff])
          # )
          fd[x, y] = min(
            fd[x-1, y] + remove_cost(An[[x+ioff]]),
            fd[x, y-1] + insert_cost(Bn[[y+joff]]),
            fd[x-1,y-1] + update_cost(An[[x+ioff]], Bn[[y+joff]])
          )
          # print(paste("x+ioff: ", x+ioff))
          # print(paste("y+joff: ", y+joff))
          # treedists[x+ioff][y+joff] = fd[x][y] FIXME
          treedists[x+ioff, y+joff] = fd[x, y]
        } else {
          #                   +-
          #                   | delta(l(i1)..i-1, l(j1)..j) + gamma(v -> lambda)
          # delta(F1 , F2 ) = min-+ delta(l(i1)..i , l(j1)..j-1) + gamma(lambda -> w)
          #                   | delta(l(i1)..l(i)-1, l(j1)..l(j)-1)
          #                   |                     + treedist(i1,j1)
          #                   +-
          # p = Al[x+ioff]-1-ioff
          # q = Bl[y+joff]-1-joff
          p = Al[[x+ioff]]-1-ioff
          q = Bl[[y+joff]]-1-joff
          # print(paste("dim(fd): ", dim(fd)))
          # print(paste("p: ", p))
          # print(paste("q: ", q))
          # print(paste("fd[p, q]: ", fd[p, q]))
          #print (p, q), (len(fd), len(fd[0]))
          # fd[x][y] = min(
          #   fd[x-1][y] + remove_cost(An[x+ioff]),
          #   fd[x][y-1] + insert_cost(Bn[y+joff]),
          #   fd[p][q] + treedists[x+ioff][y+joff]
          # )
          fd[x, y] = min(
            fd[x-1, y] + remove_cost(An[[x+ioff]]),
            fd[x, y-1] + insert_cost(Bn[[y+joff]]),
            fd[p, q] + treedists[x+ioff, y+joff] #FIXME
          )
        }
      }
    }
    
    return(treedists)
  }

  for (i in A$keyroots)
    for (j in B$keyroots) {
      treedists <- treedist(i,j, treedists)
    }
  
  # return (treedists[-1][-1])
  # return(treedists)
  return (treedists[dim(treedists)[1], dim(treedists)[2]])
}

AnnotatedTree <- function(root, get_children){
  require(rstackdeque)
  
  me = list(
    get_children= get_children,
    root= root,
    nodes=list(),  # a post-order enumeration of the nodes in the tree
    ids=list(), # a matching list of ids
    lmds=list(), # left most descendents
    keyroots=NULL
  )
  
  # k and k' are nodes specified in the post-order enumeration.
  # keyroots = {k | there exists no k'>k such that lmd(k) == lmd(k')}
  # see paper for more on keyroots
  
  stack = rstack()
  pstack = rstack()
  # stack.append(c(root, deque()))
  stack <- insert_top(stack, list(root, rdeque()))
  j = 1
  while (length(stack) > 0){
    #(n, anc) = stack.pop();
    x = peek_top(stack)
    n = x[[1]]
    anc = x[[2]]
    stack <- without_top(stack)
    
    nid = j
    
    for (c in me$get_children(n)){
      # print(paste("child: ", c))
      a = as.rdeque(anc)

      # a.appendleft(nid)
      a <- insert_back(a,nid)
      
      # stack.append(c(c, a))
      stack <- insert_top(stack, list(c, a))
    }
    
    # pstack.append(c(c(n, nid), anc))
    pstack <- insert_top(pstack, list(list(n, nid), anc))
    j <- j+ 1
  }
  # lmds = dict()
  # keyroots = dict()
  lmds = list()
  keyroots = list()
  i = 1
  while (length(pstack) > 0) {
    # (n, nid), anc = pstack.pop()
    x = peek_top(pstack)
    n = x[[1]][[1]]
    nid = x[[1]][[2]]
    anc = x[[2]]
    pstack <- without_top(pstack)
    
    me$nodes <- append(me$nodes, n)
    me$ids <- append(me$ids, nid)
    
    #print n.label, [a.label for a in anc]
    #if not me.get_children(n)
    if (length(me$get_children(n)) == 0) {
      lmd = i

      # for (a in anc) {
      while(length(anc) > 0){
        a <- peek_back(anc)
        anc <- without_back(anc)
        
        # print(paste("a: ", a))
        #if a not in keys of lmds
        if (is.null(lmds[a][[1]])){
             lmds[a] = i
        }
        else break
      }
    } else {
        # try: lmd = lmds[nid]
      lmd = lmds[[nid]]
      # except:
      #   import pdb
    }
    # pdb.set_trace()
    me$lmds <- append(me$lmds, lmd)
    keyroots[lmd] = i
    i <- i + 1
  }
  
  values <- unlist(lapply(keyroots, function(keyValue) keyValue[[1]]))
  
  me$keyroots = sort(values)
  
  ## Set the name for the class
  class(me) <- append(class(me),"AnnotatedTree")
  return(me)
}