# Read in a dependency XMI file
# Input: filename
# Calls: XML package (available from CRAN)
#        strip.text()
# Output: vector of character strings, giving the words in order
read.doc <- function(filename) {
  
  fulltext <- readLines(filename, encoding="UTF-8")
  # ASSUMES: this should be a SINGLE character string
  text <- strip.text(fulltext) # Turn into a vector of strings
  return(text)
}

#spaces

# Exhaustive search to group together objects based on an objective function
# Input: matrix representing the relationships between objects
# an objective function to maximize (eval.fun)
# Output: a list of vectors representing grouped objects
#FIXME redo for a list of matrices, instead of cfg only
exhaustive.search <- function (mydata, names, eval.fun) 
{
  mydata1 = mydata$cfg
  
  require(blockmodeling)
  dims = sapply(mydata, function(d){dim(d)})
  
  #is it one or two
  if (!all(apply(dims, 1, function(r) unique(r))))
    stop("the matrices have different dimensions")
  
  if (!all(apply(dims, 2, function(c) unique(c))))
    stop("the dimensions do not match")
  
  #a trick
  len <- dims[1,1]
  if (len == 0) 
    stop("Empty dataset")
  
  eval.fun = match.fun(eval.fun)
  
  #Remove all self-dependencies
  #mydata <- remove.diagonal(mydata)
  
  
  best = list(result = -Inf, groups=rep(0, len))
  for (size in 1:dim(mydata1)[1]) {
      child_comb = nkpartitions(len, size, exact = TRUE, print = FALSE)
      for (i in 1:dim(child_comb)[1]) {
        clusters = child_comb[i,]
        names(clusters) <- names
        result = eval.fun(mydata, clusters)
        if (result > best$result) {
          best$result = result
          best$groups = clusters
        }
        
    }
    
  } 
  return(best)
}

#Input: an adjacency matrix
remove.diagonal <- function(d) {
  
  
  len = dim(d)[1]
  
  if (len == 0)
    stop("Empty dataset")
  
  for (i in 1:len) {
    d[i,i] <- 0
  }
  
  return(d)
}

Cpp.MQ.evaluator <- function(mydata, group) {
  
  if (is.null(mydata$cfg))
    stop("There is no adjacency matrix denoting the CFG")
  
  return(evaluateMQ(mydata$cfg, group))
  
}

Cpp.LQ.evaluator <- function(mydata, group) {
  
  if (is.null(mydata$lexsim))
    stop("There is no lexical similarity/dissimilarity difference for the system")
  
  return(evaluateLQ(mydata$lexsim, group))
  
}

#Input: named cluster vector
MQ.evaluator <- function(mydata, group) {
  
  #print(group)
  
  compute_cluster_factor <- function(m, n) {
    if (m > 0) return((2 * m) / ((2 * m) + n))
    return(0)
  }
  
  if (is.null(mydata$cfg))
    stop("There is no adjacency matrix denoting the CFG")
  
  #k is the number of clusters  
  k <- max(group)
  
  if (k <= 0)
    stop("wrong number of clusters")
    
  #cluster factors scores
  cluster.factor <- c(rep(0, k))
  
  grouped.nodes <- lapply(1:k, function(x) {group[group %in% x]})
  
  intra.edges <- sapply(grouped.nodes, function(x){ sum(mydata$cfg[names(x), names(x)])})
  
  all.edges <- sapply(grouped.nodes, function(x){ sum(mydata$cfg[names(x), ], mydata$cfg[, names(x)])})
  
  inter.edges <- all.edges - 2 * intra.edges
   
  cluster.factor <- mapply(compute_cluster_factor, intra.edges, inter.edges)
  
  return(sum(cluster.factor))

  }

#Lexical Quality
#Input: similarity matrix
# a dissimilarity matrix is derived
LQ.evaluator <- function(mydata, group) {
  
  compute_cluster_factor <- function(m, n) {
    if ((m + n) != 0) return((m - n) / abs(m + n))
      return(0)
  }
  
  if (is.null(mydata$lexsim))
    stop("There is no lexical similarity/dissimilarity difference for the system")

  #create dissimilarity matrix
  #lexdis <- 1 - mydata$lexsim
  #k is the number of clusters  
  k <- max(group)
  
  if (k <= 0)
    stop("wrong number of clusters")
  #cluster factors scores
  cluster.factor <- c(rep(0, k))
  
  grouped.nodes <- lapply(1:k, function(x) {group[group %in% x]})
  
  intra.lex.sim <- sapply(grouped.nodes, function(x){sum(mydata$lexsim[names(x), names(x)])})

  all.lex.sim <- sapply(grouped.nodes, function(x){sum(mydata$lexsim[names(x), ])})

  inter.lex.sim.dis <- all.lex.sim - intra.lex.sim
  
  cluster.factor <- mapply(compute_cluster_factor, intra.lex.sim, inter.lex.sim) 

  return(sum(cluster.factor))
  
}

make.random.string <- function(n=1, length=12)
{
  randomString <- c(1:n)                  # initialize vector
  for (i in 1:n)
  {
    randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
                                    length, replace=TRUE),
                             collapse="")
  }
  return(randomString)
}

apply.weighted.hybrid.evaluator <- function(ws, mydata, ncpus=8, size=24, neighbors=.6, tries=1) {
 return(sapply(ws, 
    function(w) { 
      evaluator <- function(d, g) Cpp.hybrid.linear.MQ.LQ.evaluator(d, g, w)
    
      best=list(result = -Inf, groups=rep(0, dim(mydata$cfg)[1]))
      
      for (i in 1:tries) {
        
        current = 
          hill.climbing.search(mydata, rownames(mydata$cfg), evaluator, ncpus=ncpus, size=size, neighbors=neighbors)
    
        if (current$result > best$result)
          best <- current
        
      }
      
      return(best)
      }
    ))
}

#Input: mydatas: a list of two matrices, one for CF and the other for Lexical representation
# weight: the lambda denoting the coefficient of the MQ evaluator
hybrid.linear.MQ.LQ.evaluator <- function(mydata, group, weight) {
  
 funclist = list(mq=MQ.evaluator(mydata, group), lq =LQ.evaluator(mydata, group))
 
 return(hybrid.linear.evaluator(funclist, c(weight, 1-weight)))
  
}

Cpp.hybrid.linear.MQ.LQ.evaluator <- function(mydata, group, weight) {
  
  funclist = list(mq=Cpp.MQ.evaluator(mydata, group), lq =Cpp.LQ.evaluator(mydata, group))
  
  return(hybrid.linear.evaluator(funclist, c(weight, 1-weight)))
  
}


#Input: eval.fun.list: a list of objective functions with one-to-one correspondence with mydata
# weights: a vector representing the weights for each objective function
hybrid.linear.evaluator <- function(eval.fun.list, weights) {
  
  if (sum(weights) != 1)
    stop("The sum of weights not equal to zero")
  
  #Enforce lazy evaluation
  results <- mapply(function(w,f){
    if (w==0)
      return(0)
    return(w*f)}, weights, eval.fun.list)
  
  return(sum(results))
  
}

#Input: neighbors: a fraction indicating the percentage of the neighbors that has to be evaluated

hill.climbing.search <- function(mydata, names, eval.fun, ncpus = 8, size=5, neighbors = .20) {
  
  require(snowfall)
  require(gelato)
  
  
  dimensions <- sapply(mydata, function(d){dim(d)})
  
  #is it one or two
  if (!all(apply(dimensions, 1, function(r) unique(r))))
      stop("the matrices have different dimensions")
      
  if (!all(apply(dimensions, 2, function(c) unique(c))))
    stop("the matrices are not square")

  #a trick
  len <- dimensions[1,1]
  if (len == 0) 
    stop("Empty dataset")
  
  eval.fun = match.fun(eval.fun)

  best=list(result = -Inf, groups=rep(0, len))

  #Make a population for k=number of clusters
  make_population <- function(k) {
    partition = list(
      group = NULL,
      result = NULL
    )
    
   # library(gelato)
    
    partition$group <- gelato::normalizeVector(sample(c(1:k), len, replace = TRUE))
    
    names(partition$group) = names
    
    partition$result <- eval.fun(mydata, partition$group)
    
    return(partition) 
  }
  
  
  #no of clusters
  #noc <- seq(len, 1, by = -1 * ceiling(len/size))  
  
  noc <- sample(2:(len-1), size)
  
  popl <- lapply(noc, make_population)
  
  find_best_partition <- function(partition) {
    
    eval_state <- function(state) {
      
      return(list(group= state, result = eval.fun(mydata, state)))
    }
  
    repeat {
      
      children <- climb.hill(partition$group, neighbors=neighbors)
      
      if (is.null(children)) 
        (break)()

      
      children_evaluated <- apply(children, 1, function(vec) {eval_state(vec)})
      
      print(length(children_evaluated))
      
      children_results = sapply(children_evaluated, function(x) x$result)
      
      children_groups = sapply(children_evaluated, function(x) x$group)
      
      
      local_best <- find.best(children_results)
      if (local_best$result > partition$result) {
        partition$result = local_best$result
        partition$group = children_groups[,local_best$idx]
      }
      else {
        (break)()
      }
    }
    
    return(partition)
    
  }
    
  sfInit(parallel=TRUE, cpus=ncpus)
  
  # 3. Wrapper, which can be parallelised. 
 
 sfExport("climb.hill", "eval.fun", "neighbors", "mydata", "find.best", "hybrid.linear.evaluator",
       "Cpp.hybrid.linear.MQ.LQ.evaluator", "Cpp.MQ.evaluator", "Cpp.LQ.evaluator")
  
  popl <- sfClusterApplyLB(popl, find_best_partition)

  print(length(popl))

  #popl <- lapply(popl, find_best_partition)
  
  #print(extensive_popl)
  
  # 7. Stop snowfall 
  sfStop()


  popl_results <- sapply(popl, function(x) x$result)
  popl_groups <- sapply(popl, function(x) x$group)

  best$result <- max(popl_results)
  
  best_result_index <- which.max(popl_results)
  
  #TODO I m worried about multiple best indices
  best$groups = popl_groups[,best_result_index]


  return(best)
}


#find unexplored zones/clusters
find.unexplored.zones <- function(noc, best_noc, size, len) {
  
  explored = c()
  
  all.zones = 1:len
  
  for (i in 1:size) {
   
    print(noc[i])
    
    print(best_noc[i])
    
    explored <- c(explored, noc[i]:best_noc[i])
  
    
  }
  return(all.zones[which(!(all.zones %in% explored))])
  
}

# returns indicies
find.subset <- function(subsets.matrix, subset) {
  
  subset = as.vector(subset)
  len = length(subset)
  
  if(len == 0)
    stop("Empty atrributes subset.")
  if(dim(subsets.matrix)[2] != len)
    stop("Incorrect dimensions.")
  
  if(dim(subsets.matrix)[1] == 0)
    return(as.integer(NULL))
  
  cond = rep(TRUE, dim(subsets.matrix)[1])
  for(i in 1:len)
    cond = cond & (subsets.matrix[,i] == subset[i])
  

  return(which(cond))
}


climb.hill <- function(parent, neighbors) {
  
  require(gelato)
  
  cols <- length(parent)
  # if(cols <= 0)
  #  stop("Parent groups set cannot be empty.")
  
  k <- max(parent)
 # if(k <= 0)
  #  stop("Parent groups not partitioned correctly.")
  
  #It's a single cluster, no neighbours can be derived
  if (k==1)
    return(NULL)

 # if (neighbors < 1)
  m <- gelato::computeRandomNeighbors(parent, cols, k, ceiling((cols * k) * neighbors))
 # else
 #   m <- compute.all.neighbors(parent, cols, k)
  
 return(m)
  
}

compute.random.neighbors <- function(parent, cols, k, n) {
  
  cols.indices <- sample(c(1:cols), n, replace= TRUE) 
  
  cols.values <- sample(c(1:k), n, replace = TRUE)
  
  m = matrix(ncol = cols, nrow = 0, byrow = TRUE, dimnames=list(c(),names(parent)))

  for(i in 1:n) {
    p <- parent
    p[cols.indices[i]] <- cols.values[i]
    p <- normalize.sample(p)
    m <- rbind(m, p)
  }
  
  return(m)
  
}

#TODO need to sort these by cols
compute.all.neighbors <- function(parent, cols, k) {
  #optimized to give faster result
  contains_matches <- function(x, e) {
    if (dim(x)[1] == 0)
      return(FALSE)
    
    for(i in dim(x)[1]:1)
      if (all(x[i,] == e))
        return(TRUE)
    return(FALSE)
  }
  
  m = matrix(ncol = cols, nrow = 0, byrow = TRUE)
  
  idx <- c(0)
  
  for (i in 1:cols)   
  {
    j<- c(1)
    while(j <= (k-1))
    {
      v <- parent
      
      #rewrite this
      if (((v[i] + j) %% (k+1)) != 0)
        v[i]<- (v[i] + j) %% (k+1)
      else
        v[i] <- 1
      
      s <- normalize.sample(v)
      
      if (!contains_matches(m,s))
        m <- rbind(m, s)
      
      j <- j+1
      
      idx <- idx+1
    }
  }
  
  return(m)
  
}

find.best <- function(results, subset = rep(TRUE, length(results))) {
  best = list(
    result = NULL,
    idx = NULL
  )
  
  w = which(subset)
  if(length(w) > 0) {
    children_results = results[w]
    best$result = max(children_results)
    best$idx = w[which.max(children_results)]
  }
  return(best)
}

#Precondition: it is a document term matrix
remove.documents <- function(bow, doc.names) {
  
  current.doc.names = rownames(bow)
  
  for (i in dim(bow)[1]:1)
    if (!(current.doc.names[i] %in% doc.names)) 
      bow <- bow[-i,]
  
  return(bow)
}


#precondition: dimenstions of both representations are both named
# rownames(mydatai) = colnames(mydatai) => square matric
make.compatible <- function(mydata) {
  mydata <- lapply(mydata, function(x) {x[order(rownames(x)), order(colnames(x))]})
  
  for (i in 1:length(mydata)) {
    iter <- 0
    
    for (j in 1:length(mydata)) {
      if (i==j)
        next
    
      for (k in dim(mydata[[i]])[1]:1) {
        name <- rownames(mydata[[i]])[k]

        if (!(name %in% rownames(mydata[[j]]))) {
          print(name)
          mydata[[i]] <- mydata[[i]][-k,]
          mydata[[i]] <- mydata[[i]][,-k]
          iter <- iter +1
          print(iter)
          
        }       
      }
    }
    
    cat("The number of rows and cols omitted in mydata1 ", iter, "\n")
    
  }
  
 return(mydata)
}


#precondition: a list of positive numbers 
#FIXME: need to fix the index to max_seen

normalize.sample <- function(vec, lb=1) {    
  if (length(vec) <= 0) 
    stop("Empty sample")    
  idx <- lb
  mv <- max(vec)
  for (i in 1: length(vec)) {
    x <- vec[i]
    
    if (x >= idx) {
      
      vec <- replace(vec, vec==idx, mv+1)
      
      mv <- mv+1
      vec <- replace(vec, vec==vec[i], idx)        
      idx <- idx+1       
    }
  }   
  return(vec)
}

import.xmi.call.graph <- function(filename) {
  
  require(XML)
  doc <- xmlRoot(xmlParse(filename))
  
  node.set = getNodeSet(doc, "//programunits")
  
  node.set.name = sapply(node.set, xmlGetAttr, "name")
  
  punits = unique(node.set.name)
  
  len = length(punits)
  
  m = matrix(c(0), dimnames= list(punits, punits), ncol=len, nrow=len, byrow = TRUE)

  
  #callers
  caller.set = getNodeSet(doc, "//programunits[@calls]")
  
  for(i in 1:length(caller.set)) {
    
    src.node = xmlGetAttr(caller.set[[i]], "name")
    
    tgt.ref <- xmlGetAttr(caller.set[[i]], "calls")
    
    #TODO need to make this the substring, starting from index of last '.'
    #Also check to ensure it is programuntis typed nodes
    
    tgt.idx <- substring(tgt.ref, 17)
    
    m[src.node, node.set.name[as.numeric(tgt.idx)+1]] <- 1
    
  }
  return(m)
}


find.singletons <- function(m){
  
  nrows = dim(m)[1]
  
  sumrows <- apply(m, 1, sum)
  sumcols <- apply(m, 2, sum)
  
  
  all <- sumrows+ sumcols
  
  for (i in 1:length(all)) 
    if (all[i] == 0)
    print(all[i])
  
  return(names(all[all == 0]))
  
}

export.matrix.bunch <- function(x) {

  require(reshape)

  
  m <- melt(x)
  output <- subset(m[m$value == 1,], select = c(X1,X2))
  
  print(output)
  write.table(output, "export.txt", sep=" ", quote = FALSE, eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE) 
}

import.bunch.matrix <- function(filename, exclude.ext=c()) {
  
  fulltext <- readLines(filename, encoding="UTF-8")

  if (length(exclude.ext) > 0) {

    len <-length(fulltext)
    
    for (i in len:1) {
      
      src.tgt <- unlist(strsplit(fulltext[i], " "))
      nodes <- cbind()
      if (exclude.scope.ext(src.tgt[2], exclude.ext))        
          fulltext <- fulltext[-i]

    }
  }
  
  
  nodes <- unique(unlist(lapply(fulltext, function(x) {unlist(strsplit(x, " "))})))
  
  m = matrix(0, nrow = length(nodes), ncol= length(nodes), dimnames = list(as.vector(nodes), as.vector(nodes)))
  
  print(length(nodes))
  
  for (i in 1:length(fulltext)) {

    src.tgt <- unlist(strsplit(fulltext[i], " "))
    m[src.tgt[1], src.tgt[2]] <-m[src.tgt[1], src.tgt[2]] +1      
  }

  return(as.matrix(m))
  
}

exclude.scope.ext <- function(text, extensions) {
  len <- length(extensions)  
  cond = rep(FALSE, len)
  for(i in 1:len)
    cond = cond || (grepl(extensions[i], text))
  
  return(cond)
  
}


import.call.graphs <- function(filenames, exclude.ext=c()){

  cfgs <- lapply(filenames, import.bunch.matrix , exclude.ext)
  
  
  names = c()
  
  for (i in 1:length(cfgs))
    
    names <- union(names, rownames(cfgs[[i]]))

  print(length(names))
  
  m = matrix(0, nrow = length(names), ncol= length(names), 
             dimnames = list(names, names))
  
  for(k in 1:length(cfgs)) 
    for (i in 1:dim(m)[1])
      if (names[i] %in% rownames(cfgs[[k]]))
      for (j in 1:dim(m)[2]) 
        if  (names[j] %in% rownames(cfgs[[k]]))
          m[i, j] <- m[i, j] + cfgs[[k]][names[i], names[j]]
      
    
  return(m)
  
  
}


compute.rank <- function(n , m) {
  #r = (n * m) ^ 0.2
  return(ceiling((n * m) ^ 0.2))
}

apply.hill.climbing.search <- function(prname, neighbors, size , ncpus, tries=1) {
  
  
  setwd("~/workspace")
  
  require(igraph)
  require(Rcpp)
  require(lsa)
  
  #Load the priori decomposition
  decomposition <- read.csv(paste("benchmark", prname ,"decomposition.csv", sep="/"), sep=",",  header = TRUE)
  priori.decomp <- decomposition$x
  names(priori.decomp) <- decomposition$X
  priori.decomp <- normalizeVector(priori.decomp)
  
  
  #Load the adjacency matrix
  extensions= c("java/", "org/xml/", "javax/")
  cfg <- import.bunch.matrix(paste("benchmark", prname ,"dep_graph.txt", sep="/"), exclude.ext=extensions)
  #cfg <- read.table("benchmark/jedit-5.1.0/cfg.csv", sep=",", row.names = 1, header = TRUE, check.names = FALSE)
  #cfg <- unweight.adjacency(cfg)
  cfg <- cfg[intersect(rownames(cfg),names(priori.decomp)), intersect(colnames(cfg),names(priori.decomp))]
  
  #Find intersection with priori decomp
  cfg <- cfg[intersect(rownames(cfg),names(priori.decomp)), intersect(colnames(cfg),names(priori.decomp))]
  
#   return(cfg)
  #Load and compute the lexical sim matrix
  bow <- read.table(paste("benchmark", prname , "mydata-idf-BoW-matrix.csv", sep="/"), sep=",", row.names = 1, header = TRUE, check.names = FALSE)
  bow <- remove.documents(bow, rownames(cfg))
  write.table(bow, file = paste("benchmark", prname , "mydata-compatible-idf-BoW-matrix.csv", sep="/"),row.names=TRUE, col.names=NA,sep=",", quote=FALSE)
  
  #Apply LSA
  ndims = compute.rank(dim(bow)[1], dim(bow)[2])
  print("No of Dimensions:")
  print(ndims)
  
  space1 = lsa(bow, ndims)
  lsa.bow = as.textmatrix(space1)
  lexsim = cosine(t(lsa.bow))
  
  #write.table(lexsim, file ="mydata-LSA-Cos-similarity-matrix.csv",row.names=TRUE, col.names=NA,sep=",", quote=FALSE)
  #lexsim <- read.table(paste("benchmark", prname , "mydata-Cos-similarity-matrix.csv", sep="/"), sep=",", row.names = 1, header = TRUE, check.names = FALSE)
  
  
  
  #Build mydata
  mydata = list(cfg=as.matrix(cfg), lexsim=as.matrix(lexsim))
  mydata <- make.compatible(mydata)
  mydata$cfg <- remove.diagonal(mydata$cfg)
  mydata$lexsim <- remove.diagonal(mydata$lexsim)  
  
  N <- dim(mydata$cfg)[1]
  
  #Perform local search from weights 0 to 1
  ws <- c(seq(0, 0.5, 0.2), 0.5, seq(0.6, 1, 0.2))
  
  clusters <- apply.weighted.hybrid.evaluator(ws, mydata, neighbors=neighbors, size = size, ncpus=ncpus, tries)
  
  #Compute distance  
  priori.decomp <- find.intersection(priori.decomp, clusters[,1]$groups)
  priori.decomp <- normalizeVector(priori.decomp)
  
  mojo <- apply(clusters, 2, function(c) compute.MoJo(c$groups, priori.decomp))
  
  mojosim <- sapply(mojo, function(m) 1 - (m/N))
  
  nmis <- apply(clusters, 2, function(c) compute.NMI(c$groups, priori.decomp))
  
  purities <- apply(clusters, 2, function(c) compute.cluster.purity(c$groups, priori.decomp))
  
  
  #Save Results
  results = matrix(nrow = 0, ncol = length(ws))
  score <- clusters["result",]
  results <- rbind(results, score)
  results <- rbind(results, mojo)
  results <- rbind(results, mojosim)
  results <- rbind(results, nmis)
  results <- rbind(results, purities)
  colnames(results) <- ws
  
  write.table(results, file =paste("benchmark", prname ,"singleobj-results.csv", sep="/"), row.names=TRUE, col.names=NA,sep=",", quote=FALSE)
  
  #Store groups
  groups = matrix(nrow = 0, ncol = dim(mydata$cfg)[1])
  for (i in 1:length(ws))
    groups <- rbind(groups, clusters[,i]$groups)
  
  rownames(groups) <- ws  
  colnames(groups) <- rownames(mydata$cfg)
  print(groups)
  write.table(groups, file =paste("benchmark", prname ,"singleobj-groups.csv", sep="/"),row.names=TRUE, col.names=NA,sep=",", quote=FALSE)
  
  #Save Objective Scores
  scores = matrix(nrow = length(ws), ncol = 2, dimnames= list(ws, c("MQ", "LQ")))
  for (i in 1:length(ws)) {
    group <- clusters[,i]$groups
    
    scores[i,1] <- Cpp.MQ.evaluator(mydata, group)
    scores[i,2] <- Cpp.LQ.evaluator(mydata, group)
    
  }
  print(scores)
  write.table(scores, file =paste("benchmark", prname ,"singleobj-scores.csv", sep="/"),row.names=TRUE, col.names=NA,sep=",", quote=FALSE)
  
}

analyze.singleobj.results <- function(prname) {
  
  require(scatterplot3d)
  
  setwd("~/workspace")
  
  #Load the results
  results <- read.csv(paste("benchmark", prname ,"singleobj-results.csv", sep="/"), sep=",",  header = TRUE)
  
  rownames(results) <- results$X
  
  results <- results[,-1]
  
  weights = c(0, 0.2, 0.4, 0.5, 0.6, 0.8, 1)
  
  colnames(results) <-  weights 
  
  
  scores <- read.csv(paste("benchmark", prname ,"singleobj-scores.csv", sep="/"), sep=",",  header = TRUE)
  
  
  if(dim(results)[2] != dim(scores)[1])
    stop("incompatible dataset")
  

  results <- rbind(results, scores$MQ)
  results <- rbind(results, scores$LQ)
  
  rownames(results[4,]) = "MQ"
  rownames(results[5,]) = "LQ"
  
  groups <- read.csv(paste("benchmark", prname ,"singleobj-groups.csv", sep="/"), sep=",",  header = TRUE)
  
  rownames(groups) <- groups$X
  groups <- groups[,-1]
  
  
  average_case <- function(mydata, groups) {
    case = list(noc= 0, mojosim = 0)
    
    case$mojosim = mean(unlist(mydata[3,]))
    case$noc = mean(apply(groups, 1, function(x) max(normalizeVector(x))))
    
    return(case)
    
  }
  
  find_case <- function(mydata, groups, fun) {
    case = list(lambda = 0, noc= 0, mojosim = 0)
    
    case$mojosim = fun(unlist(mydata[3,]))
    
    for (i in 1:dim(mydata)[2])
      if(mydata[3,i] == case$mojosim) {
        
        case$lambda <- colnames(mydata)[i]

        case$noc <- max(normalizeVector(unlist(groups[i,])))
        break()
      }
    
    return(case)
    
  }
  
  best_case <- find_case(results, groups, max)
  worst_case <- find_case(results, groups, min)
  average_case <- average_case(results, groups)
 
  LQ_output <- paste(max(normalizeVector(unlist(groups[1,]))), 
                     round(results[3,1], digits = 2), sep = "&")
  MQ_output <- paste(max(normalizeVector(unlist(groups[7,]))), 
                     round(results[3,7], digits = 2), sep = "&")
  output_50 <- paste(max(normalizeVector(unlist(groups[4,]))), 
                     round(results[3,4], digits = 2), sep = "&")
  
  best_output <- paste(best_case$lambda, best_case$noc, round(best_case$mojosim, digits = 2), sep = "&")
  worst_output <- paste(worst_case$lambda, worst_case$noc, round(worst_case$mojosim, digits = 2), sep = "&")
  average_output <- paste(round(average_case$noc, digits = 2), round(average_case$mojosim, digits = 2), sep = "&")
  
  
  output <- paste(LQ_output, MQ_output, best_output, worst_output, average_output, "&", sep="&")
  
  print("printing All results")
  print(output)
  
  print("printing 50% results")
  print(output_50)
  
  
  # with(results, {
  #mojosim <- unlist(results[3,])
  # plot(weights, mojosim , type="l", lwd=3)       # x and y axis
  
  #polynomial of degree two
  #fit3 <- lm(mojosim ~ poly(weights, 10, raw=TRUE))
  #print(summary(fit3))
  #pol3 <- function(x) fit3$coefficient[4]*x^3 + fit3$coefficient[3]*x^2 + fit3$coefficient[2]*x + fit3$coefficient[1]
  #curve(pol3, col="red", lwd=2)
  #points(weights, predict(fit), type="p", col="red", lwd=2)
  #points(weights, mojosim, type="p", lwd=3)
  #points(weights, predict(fit3), type="l", lwd=2)
  #fit <- lm(MoJoSim ~ MQ+LQ)
  #s3d$plane3d(fit)
  
  #  text(s3d.coords$x, s3d.coords$y,             # x and y coordinates
  # labels=row.names(mtcars),               # text to plot
  #       cex=.5, pos=4)           # shrink text 50% and place to right of points)
  # })
  
  #dev.off()
  
  return(results)  
  
}

plot.singleobj.multiobj <- function(r, prname) {
  t1<-theme(                              
    plot.background = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    panel.background = element_blank(),
    axis.line = element_line(size=.4)
  )
  
  setwd("~/workspace/sci.hage0101.GELATO")
  require("splines")
  
  
}


predict.singleobj <- function(r, prname) {
  setwd("~/workspace")
  require(splines)
  require(ggplot2)
  
  t1<-theme(                              
    plot.background = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    panel.background = element_blank(),
    axis.line = element_line(size=.4),
    axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold"),
    legend.title=element_text(size=10, vjust=-12)
    #     guide_colourbar.title = element_text(draw.ulim = FALSE, draw.llim = FALSE)
    #    legend.key.height=unit(3,"line"),
    #    legend.key.width=unit(3,"line")
  )
  
  #print(y)
  #Single-Objective
    
  Weight = c(0, 0.2, 0.4, 0.5, 0.6, 0.8, 1)
  MoJoSim = unlist(r[3,])
  so.results = data.frame(Weight, MoJoSim)
  
  
  #Multi-Objective
  results <- read.csv(paste("benchmark", prname ,"multiobj-results.csv", sep="/"), sep=",",  header = TRUE)
  
  results <- results[,-1]
  
  scores <- read.csv(paste("benchmark", prname ,"multiobj-scores.csv", sep="/"), sep=",",  header = TRUE)
  
  scores <- scores[,-1]
  
  if(dim(results)[1] != dim(scores)[1])
    stop("incompatible dataset")
  
  results <- cbind(results, scores)
  results <- unique( results[ , 1:5 ] )
  
  
  #Project to One-dim
  d <- dist(results[,4:5]) # euclidean distances between the COLS
  fit <- cmdscale(d,eig=TRUE, k=1) # k is the number of dim
  x <- fit$points[,1]
  
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  
  Weight <- range01(x)
  MoJoSim = results[,3]
  
  mo.results = data.frame(Weight, MoJoSim)
  
  plot <- ggplot(NULL, aes(Weight, MoJoSim)) + t1 #+ scale_y_continuous(limits=c(0.42, 0.51))
  plot <- plot + geom_point(data=so.results) +
    stat_smooth(data = so.results, method = "lm", se= FALSE, formula = y ~ ns(x,6), colour = "blue")
  
#   plot <- plot +
#     stat_smooth(data = mo.results, method = "lm", se= FALSE, formula = y ~ ns(x,15), colour = "red")
  
  
  
  ggsave(filename= paste("benchmark", prname ,"single_scatterplot.png", sep="/"), plot=plot, pointsize = 15, width = 10, height = 10)
  
  dev.off()
}

run.experiments <- function() {
  
 # apply.hill.climbing.search("apache-ant-1.9.3", 0.8, 40, 8, 5)
 ## apply.hill.climbing.search("jhotdraw-7.0.6", 0.8, 40, 8, 5)
 # apply.hill.climbing.search("junit-4.12", 0.8, 40, 8, 5)
  #apply.hill.climbing.search("jfreechart-1.2.0", 0.8, 40, 8, 3)
  apply.hill.climbing.search("jdom-2.0.5", 0.8, 40, 8, 5)
  apply.hill.climbing.search("apache-log4j-1.2.17", 0.8, 32, 8, 5)
  #apply.hill.climbing.search("jedit-5.1.0", 0.4, 22, 8, 1)
  #apply.hill.climbing.search("hadoop-0.20.2", 0.4, 20, 8, 2)
  
  #apply.hill.climbing.search("", 0.80, 48, 8, 10)
}
