
#TESTED
load_BoF <- function(prname, id2t = c(T, F) ){
  
  
  stopifnot(id2t[1] || id2t[2])
  
#   if (!id2t[1] && !id2t[2])
#     stop("one of the two identifiers or types must be set to TRUE")
  
  setwd("~/workspace")
#   myBoF = read.csv(paste("benchmark", prname , "BoF", paste(prname, "BoF.csv", sep="-"), sep="/"),  sep = ",")
#   rownames(myBoF) <- myBoF[,1]
#   myBoF <- myBoF[,-1]
#   myBoF <- data.matrix(myBoF)
  
  myBoF = read.csv(paste("benchmark", prname , "BoF", paste(prname, "BoF.csv", sep="-"), sep="/"),  sep = ",", header=FALSE, stringsAsFactors = FALSE)
  

  identifiers <- unlist(myBoF[1,2:dim(myBoF)[2]])
  types <- unlist(myBoF[2,2:dim(myBoF)[2]])

  myBoF <- myBoF[-c(1,2),]


  rownames(myBoF) <- myBoF[,1]
  myBoF <- myBoF[,-1]
  myBoF <- data.matrix(myBoF)
  
  #combine the idetifiername and type name
  if (id2t[1] && id2t[2]){
    
    combined_names <- paste(identifiers, types, sep = ":")
    
    colnames(myBoF) <- combined_names
    
    return(myBoF)
  }

  #Only return the identifiers
  if (id2t[1]){
    unique_identifiers <- unique(identifiers)
    
    m <- matrix(0, nrow=dim(myBoF)[1], ncol = length(unique_identifiers))
    
    
    for (i in 1:length(unique_identifiers)){
      
      
      indices <- which(identifiers == unique_identifiers[i])
      
      
      s <- myBoF[,indices]
      
      if (is.vector(s))
        m[,i] <- s
      else
        m[,i] <- rowSums(s)
      
      
    }
    
    colnames(m) <- unique_identifiers
    rownames(m) <- rownames(myBoF)
    
    return(m)
    
  }
  #Only return the types
  if (id2t[2]){
    unique_types <- unique(types)
    
    m <- matrix(0, nrow=dim(myBoF)[1], ncol = length(unique_types))
    
    
    for (i in 1:length(unique_types)){
      
      
      indices <- which(types == unique_types[i])
      
      
      s <- myBoF[,indices]
      
      if (is.vector(s))
        m[,i] <- s
      else
        m[,i] <- rowSums(s)
      
      
    }
    
    colnames(m) <- unique_types
    rownames(m) <- rownames(myBoF)
    
    return(m)
    
  }

return(myBoF)
}

compute_semantic_similarity <- function(Adj, eval.fun, weights){
  #remove Adj of types to identifier names
  #find the first occurence of a type name in the colnames(ADj), set everything else to 0
  names <- colnames(Adj)
  #   names <- rownames(Adj)
  
  #FIND THE STARTING INDEX OF TYPE NAMES  
  classTypeIndex <- which(unlist(gregexpr(pattern ="\\.",names)) > 0)[1]
  primitiveTypeIndices <- which(names %in% c("float", "int", "char", "byte", "void", "double", "boolean"))
  startIndex <- min(c(classTypeIndex, primitiveTypeIndices))
  
  
  print("starting Index")
  print(startIndex)
  
  #size of dictionary
  D <- startIndex -1
  #size of types
  V <- dim(Adj)[2] - startIndex +1
  
  #Also remove contaninmet dependencies from method names to local variables
  #   for (i in 1:dim(Adj)[1])
  #     for (j in 1:(startIndex-1))
  #       Adj[i,j] <- 0
  
  S <- Adj[startIndex:dim(Adj)[1], startIndex:dim(Adj)[2]] 
  # dimnames(S) <- dimnames(Adj[startIndex:dim(Adj)[1], startIndex:dim(Adj)[2]])

  print("printing weights")
  print(weights)

  # Get the type contents
  C <- Adj[startIndex:dim(Adj)[1], 1:(startIndex-1)] 
  # dimnames(C) <- dimnames(Adj[startIndex:dim(Adj)[1], 1:(startIndex-1)])
  
  #DONE make this a higher function argument
  type_sim <- eval.fun(S, C)
  
  if (sum(weights) > 1){
    r <- compute_combined_IPO_ISA_SN(Adj, startIndex, type_sim)
    
    return(list(sim=r$sim, normalized_sim=r$normalized_sim))
  }
  
  if (weights[1] > 0){
    ISA <- compute_ISA_SN(Adj, startIndex, type_sim)

  }
  else
    ISA <- list(sim=0, normalized_sim=0)
  
  
  
  if (weights[2] > 0)
    IPO <- compute_IPO_SN(Adj, startIndex, type_sim)
  else
    IPO <- list(sim=0, normalized_sim=0)
  
  #then combine the two Synsets and get one similarity matrix
 
  
  list(sim=(weights[1] * ISA$sim + weights[2] * IPO$sim), normalized_sim = (weights[1] * ISA$normalized_sim + weights[2] * IPO$normalized_sim) )
}

#combine combined ISA-IPO semantic similarity between terms for SN
#eval.fun is the evaluation function
compute_combined_IPO_ISA_SN <- function(Adj, startIndex, type_sim)
{
  names <- rownames(Adj)

  outIndices <- list()
  outISA_IPO <- list()
  outdegree <- c()
  
  for (i in 1:(startIndex-1)){
    
    outIndices[[i]] <- c(which(Adj[startIndex:dim(Adj)[1], i]>0), which(Adj[i, startIndex:dim(Adj)[2]]>0))
    
    outISA_IPO[[i]] <- c(Adj[names(outIndices[[i]]),i], Adj[i,names(outIndices[[i]])])
    
    outdegree[i] <- sum(outISA_IPO[[i]])        
  }
  
  result <- compute_term_similarity(outIndices, outISA_IPO, outdegree, type_sim)
  
  dimnames(result$sim) <- dimnames(Adj[1:(startIndex-1), 1:(startIndex-1)])
  dimnames(result$normalized_sim) <- dimnames(Adj[1:(startIndex-1), 1:(startIndex-1)])
  
  result
}

#compute ISA semantic similarity between terms for SN
#eval.fun is the evaluation function
compute_IPO_SN <- function(Adj, startIndex, type_sim){
  
  names <- rownames(Adj)
  #   names <- rownames(Adj)
  
  outIndices <- list()
  outIPO <- list()
  outdegree <- c()
  
  for (i in 1:(startIndex-1)){
    
    outIndices[[i]] <- which(Adj[startIndex:dim(Adj)[1], i]>0)
    
    outIPO[[i]] <- Adj[names(outIndices[[i]]),i]
    
    outdegree[i] <- sum(outIPO[[i]])        
    
  }
  
  result <- compute_term_similarity(outIndices, outIPO, outdegree, type_sim)
  
  dimnames(result$sim) <- dimnames(Adj[1:(startIndex-1), 1:(startIndex-1)])
  dimnames(result$normalized_sim) <- dimnames(Adj[1:(startIndex-1), 1:(startIndex-1)])

  result
}

# The initial implementation set the similarity score between any two identifiers to 1, 
# as long as they had a single edge pointing to the same type
# In this version, this is no longer the case, rather the max of all possible combinations is computed.  
compute_term_similarity  <- function(outIndices, outlinks, outdegree, type_sim){
  
  #length of all parameters must mach
  len <- length(outIndices)
  
  #Make sure the diagonal of type_sim is 1
  diag(type_sim) <- 1
  stopifnot(all(type_sim) <= 1)
  
  sim <- matrix(0, nrow=len, ncol=len)
  normalized_sim <- matrix(0, nrow=len, ncol=len)
  
  for (i in 1:len){
    outIndices1 <- outIndices[[i]]
    if (length(outIndices1) < 1)
      next
    
    outlinks1 <- outlinks[[i]]
    outdegree1<- outdegree[i]
    
    for (j in i:len){
      if (i == j) #skip if i==j
        next
      
      outIndices2 <- outIndices[[j]]
      if (length(outIndices2) < 1)
        next
      
      outlinks2 <- outlinks[[j]]
      outdegree2<- outdegree[j]
      
      # READ the implementation change in the function header comment
      # if (length(intersect(outIndices1,outIndices2))>0){
      #   sim[i,j] <- 1
      #   normalized_sim[i,j] <- 1
      #   next
      # }
      
      max_sim <- 0 
      max_normalized_sim <- 0
      
      for(k in 1:length(outIndices1)){
        for(l in 1:length(outIndices2)){
          
          sim_val <- type_sim[outIndices1[k], outIndices2[l]]
          
          normalized_sim_val <- 0.5 * sim_val * ((outlinks1[k]/outdegree1) + (outlinks2[l]/outdegree2))
          
          if (sim_val > max_sim)
            max_sim <- sim_val
          
          if (normalized_sim_val > max_normalized_sim)
            max_normalized_sim <- normalized_sim_val
          
        }
      }
      sim[i,j] <- max_sim
      normalized_sim[i,j] <- max_normalized_sim
    }
    
  }
  
  
  sim <- fill_lower_diagonal(sim)
  normalized_sim <- fill_lower_diagonal(normalized_sim)
  diag(normalized_sim) <- 1
  
  return(list(sim=sim, normalized_sim=normalized_sim))
}

#compute ISA conceptual density between terms for SN
#eval.fun is the evaluation function
compute_ISA_SN <- function(Adj, startIndex, type_sim){
  
  names <- rownames(Adj)

  outIndices <- list()
  outISA <- list()
  outdegree <- c()
  
  D <- startIndex-1
  
  for (i in 1:(startIndex-1)){
    
    outIndices[[i]] <- which(Adj[i, startIndex:dim(Adj)[2]]>0)
    
    outISA[[i]] <- Adj[i,names(outIndices[[i]])]
    
    outdegree[i] <- sum(outISA[[i]])        
      
  }
  
  result <- compute_term_similarity(outIndices, outISA, outdegree, type_sim)
  
  dimnames(result$sim) <- dimnames(Adj[1:(startIndex-1), 1:(startIndex-1)])
  dimnames(result$normalized_sim) <- dimnames(Adj[1:(startIndex-1), 1:(startIndex-1)])
  
  result
}

compute_normalized_Resnik_similarity <- function(S,C) {
  compute_Resnik_similarity(S, C, normalize_corpus=T)
}

#Input: 
# S is type hypernymy graph
# C is type content
compute_Resnik_similarity <- function(S, C, normalize_corpus=F){
  
  stopifnot(rownames(S) == rownames(C))
  
  if(normalize_corpus)  
    C <- compute_normalize_corpus(C)
  
  sum_BoT <- rowSums(C)
  
  information_content <- -log(sum_BoT /sum(sum_BoT))
  
  information_content <- unlist(lapply(information_content, function(ic) if (is.infinite(ic)) 0 else ic))
  
  sim_matrix<- matrix(0, nrow=dim(S)[1], ncol=dim(S)[2])
  
  for (i in 1:dim(S)[1])
    for(j in i:dim(S)[2])
      if (i != j)
      {    
        #check if one is subtype of other
        LCHs <- compute_least_common_hypernym(S, i, j)
        
        if (length(LCHs) < 1)
          next
        
        maxSim <- 0
        
        for (k in 1 : length(LCHs)){
          
          LCH <- LCHs[k]
          
          sim <- information_content[LCH]
          
          if (sim > maxSim)
            maxSim <- sim       
        } 
        
        sim_matrix[i,j] <- maxSim       
      }
  
  sim_matrix <- fill_lower_diagonal(sim_matrix)  
  dimnames(sim_matrix) <- dimnames(S)
  
  sim_matrix  
}


compute_normalize_corpus <- function(C, nat.langs="english", prog.langs="java"){
  docs =c()
  identifierNames <- colnames(C)
  
  for (i in 1:dim(C)[1]){
    d <- ""
    for (j in 1:dim(C)[2]){
      occur <- C[i,j]
      identifierName <- identifierNames[j]
      if (occur > 0){
        for (k in 1:occur)
          d <- paste(d, identifierName, sep=" ")
        
      }
    }
    docs[i] <- d
  }
  
  mydata <- lapply(docs, function(d) strip.java.text(d))
  
  mydata <- lapply(mydata, function(d) if (is.na(d)) "" else d )
  
  names(mydata) <- rownames(C)
  
  for (i in 1:length(prog.langs))
    mydata <- prepare.prog.lang.list(mydata, prog.langs[i])
  
  for (i in 1:length(nat.langs)) 
    mydata <- prepare.natural.lang.list(mydata, nat.langs[i])
  
  mydata.BoW.list <- make.BoW.list(mydata)
  
  mydata.BoW.frame <- make.BoW.frame(mydata.BoW.list, names(mydata.BoW.list))
  
  mydata.BoW.frame
}

compute_normalized_Lin_similarity <- function(S,C) {
  compute_Lin_similarity(S, C, normalize_corpus=T)
}

#Input: 
# S is type hypernymy graph
# C is type content
compute_Lin_similarity <- function(S, C, normalize_corpus=F){
  stopifnot(rownames(S) == rownames(C))
  
  if(normalize_corpus)   
    C <- compute_normalize_corpus(C)
    
   
  sum_BoT <- rowSums(C)
    
  information_content <- -log(sum_BoT /sum(sum_BoT))
  
  information_content <- unlist(lapply(information_content, function(ic) if (is.infinite(ic)) 0 else ic))


  sim_matrix<- matrix(0, nrow=dim(S)[1], ncol=dim(S)[2])
  
  
  
  for (i in 1:dim(S)[1])
    for(j in i:dim(S)[2])
      if (i != j)
      {
        if (information_content[i] + information_content[j] <= 0)
          next
        
        
        #check if one is subtype of other
        LCHs <- compute_least_common_hypernym(S, i, j)
        
        if (length(LCHs) < 1)
          next
        
        maxSim <- 0
        
        for (k in 1 : length(LCHs)){
          
          LCH <- LCHs[k]
          
          sim <- 2 * information_content[LCH] / (information_content[i] + information_content[j])
          
          if (sim > maxSim)
            maxSim <- sim
          
        }        
        sim_matrix[i,j] <- maxSim
        
      }

  sim_matrix <- fill_lower_diagonal(sim_matrix)  
  dimnames(sim_matrix) <- dimnames(S)

  sim_matrix    
}


compute_Wu_Palmer_similarity <- function(S, C, rootNode ="java.lang.Object"){
  require(compiler)
  require(igraph)
  
  g <- graph.adjacency(S, weighted=TRUE, mode="directed", diag=F)
  
  #this could be average, max and min
  SP <- shortest.paths(g, mode="all")
  
  SP[is.infinite(SP)] <- 0 
  
  sim_matrix <- matrix(0, nrow=dim(S)[1], ncol=dim(S)[2])
  
  root <- which(rownames(S)==rootNode)
  
  V <- dim(S)[1]
  roots <- unlist(lapply(1:V, function(v) compute_depth(S, root, v))) 
  
  CLCH <- cmpfun(compute_least_common_hypernym)

  for (i in 1:dim(S)[1])
    for(j in i:dim(S)[2])
      if (i != j)
      {
        
        #check if one is subtype of other
        LCHs <- CLCH(S, i, j)
        
        if (length(LCHs) < 1)
          next
        
        maxSim <- 0
        
        for (k in 1 : length(LCHs)){
          
          LCH <- LCHs[k]
          depth <- roots[LCH]

          distance1 <- SP[i, LCH]
          distance2 <- SP[j, LCH]

          sim <- 2*depth/(distance1 + distance2 + 2*depth)
          
          if (sim > maxSim)
            maxSim <- sim
        }
        
        sim_matrix[i,j] <- maxSim
        
      }
  
  
  sim_matrix <- fill_lower_diagonal(sim_matrix)
  
  dimnames(sim_matrix) <- dimnames(S)
  sim_matrix  
}


compute_Leacock_Chodorow_similarity <- function(S, C, rootNode ="java.lang.Object"){
  
  require(igraph)
  
  g <- graph.adjacency(S, weighted=TRUE, mode="directed", diag=F)
  
  #this could be average, max and min
  SP <- shortest.paths(g, mode="all")
  
  SP[is.infinite(SP)] <- 0 
  
  root <- which(rownames(S)==rootNode)
  
  sim_matrix <- matrix(0, nrow=dim(S)[1], ncol=dim(S)[2])
  
  V <- dim(S)[1]
  
  roots <- unlist(lapply(1:V, function(v) compute_depth(S, root, v)))
  
  for (i in 1:dim(S)[1])
    for(j in i:dim(S)[2])
      if (i != j)
      {
          
        max_depth <- max(roots[i], roots[j])
        distance <- SP[i, j]
          
        sim <- -log(distance/(2*max_depth))

        sim_matrix[i,j] <- sim     
      }
  
  sim_matrix <- fill_lower_diagonal(sim_matrix)
  
  dimnames(sim_matrix) <- dimnames(S)
  sim_matrix  
}

compute_inverted_path_length <- function(S, C, alpha=1){
 
  require(igraph)
  
  g <- graph.adjacency(S, weighted=TRUE, mode="directed", diag=F)
  
  #this could be average, max and min
  SP <- shortest.paths(g, mode="all")
  
  SP[is.infinite(SP)] <- 0 
  
  1 / (1 + SP^alpha)
}

#Input: Adj is a directed graph

#TODO what to do with subtypes
compute_conceptual_density <- function(S, C){
  library(compiler)
  
  branching_factor <- compute_branching_factor(S)
  
  CD <- matrix(0, nrow=dim(S)[1], ncol=dim(S)[2])
  
  CLCH <- cmpfun(compute_least_common_hypernym)
  CCD <- cmpfun(calculate_conceptual_density)

  for (i in 1:dim(S)[1])
    for(j in i:dim(S)[2])
      if (i != j)
        {
      
        #check if one is subtype of other
        LCHs <- CLCH(S, i, j)
        
        #If no common hypernym
        if (length(LCHs) < 1)
          next
        
        maxCD <- 0
        
        for (k in 1 : length(LCHs)){
          
          LCH <- LCHs[k]
          
          cd <- CCD(LCH, branching_factor)
          
          if (cd > maxCD)
            maxCD <- cd
        }   
      
        CD[i,j] <- maxCD
        
      }
  
  CD <- fill_lower_diagonal(CD)

  CD
}


fill_lower_diagonal <- function(sim_matrix){
  sim_matrix[lower.tri(sim_matrix)] <- t(sim_matrix)[lower.tri(sim_matrix)]
  
  sim_matrix  
}


#branching factor is a list ()

calculate_conceptual_density <- function(LCH, branching_factor){
    
  total_branching_factor <- branching_factor[[LCH]]$total_branching
#   print(total_branching_factor)
  total_children <- branching_factor[[LCH]]$total_children
#   print(total_children)
  avg_branching_factor <- total_branching_factor / total_children
  
  #Output: h for calculating the conceptual density
  estimate_depth <- function(avg_branching_factor){
    if (avg_branching_factor== 1)
      return(2)
    
    floor(log(2, base = avg_branching_factor))
    
  }
  
  
  h <- estimate_depth(avg_branching_factor)
  
  sum_squared_avg_branching_factor <- 0
  
  for (i in 0:h){
    
    sum_squared_avg_branching_factor <- sum_squared_avg_branching_factor + avg_branching_factor^i
  }
  
  sum_squared_avg_branching_factor/total_branching_factor
  
}


#Input: root, x, y should be index in adj
#Output: a list of nodes, one for each node, excluding 
compute_subhierarchy <- function(adj, root, x, y){
  require(igraph)

  g <- graph.adjacency(adj, mode='directed', diag=FALSE)  
  
  p <- get.shortest.paths(g, root, to = c(x,y), mode = c("in"))
  
  unique(unlist(p$vpath))
  
}

#FIXME: there could be another way of implementing this by finding all the reachable paths from a node, then finding the one
#minimal distance that reaches one possible root Node
#In case where a unique top node exists, the following algorithm: the shortest path from root node to the ther node suffices
#Input: Adjacency, 
#Output: 
compute_depth <- function(adj, root, x){
  
  require(igraph)
  
  g <- graph.adjacency(adj, mode='directed', diag=FALSE)
  
  shortest.paths(g, v=root, to=x, mode="in")
}


#Input: an adjacency graph Adj, LCA(x,y)
#Output single least common ancestors of x and y
compute_least_common_hypernym <- function(adj, x, y){
  
  #Let say you want to compute LCA(x,y) with x and y two nodes. Each node must have a value color and count, resp. initialized to white and 0.
  
  #Color all ancestors of x as blue (can be done using BFS)
  #Color all blue ancestors of y as red (BFS again)
  #For each red node in the graph, increment its parents' count by one
  #Each red node having a count value set to 0 is a solution.
  
  require(igraph)
  
  g <- graph.adjacency(adj, mode='directed', diag=FALSE,weighted=T)  
  
  V <- dim(adj)[1]
  
  blue <- subcomponent(g, x, "out")
  
  temp <- subcomponent(g, y, "out")
  
  
  red <- blue[blue %in% temp]
  
#   indices <- which(red)
  count <- rep(0, V)
  
  
  for (i in 1:length(red))
    count[which(adj[red[i],]>0)] <- count[which(adj[red[i],]>0)] + 1
  
#FIXME change this
if (length(red) < 1)
  return (c())

for (i in 1:length(red))
  if (count[red[i]] == 0 )
    return(which(rownames(adj) == red[i]$name))
}


#Input: adjacency matrix of a directed acyclic graph
#Output: computes branching factor for each node (type)
compute_branching_factor <- function(adj){
  
  require(igraph)
  
  g <- graph.adjacency(adj, mode='directed', diag=FALSE,weighted=T)  
  
  V <- dim(adj)[1]
  
  branching <- colSums(adj)
  
  r <- c()
  
  for (i in 1:V){
    reachable_nodes <- subcomponent(g, i, "in")
    
    
    total_children <- length(reachable_nodes)
    
    total_branching <- sum(branching[reachable_nodes])
    
    r[[i]] <- list(total_branching=total_branching, total_children=total_children)
    
  }
  
  r
    
}

#prnames <- c("")
run_in_batch_mode <- function(prnames){
  
  for(i in 1:length(prnames)){
    run_in_batch_mode_each_project(prnames[[i]])
    
  }
  
}

run_in_batch_mode_each_project <- function(prname, size=0.25){
  
  run_each_setting <- function(weights, eval.fun=NULL, normalized=T, lex.fun =NULL,  Adj, myBoF){
    require(igraph)
    require(GeLaToLab)
    
    sim_kernel <- NULL
    
    if (!is.null(eval.fun)){
    
      r <- compute_semantic_similarity(Adj, eval.fun, weights)
      
      if (normalized)
        sim_kernel <- r$normalized_sim
      else
        sim_kernel <- r$sim
    }
    
    #Calculate the identifier names set
    names <- colnames(Adj)
      
    #FIND THE STARTING INDEX OF TYPE NAMES  
    classTypeIndex <- which(unlist(gregexpr(pattern ="\\.",names)) > 0)[1]
    primitiveTypeIndices <- which(names %in% c("float", "int", "char", "byte", "void", "double", "boolean"))
    startIndex <- min(c(classTypeIndex, primitiveTypeIndices)) 
      
    print("starting Index")
    print(startIndex)
      
    #size of dictionary
    D <- startIndex -1
      
    SN_identifier_names <- tolower(rownames(Adj)[1:D])
    BoF_identifier_names <- tolower(colnames(myBoF))

    #Find common identifiernames between BoF and the Semantic Network
    identifierNames <- intersect(BoF_identifier_names, SN_identifier_names)
    
    #Eliminate identifier names in BoF that are not in SN
    identifierIndices <- which(tolower(colnames(myBoF)) %in% identifierNames)    
    myBoF <- myBoF[,identifierIndices]
    
    #remove classes with no identifiers, when combined with the semantic network
    
    myBoF <- myBoF[which(!apply(myBoF,1,FUN = function(x){all(x == 0)})),] 
    myBoF <- myBoF[order(rownames(myBoF)), order(colnames(myBoF))]

    if (!is.null(sim_kernel)){
      sim_kernel <- sim_kernel[order(rownames(sim_kernel)), order(colnames(sim_kernel))]
      
      #Eliminate identifier names in sim_kernel that are not in BoF
      identifierIndices <- which(tolower(colnames(sim_kernel)) %in% identifierNames) 
      
      print("some info:")
      print(dim(sim_kernel))
      print(length(identifierNames))      
      sim_kernel <- sim_kernel[identifierIndices,identifierIndices]
    }
    else
      sim_kernel <- matrix(1, nrow=length(identifierNames), ncol=length(identifierNames), dimnames=list(identifierNames, identifierNames))

    sim_kernel <- sim_kernel[order(rownames(sim_kernel)), order(colnames(sim_kernel))]
    #check colnames(myBoF) == rownames|colnames(sim_kernel)
    stopifnot(all(tolower(colnames(myBoF)) == tolower(rownames(sim_kernel))))
    
    #remove classes with no identifiers, when combined with the semantic network
    
    if (!is.null(lex.fun)){
        names(identifierNames) <- identifierNames
        string_kernel <- lex.fun(identifierNames)
        string_kernel <- string_kernel[order(rownames(string_kernel)), order(colnames(string_kernel))]      
    }
    else
      string_kernel <- matrix(1, nrow=length(identifierNames), ncol=length(identifierNames))
    
    
    return(list(priori.decomp=priori.decomp, myBoF=myBoF, string_kernel=string_kernel, sim_kernel=sim_kernel))
  }
  
  library(igraph)
  library(GeLaToLab)
  library(foreach)
  library(doParallel)
  library(compiler)
  
  setwd("~/workspace")
  # Read the authoritative decomposition
  decomposition <- read.csv(paste("benchmark", prname ,"decomposition.csv", sep="/"), sep=",",  header = TRUE)
  priori.decomp <- decomposition$x
  names(priori.decomp) <- decomposition$X
  priori.decomp <- normalizeVector(priori.decomp)
  
  #Bag of Features
  myBoF <- load_BoF(prname, c(T,F)) 
  myBoF <- merge_names_by_lower_case(myBoF, 2)
  
  #Get the sample src code units
  src.code.units <- intersect(rownames(myBoF), names(priori.decomp))
  myBoF <- myBoF[src.code.units,]
  priori.decomp <- priori.decomp[src.code.units]
  
#   if (size < 1) #NOW LOOKING FOR ALL DOCS IN PACKAGES WITH 5 OR MORE ELEMENTS
  myBoF <- myBoF[get_sample_docs(prname, priori.decomp, size),]
  
  #Remove unused identifiernames 
  myBoF <- myBoF[,which(!apply(myBoF,2,FUN = function(x){all(x == 0)}))]
  
  #Filter out names shorter than 7
  identifierNames <- colnames(myBoF)
  identifierNames <- identifierNames[which(unlist(lapply(identifierNames, nchar))>=7)]
  myBoF <- myBoF[,identifierNames]

  #Remove empty classes/interfaces
  myBoF <- myBoF[which(!apply(myBoF,1,FUN = function(x){all(x == 0)})),]

 Adj <- load_SN(prname,make_symmetric = F, makeTopNode=T, identifiers=identifierNames)
    
  #populate different settings
  
  eval.funs <- c(compute_inverted_path_length, compute_Wu_Palmer_similarity, compute_Leacock_Chodorow_similarity, 
                 compute_conceptual_density, compute_Lin_similarity,  compute_Resnik_similarity )
  
  # Compile the functions used in the foreach loop
  
  compiled_eval.funs <- lapply(eval.funs, function(f) cmpfun(f))
  compiled_lex.fun <- cmpfun(normalized_LCU_kernel)
  compiled_run_each_setting <- cmpfun(run_each_setting)
  compiled_compute_semantic_similarity_clustering <- cmpfun(compute_semantic_similarity_clustering)

#TRIED AGAINST NORMALIZED IC-BASED MEASUREMENTS: NO DIFFERENCE compute_normalized_Lin_similarity & compute_normalized_Resnik_similarity
  weights <- c(1,0) #NOW ONLY CONSIDERING AS ISA (instance-of) relationship, forming a SYNSET
  
    # for (j in 1:length(compiled_eval.funs)){

  no_cores <- detectCores() - 1
  cl<-makeCluster(no_cores)
  registerDoParallel(cl)
  
  foreach(j = 1:length(compiled_eval.funs)) %dopar% {
      sim <- compiled_run_each_setting(weights, compiled_eval.funs[[j]], T, compiled_lex.fun, Adj, myBoF)
      
      #Unpack result
      priori.decomp <- sim$priori.decomp
      myBoF <- sim$myBoF
      string_kernel <- sim$string_kernel
      sim_kernel <- sim$sim_kernel
      
      #Just check to for two of the eval.funs to ensure they are working on equivalent filtered SN
      if (j %in% c(1,4)) {
        semantic <- diag(dim(sim_kernel)[2])
        r <- compiled_compute_semantic_similarity_clustering(semantic, myBoF, priori.decomp)
        print_clustering_results(prname, r, txt.file = paste("TYPE_BoF_", j, "SN_FILTERED_3_SPHERICAL_KMEANS_7_CHARS.txt", sep=""))
      }
      
      semantic <- sim_kernel
      r <- compiled_compute_semantic_similarity_clustering(semantic, myBoF, priori.decomp)
      print_clustering_results(prname, r, txt.file = paste("TYPE_ISA_", j, "SN_FILTERED_3_SPHERICAL_KMEANS_7_CHARS.txt", sep=""))
      
      semantic <- string_kernel * sim_kernel
      r <- compiled_compute_semantic_similarity_clustering(semantic, myBoF, priori.decomp)
      print_clustering_results(prname, r, txt.file = paste("TYPE_ISA_STRING_", j, "_SN_FILTERED_3_SPHERICAL_KMEANS_7_CHARS.txt", sep=""))
      
  }
  stopCluster(cl)
  

#   lex.funs <- c(normalized_LCS_kernel, normalized_LCU_kernel, constant.string.kernel)
# 
#   for (i in 1:length(lex.funs)){
#     
#     sim <- run_each_setting(c(1,0), NULL, T, lex.funs[[i]], Adj, myBoF)
# 
#     #Unpack result
#     priori.decomp <- sim$priori.decomp
#     myBoF <- sim$myBoF
#     string_kernel <- sim$string_kernel
#     sim_kernel <- sim$sim_kernel
#     
#     semantic <- string_kernel
#     r <- compute_semantic_similarity_clustering(semantic, myBoF, priori.decomp)
#     print_clustering_results(prname, r, txt.file = paste("SN_FILTERED_BY_3_STRING", i, ".txt", sep=""))
#   
#   }

}

compute_semantic_similarity_clustering <- function(semantic, myBoF, priori.decomp){
  require(skmeans)
  #SVD to compute to USU^T
  USUt <- svd(semantic)
  S <- USUt$u %*% diag(sqrt(USUt$d))
  
  #diagonal matrix for term weighings
  #TODO CHECK if this is correct
  doc.freq <- colSums(myBoF>0)
  doc.freq[doc.freq == 0] <- 1
  
#   term.freq <- rowSums(myBoF)
#   term.freq[term.freq == 0] <- 1
  
#   w <- 1/log(nrow(myBoF)/doc.freq)
  w <- log(doc.freq/nrow(myBoF))
  R <- diag(w)
  
  #Compute cosine similarity
  Phi_d <- myBoF %*% R %*% S
  
  dimnames(Phi_d) <- dimnames(myBoF)
  Phi_d <- Phi_d[order(rownames(Phi_d)),]
  
  #TODO REMOVE duplicated rows
  
  #ADD LATER IF SPHERICAL k-means doesn't work
#   kernel <- compute_cosine_kernel(Phi_d)
#   kernel <- kernel[order(rownames(kernel)), order(colnames(kernel))]
  
#   if (max(kernel) > 1)
#     stop("wrong similarity matrix!")
#   
#   #compute distance from kernel
#   dist <- 1 - kernel
#   dimnames(dist) <- dimnames(kernel)
  
  #Fix priori decomposition 
  dummy_v <- rep(0, dim(Phi_d)[1])
  names(dummy_v) <- rownames(Phi_d)
  
  priori.decomp <- find.intersection(priori.decomp, dummy_v)
  priori.decomp <- normalizeVector(priori.decomp)
  priori.decomp <- priori.decomp[order(names(priori.decomp))]
  
  if(!all(rownames(Phi_d) == names(priori.decomp)))
    stop("names don't match!")
  
  #K number of clusters
  noc <- max(priori.decomp)
  
  print("printing the numer of groups:")
  print(max(priori.decomp))
  
  
  #find the intersection with the available classes
#Spectral clustering
#   L <- laplacian(kernel, TRUE)
# clusters <- spectral.clustering(L, noc)
# names(clusters) <- rownames(L)

# Spherical K-mean
  clusters <- skmeans(Phi_d, noc, control = list(maxiter = 250, popsize=40, nruns=3))$cluster

  clusters <- normalizeVector(clusters)

  if(!all(names(clusters) == names(priori.decomp)))
    stop("names don't match!")
  
  # precision <- compute.precision(clusters, priori.decomp)
  # recall <- compute.recall(clusters, priori.decomp)
  f1.score <- compute.f1(clusters, priori.decomp)  
  adjustedRI <- compute.AdjRI(clusters, priori.decomp)
  mojosim <- compute.MoJoSim(clusters, priori.decomp)

  #   write(result, file = paste("benchmark", prname ,"DIFF.txt", sep="/"))
  
  return(list(mojosim = mojosim, f1.score=f1.score, adjustedRI=adjustedRI))
  
}

print_clustering_results <- function(prname, r, txt.file){
  setwd("~/workspace")
  
  #Prepare the result for printing to file b rounding to 3 decimal places
  print.mojosim <- round(r$mojosim, 3)
  print.f1.score <- round(r$f1.score, 3)
  print.adjustedRI <-round(r$adjustedRI, 3)
  result <- paste(print.mojosim, print.f1.score , print.adjustedRI, sep="&")
  
  write(result, file = paste("benchmark", prname, txt.file, sep="/"))
}

run_each_project_with_diffusion_kernel <- function(prname){

compute_string_semantic_diffusion_kernels <- function(prname, beta=0.5){
  require(igraph)
  require(GeLaToLab)
  setwd("~/workspace")
  # Read the authoritative decomposition
  decomposition <- read.csv(paste("benchmark", prname ,"decomposition.csv", sep="/"), sep=",",  header = TRUE)
  priori.decomp <- decomposition$x
  names(priori.decomp) <- decomposition$X
  priori.decomp <- normalizeVector(priori.decomp) 
  
  #Bag of Features
  myBoF <- load_BoF(prname, c(T,F))  
  myBoF <- merge_names_by_lower_case(myBoF, 2)
  
  #Get the sample src code units
  src.code.units <- intersect(rownames(myBoF), names(priori.decomp))
  
  myBoF <- myBoF[src.code.units,]
  priori.decomp <- priori.decomp[src.code.units]
  
  #   if (size < 1) #NOW LOOKING FOR ALL DOCS IN PACKAGES WITH 5 OR MORE ELEMENTS
  myBoF <- myBoF[get_sample_docs(prname, priori.decomp, size=0.5),]
  
  #Remove unused identifiernames
  myBoF <- myBoF[,which(!apply(myBoF,2,FUN = function(x){all(x == 0)}))]
  
  #Filter out names shorter than 5
  identifierNames <- colnames(myBoF)
  identifierNames <- identifierNames[which(unlist(lapply(identifierNames, nchar))>=4)]
  
  myBoF <- myBoF[, identifierNames]
  
  # myBoF <- myBoF[which(rowSums(myBoF) >= 4),]
  myBoF <- myBoF[which(apply(myBoF, 1, function(x) length(which(x!=0))) >= 3),]
  
  #Remove unused identifiernames
  myBoF <- myBoF[,which(!apply(myBoF,2,FUN = function(x){all(x == 0)}))]
  
  #FIXME does it make sense to pass in the identifierNames
  Adj <- load_SN(prname, make_symmetric = T, makeTopNode=T, identifiers=c())
  r <- process_All_SN(Adj, beta)
  
  K <- r$kernel
  K_names <- rownames(K)
  
  startIndex <- get.start.index.of.types(K_names)
  
  #size of dictionary
  D <- startIndex -1
  sim_kernel <- K[1:D, 1:D]
  
  #Lower case the identifer Names
  SN_identifier_names <- tolower(rownames(sim_kernel))
  BoF_identifier_names <- tolower(colnames(myBoF))
  
  #Find common identifiernames between BoF and the Semantic Network
  identifierNames <- intersect(BoF_identifier_names, SN_identifier_names)
  
  identifierIndices <- which(tolower(colnames(myBoF)) %in% identifierNames)
   
  myBoF <- myBoF[,identifierIndices]
  
  #remove classes with no identifiers, when combined with the semantic network
  
  myBoF <- myBoF[which(!apply(myBoF,1,FUN = function(x){all(x == 0)})),] 
  myBoF <- myBoF[order(rownames(myBoF)), order(colnames(myBoF))]

  sim_kernel <- sim_kernel[order(rownames(sim_kernel)), order(colnames(sim_kernel))] 
  #Eliminate identifier names in sim_kernel that are not in BoF
  identifierIndices <- which(tolower(colnames(sim_kernel)) %in% identifierNames) 
  sim_kernel <- sim_kernel[identifierIndices,identifierIndices]  
  
#   sim_kernel <- sim_kernel[identifierNames, identifierNames]   #IS THIS NECESSARY?!? NOT
#   sim_kernel <- sim_kernel[order(rownames(sim_kernel)), order(colnames(sim_kernel))] 
  
#   identifierNames <- colnames(myBoF) NOT NECESSARY
  string_kernel <- normalized_LCU_kernel(colnames(myBoF))
  string_kernel <- string_kernel[order(rownames(string_kernel)), order(colnames(string_kernel))]

  stopifnot(all(tolower(rownames(sim_kernel)) == tolower(colnames(myBoF))))
  stopifnot(all(tolower(rownames(string_kernel)) == tolower(colnames(myBoF))))
  stopifnot(all(tolower(rownames(sim_kernel)) == tolower(rownames(string_kernel))))

  return(list(priori.decomp=priori.decomp, myBoF=myBoF, string_kernel=string_kernel, sim_kernel=sim_kernel))
}

  sim <- compute_string_semantic_diffusion_kernels(prname)
  #Unpack result
  priori.decomp <- sim$priori.decomp
  myBoF <- sim$myBoF
  string_kernel <- sim$string_kernel
  sim_kernel <- sim$sim_kernel

  semantic <- diag(dim(sim_kernel)[2])
  r <- compute_semantic_similarity_clustering(semantic, myBoF, priori.decomp)
  print_clustering_results(prname, r, txt.file = "DIFF_BoF_WITH_SN_WITH_SMALL_ELIMINATED.txt")
# 
#   semantic <- sim_kernel
#   r <- compute_semantic_similarity_clustering(semantic, myBoF, priori.decomp)
#   print_clustering_results(prname, r, txt.file = "DIFF_WITH_SN_WITH_SMALL_ELIMINATED.txt")

  semantic <- string_kernel * sim_kernel
  r <- compute_semantic_similarity_clustering(semantic, myBoF, priori.decomp)
  print_clustering_results(prname, r, txt.file = "DIFF_STRING_WITH_SN.txt")
  
}

    projects <- list("apache-ant-1.9.3", "hadoop-0.20.2", "apache-log4j-1.2.17",
"jdom-2.0.5", "jedit-5.1.0", "jfreechart-1.2.0", "jhotdraw-7.0.6", "junit-4.12" ,"weka-3.6.11", "eclipse-jdt-core-3.8")