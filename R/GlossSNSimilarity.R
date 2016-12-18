

convert_SN_prop_method_names_into_fulltext<- function(Adj, startIndex){
  D <- startIndex - 1
  subAdj <- Adj[startIndex:dim(Adj)[1], 1:D]
 
  identifierNames <- colnames(subAdj) 
  txts = list()
  
  for (i in 1:dim(subAdj)[1]){
    indices <- which(subAdj[i, ] > 0)
    
    txt = c()
    if (length(indices > 0)) { ####This may cause an error, since it is an empty vector
      for (j in 1:length(indices))
        txt <- c(txt, rep(identifierNames[indices[j]], subAdj[i, indices[j]]))
    
    }
    txts[[i]] <- paste(txt, collapse = " ")
  }
  names(txts) <- rownames(subAdj)
  
  return(txts)
}


convert_SN_to_Bow <- function(Adj, startIndex, nat.langs="english", prog.langs="java"){
  txts <- convert_SN_prop_method_names_into_fulltext(Adj, startIndex)
    
  mydata  <- lapply(txts, function(fulltext) strip.java.text(fulltext, lengthLowerBound=4)) # Turn into a vector of strings
  
  for (i in 1:length(prog.langs))
    mydata <- prepare.prog.lang.list(mydata, prog.langs[i])
  
  for (i in 1:length(nat.langs)) 
    mydata <- prepare.natural.lang.list(mydata, nat.langs[i])
  
  mydata.BoW.list <- make.BoW.list(mydata)
  
  mydata.BoW.frame <- make.BoW.frame(mydata.BoW.list, names(mydata.BoW.list))
  
  #Load Bag of Features, so remove the src code units which have no features in BoF --FOR COMPATIBILITY
  #   myBoF = read.csv(paste("BoF", paste(prname, "BoF.csv", sep="-"), sep="/"),  sep = ",")
  #   rownames(myBoF) <- myBoF[,1]
  #   myBoF <- myBoF[,-1]
  #   myBoF <- data.matrix(myBoF)
  #   
  #   #Remove empty source code units
  #   myBoF <- myBoF[which(!apply(myBoF,1,FUN = function(x){all(x == 0)})),]
  #   
  #   #Intersect BoW and Bof
  #   mydata.BoW.frame <- mydata.BoW.frame[rownames(myBoF),]
  #Write BoW
#   write.table(mydata.BoW.frame, file=paste("SN", "PM_BoW.csv", sep="/"),row.names=TRUE, 
#               col.names=NA,sep=",", quote=FALSE)
  
  return(mydata.BoW.frame)
}


compute_BoW_SN <- function(prname, beta=0.5){
  require(igraph)
  require(GeLaToLab)
  setwd("~/workspace")
  # Read the authoritative decomposition
  decomposition <- read.csv(paste("benchmark", prname ,"decomposition.csv", sep="/"), sep=",",  header = TRUE)
  priori.decomp <- decomposition$x
  names(priori.decomp) <- decomposition$X
  priori.decomp <- normalizeVector(priori.decomp) 
  
  Adj <- load_SN(prname,make_symmetric = T, makeTopNode=T)
  r <- process_All_SN(Adj, beta)
  
  K <- r$kernel
  
  names <- rownames(Adj)
  
  #FIND THE STARTING INDEX OF TYPE NAMES  
  classTypeIndex <- which(unlist(gregexpr(pattern ="\\.",names)) > 0)[1]
  primitiveTypeIndices <- which(names %in% c("float", "int", "char", "byte", "void", "double", "boolean"))
  startIndex <- min(c(classTypeIndex, primitiveTypeIndices))
  
  BoW <- convert_SN_to_Bow(Adj, startIndex)
  
  BoW_idf <- idf.weight(BoW)  
  
  #compute cosine similariry
  type_sim_kernel <- compute_cosine_kernel(BoW_idf)
  type_sim_kernel <- type_sim_kernel[order(rownames(type_sim_kernel)), order(colnames(type_sim_kernel))]

  type_adj <- Adj[startIndex:dim(Adj)[1], startIndex:dim(Adj)[2]]
  type_adj <- type_adj[order(rownames(type_adj)), order(colnames(type_adj))]
  
#   type_names <- rownames(type_adj)
#   
#   MU = 0.5
#   require(igraph)
#   
#   type_graph <- graph.adjacency(type_adj, mode='directed', diag=FALSE)
#   
#   for (i in 1:dim(type_adj)[1])
#     for (j in 1:dim(type_adj)[2]){
#       if (i==j)
#         break
#       depth <- compute_depth_nodes(type_graph, type_names[i], type_names[j])
#       if (is.finite(depth)){
#         hyp_sim <- MU ^ depth
#         type_sim_kernel[i,j] <- hyp_sim
#         type_sim_kernel[j,i] <- hyp_sim
#       }
#     }

  type_sim_kernel <- K[startIndex:dim(K)[1], startIndex:dim(K)[2]]
  
  ISA <- compute_ISA_SN(Adj, startIndex, type_sim_kernel)
  sim_kernel <- ISA$normalized_sim
  
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
  
  identifierNames <- rownames(Adj)[1:D]
  #Find common identifiernames between BoF and the Semantic Network
  identifierNames <- intersect(colnames(myBoF), identifierNames)
  
  #Filter out names shorter than 4
  identifierNames <- identifierNames[which(unlist(lapply(identifierNames, nchar))>4)]
  
  myBoF <- myBoF[,identifierNames]
  
  #remove classes with no identifiers, when combined with the semantic network
  
  myBoF <- myBoF[which(!apply(myBoF,1,FUN = function(x){all(x == 0)})),]
  myBoF <- myBoF[order(rownames(myBoF)), ] 
  
  sim_kernel <- sim_kernel[identifierNames, identifierNames]  
  sim_kernel <- sim_kernel[order(rownames(sim_kernel)), order(colnames(sim_kernel))] 
  
  
  #Element-wise product  
  #   semantic <- string_kernel * sim_kernel
  semantic <- sim_kernel
  
  #SVD to compute to USU^T
  USUt <- svd(semantic)
  S <- USUt$u %*% diag(sqrt(USUt$d))
  
  #diagonal matrix for term weighings
  #TODO CHECK if this is correct
  doc.freq <- colSums(myBoF>0)
  doc.freq[doc.freq == 0] <- 1
  w <- 1/log(nrow(myBoF)/doc.freq)
  R <- diag(w)  
  
  #Compute cosine similarity
  Phi_d <- myBoF %*% R %*% S
  
  dimnames(Phi_d) <- dimnames(myBoF)
  
  # return(Phi_d)
  
  kernel <- compute_cosine_kernel(Phi_d)
  kernel <- kernel[order(rownames(kernel)), order(colnames(kernel))]
  
  #Fix priori decomposition 
  dummy_v <- rep(0, dim(kernel)[1])
  names(dummy_v) <- rownames(kernel)
  
  
  priori.decomp <- find.intersection(priori.decomp, dummy_v)
  priori.decomp <- normalizeVector(priori.decomp)
  priori.decomp <- priori.decomp[order(names(priori.decomp))]  
  
  if(!all(rownames(kernel) == names(priori.decomp)))
    stop("names don't match!")
  
  #K number of clusters
  noc <- max(priori.decomp)
  
  print("printing the numer of groups:")
  print(max(priori.decomp))
  
  #find the intersection with the available classes
  
  #kmeans
  clusters <- kmeans(kernel, centers = noc, iter.max = 1500, nstart = 20000)$cluster
  clusters <- normalizeVector(clusters)
  
  # precision <- compute.precision(clusters, priori.decomp)
  # recall <- compute.recall(clusters, priori.decomp)
  f1.score <- compute.f1(clusters, priori.decomp)  
  adjustedRI <- compute.AdjRI(clusters, priori.decomp)
  mojosim <- compute.MoJoSim(clusters, priori.decomp)
    
  #Prepare the result for printing to file b rounding to 3 decimal places
  print.mojosim <- round(mojosim, 3)
  print.f1.score <- round(f1.score, 3)
  print.adjustedRI <-round(adjustedRI, 3)
  result <- paste(print.mojosim, print.f1.score , print.adjustedRI, sep="&")
  
  write(result, file = paste("benchmark", prname ,"NEW_BOW_SN2.txt", sep="/"))
  
  return(list(mojosim = mojosim, f1.score=f1.score, adjustedRI=adjustedRI))
}


compute_depth_nodes <- function(graph, root, x){
  shortest.paths(graph, v=root, to=x, mode="in")
}
