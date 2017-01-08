compute_information_content_similarity <- function(prname) {
  library(igraph)
  library(GeLaToLab)
  
  setwd("~/workspace")
  
  # Read the authoritative decomposition
  decomposition <- read.csv(paste("benchmark", prname ,"decomposition.csv", sep="/"), sep=",",  header = TRUE)
  priori.decomp <- decomposition$x
  names(priori.decomp) <- decomposition$X
  priori.decomp <- normalizeVector(priori.decomp)
  
  #Bag of Features
  myBoF <- load_BoF(prname, c(F,T)) 
  myBoF <- merge_names_by_lower_case(myBoF, 2)
  
  #Remove unknown type 
  # TODO move this inside load_BoF
  unknownIdx <- which(colnames(myBoF) == "Unknown")
  myBoF <- myBoF[,-unknownIdx]
  
  #Get the sample src code units
  src.code.units <- intersect(rownames(myBoF), names(priori.decomp))
  myBoF <- myBoF[src.code.units,]
  priori.decomp <- priori.decomp[src.code.units]
  
  #   if (size < 1) #NOW LOOKING FOR ALL DOCS IN PACKAGES WITH 5 OR MORE ELEMENTS
  myBoF <- myBoF[get_sample_docs(prname, priori.decomp, size),]
  
  
  #LOAD Semantic Network
  Adj <- load_SN(prname, make_symmetric = F, makeTopNode=T, identifiers=c())
  names <- colnames(Adj)
  startIndex <- get.start.index.of.types(names)
  
  S <- Adj[startIndex:dim(Adj)[1], startIndex:dim(Adj)[2]]
  # dimnames(C) <- dimnames(Adj[startIndex:dim(Adj)[1], 1:(startIndex-1)])
  
  #DONE make this a higher function argument
  lin_type_sim <- compute_Lin_similarity(S, myBoF)
  resnik_type_sim <- compute_Resnik_similarity(S, myBoF)
  
  stopifnot(all(colnames(lin_type_sim) == colnames(resnik_type_sim)))
  
  # order row and column names
  lin_type_sim <- lin_type_sim[order(rownames(lin_type_sim)), order(colnames(lin_type_sim))]
  resnik_type_sim <- resnik_type_sim[order(rownames(resnik_type_sim)), order(colnames(resnik_type_sim))]
  
    #Create results directory, if it doesn't exist
  dir.create(file.path(getwd(), paste("benchmark", prname, "Results/InformationContent", sep="/")), showWarnings = FALSE)
  
  old_priori.decomp <- priori.decomp
  
  compute_information_content_semantic_similarity_clustering(priori.decomp, lin_type_sim, "Lin_Type_Sim", prname)
  compute_information_content_semantic_similarity_clustering(priori.decomp, resnik_type_sim, "Resnik_Type_Sim.txt", prname)
}

compute_information_content_semantic_similarity_clustering <- function(priori.decomp, type_sim, fileName, prname){
  empty_rows <- which(apply(type_sim,1,FUN = function(x){all(x == 0)}))
  empty_cols <- which(apply(type_sim,2,FUN = function(x){all(x == 0)}))
  toEliminate <- intersect(empty_rows, empty_cols)
  
  if (length(toEliminate) > 0)
    type_sim <- type_sim[-toEliminate, -toEliminate]
  
  #Fix priori decomposition
  priori.decomp <- priori.decomp[colnames(type_sim)]
  priori.decomp <- priori.decomp[order(names(priori.decomp))]
  priori.decomp <- normalizeVector(priori.decomp)
  
  if(!all(rownames(type_sim) == names(priori.decomp)))
    stop("names don't match!")
  
  r <- compute_spectral_clustering(type_sim, priori.decomp)
  print_clustering_results(prname, r, txt.file = paste("Results/InformationContent/", fileName ,".txt", sep=""))
}

#Input: r: similarity matrix, priori.decomp
compute_spectral_clustering <- function(r, priori.decomp){
  #K number of clusters
  noc <- max(priori.decomp)
  
  L <- normalized.symmetric.laplacian(r)
  
  partition <- spectral.clustering(L, noc)
  partition <- normalizeVector(partition)
  
  names(partition) <- colnames(r)
  
  f1.score <- compute.f1(partition, priori.decomp)  
  adjustedRI <- compute.AdjRI(partition, priori.decomp)
  mojosim <- compute.MoJoSim(partition, priori.decomp)
  
  results <- list(mojosim = mojosim, f1.score=f1.score, adjustedRI=adjustedRI)
  return(results)
}

#Input: 
# S is type hypernymy graph
# C is type content
compute_Resnik_similarity <- function(S, BoT){
  
  stopifnot(all(rownames(S) == colnames(S)))
  
  sum_BoT <- colSums(BoT)
  
  information_content <- -log(sum_BoT /sum(sum_BoT))
  
  information_content <- unlist(lapply(information_content, function(ic) if (is.infinite(ic)) 0 else ic))
  
  noOfModules <- dim(BoT)[1]
  moduleNames <- rownames(BoT)
  
  sim_matrix<- matrix(0, nrow=noOfModules, ncol=noOfModules, dimnames = list(moduleNames, moduleNames))
  
  icNames <- names(information_content)
  
  typeNames <-rownames(S)
  
  for (i in 1:noOfModules) {
    #FOR NOW, FIX LATER
    typeNameIndex_i <- which(typeNames == convertModule2TypeName(moduleNames[i]))
    
    if (length(typeNameIndex_i) != 1)
      next
    
    for(j in i:noOfModules)
      if (i != j)
      {    
        #check if one is subtype of other
        typeNameIndex_j <- which(typeNames ==  convertModule2TypeName(moduleNames[j]))
        
        if (length(typeNameIndex_j) != 1)
          next
        
        LCHs <- compute_least_common_hypernym(S, typeNameIndex_i, typeNameIndex_j)
        
        if (length(LCHs) < 1)
          next
        
        maxSim <- 0
        
        for (k in 1 : length(LCHs)){
          
          LCH <- LCHs[k]
          
          LCH_name <- typeNames[LCH]
          
          if (!(LCH_name %in% icNames))
            next
          
          sim <- information_content[LCH_name]
          
          if (sim > maxSim)
            maxSim <- sim       
        } 
        
        sim_matrix[i,j] <- maxSim       
      }
  }
  
  sim_matrix <- fill_lower_diagonal(sim_matrix)  
  
  sim_matrix  
}

convertModule2TypeName <- function(moduleName){
  g <- regexpr("\\.[^\\.]*$", moduleName)
  
  if (g <= 0)
    return(moduleName)
  
  #Found a dot character
  typeName <- substr(moduleName, start=1, stop=g-1)
  typeName <- gsub("/", "\\.", typeName)
  return(typeName)
}

#Input: 
# S is type hypernymy graph
# BoT is the bag of types extracted from the corpu of the source code 
compute_Lin_similarity <- function(S, BoT){
  stopifnot(all(rownames(S) == colnames(S)))
  
  sum_BoT <- colSums(BoT)
  
  information_content <- -log(sum_BoT /sum(sum_BoT))
  
  information_content <- unlist(lapply(information_content, function(ic) if (is.infinite(ic)) 0 else ic))
  
  noOfModules <- dim(BoT)[1]
  moduleNames <- rownames(BoT)
  
  sim_matrix<- matrix(0, nrow=noOfModules, ncol=noOfModules, dimnames = list(moduleNames, moduleNames))
  
  icNames <- names(information_content)
  
  typeNames <-rownames(S)
  
  for (i in 1:noOfModules) {
    #FOR NOW, FIX LATER
    typeNameIndex_i <- which(typeNames == convertModule2TypeName(moduleNames[i]))
    
    if (length(typeNameIndex_i) != 1)
      next
    
    for(j in i:noOfModules)
      if (i != j)
      {    
        #check if one is subtype of other
        typeNameIndex_j <- which(typeNames ==  convertModule2TypeName(moduleNames[j]))
        
        if (length(typeNameIndex_j) != 1)
          next
        
        LCHs <- compute_least_common_hypernym(S, typeNameIndex_i, typeNameIndex_j)
        
        if (length(LCHs) < 1)
          next
        
        maxSim <- 0
        
        for (k in 1 : length(LCHs)){
          
          LCH <- LCHs[k]
          LCH_name <- typeNames[LCH]
          if (!(LCH_name %in% icNames))
            next
          
          i_name <- typeNames[typeNameIndex_i]
          if (!(i_name %in% icNames))
            next
          
          j_name <- typeNames[typeNameIndex_j]
          if (!(j_name %in% icNames))
            next
          
          
          sim <- 2 * information_content[LCH_name] / (information_content[i_name] + information_content[j_name])
          
          if (sim > maxSim)
            maxSim <- sim
          
        }        
        sim_matrix[i,j] <- maxSim
        
      }
  }
  
  sim_matrix <- fill_lower_diagonal(sim_matrix)  

  sim_matrix    
}
