# projects <- list("apache-ant-1.9.3", "hadoop-0.20.2", "apache-log4j-1.2.17", "jdom-2.0.5", "jedit-5.1.0",
#                  "jfreechart-1.2.0", "jhotdraw-7.0.6", "junit-4.12" ,"weka-3.6.11", "eclipse-jdt-core-3.8")

run_context_type_clustering <- function(projects){
  library(foreach)
  library(doParallel)
  library(GeLaToLab)
  
  setwd("~/workspace")
  
  no_cores <- detectCores() - 1
  # no_cores <- 3
  cl<-makeCluster(no_cores)
  registerDoParallel(cl)
  
  foreach(i = 1:length(projects)) %dopar% {
    GeLaToLab::compute_type_similarity_clustering(projects[[i]])
  }
  
  stopCluster(cl)
}


run_context_type_identifier_clustering <- function(projects){
  library(foreach)
  library(doParallel)
  library(GeLaToLab)
  
  setwd("~/workspace")
  
  no_cores <- detectCores()
  # no_cores <- 3
  cl<-makeCluster(no_cores)
  registerDoParallel(cl)
  
  foreach(i = 1:length(projects)) %dopar% {
    GeLaToLab::compute_identifier_type_similarity_clustering(projects[[i]])
  }
  
  stopCluster(cl)
}


compute_type_similarity_clustering <- function(prname, eval.fun = compute_conceptual_density){
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
  C <- Adj[startIndex:dim(Adj)[1], 1:(startIndex-1)] 
  # dimnames(C) <- dimnames(Adj[startIndex:dim(Adj)[1], 1:(startIndex-1)])
  
  #DONE make this a higher function argument
  type_sim <- eval.fun(S, C)
  dimnames(type_sim) <- dimnames(S)

  #Remove unused type names 
  myBoF <- myBoF[,which(!apply(myBoF,2,FUN = function(x){all(x == 0)}))]
 
  common_types <- intersect(rownames(type_sim), colnames(myBoF))
  
  #FIX AND ORDER TYPE_SIM 
  type_sim <- type_sim[common_types, common_types]
  type_sim <- type_sim[order(rownames(type_sim)),order(colnames(type_sim)) ]
  
  #FIX AND ORDER MyBoF
  #Remove empty classes/interfaces
  myBoF <- myBoF[,common_types]
  myBoF <- myBoF[which(!apply(myBoF,1,FUN = function(x){all(x == 0)})),]
  myBoF <- myBoF[order(rownames(myBoF)),order(colnames(myBoF)) ]
  
  #Create results directory, if it doesn't exist
  dir.create(file.path(getwd(), paste(prname, "Results/ContextModel", sep="/")), showWarnings = FALSE)
  
  #compute_semantic_similarity_clustering
  semantic <- diag(dim(myBoF)[2])
  r <- compute_semantic_similarity_clustering(semantic, myBoF, priori.decomp)
  print_clustering_results(prname, r, txt.file = "Results/ContextModel/Latest_Type_BoF_Sim.txt")
  
  semantic <- type_sim
  r <- compute_semantic_similarity_clustering(semantic, myBoF, priori.decomp)
  print_clustering_results(prname, r, txt.file = "Results/ContextModel/Latest_Enriched_Type_BoF_Sim.txt")
  
  return(r)
}

compute_identifier_type_similarity_clustering <- function(prname, eval.fun=compute_conceptual_density, lex.fun=normalized_LCU_kernel){
  library(igraph)
  library(GeLaToLab)
  
  setwd("~/workspace")
  
  # Read the authoritative decomposition
  decomposition <- read.csv(paste("benchmark", prname ,"decomposition.csv", sep="/"), sep=",",  header = TRUE)
  priori.decomp <- decomposition$x
  names(priori.decomp) <- decomposition$X
  priori.decomp <- normalizeVector(priori.decomp)
  
  #Bag of Features
  myBoF_data <- load_BoF(prname, c(T,T)) 
  myBoF <- myBoF_data$myBoF
  myBoF <- merge_names_by_lower_case(myBoF, 2)
  
  identifiers <- myBoF_data$identifiers

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
  C <- Adj[startIndex:dim(Adj)[1], 1:(startIndex-1)] 

  #DONE make this a higher function argument
  type_sim <- eval.fun(S, C)
  colnames(type_sim) <- tolower(colnames(S))
  rownames(type_sim) <- tolower(rownames(S))
  diag(type_sim) <- 1
  
  stopifnot(isSymmetric((type_sim)))
  
  # Compute identifier similarity
  unique_identifier_names <- unique(tolower(identifiers))
  identifier_sim <- lex.fun(unique_identifier_names)
  diag(identifier_sim) <- 1
  
  stopifnot(isSymmetric((identifier_sim)))
  
  #Remove unused identifier_type names 
  myBoF <- myBoF[,which(!apply(myBoF,2,FUN = function(x){all(x == 0)}))]
  #Remove empty classes/interfaces
  myBoF <- myBoF[which(!apply(myBoF,1,FUN = function(x){all(x == 0)})),]
  
  # Fix identifier type similarity
  used_identifier_types <- tolower(colnames(myBoF))
  no_used_identifier_types <- length(used_identifier_types)
  
  identifier_type_sim = matrix(0, nrow=no_used_identifier_types, ncol = no_used_identifier_types,
                               dimnames=list(used_identifier_types, used_identifier_types))
  
  all_type_names <- colnames(type_sim)
  all_identifier_names <- colnames(identifier_sim)
  
  for (i in 1:no_used_identifier_types){
    splitted_string <- unlist(strsplit(used_identifier_types[i], ":"))
    current_identifier <- splitted_string[1]
    current_type <- splitted_string[2]
    
    # print("RoW:")
    # print(i)
    
    if (!(current_type %in% all_type_names) || !(current_identifier %in% all_identifier_names) )
      next
    
    for (j in i:no_used_identifier_types)
      if (i != j) {
        
        # print("Col:")
        # print(j)
        
      splitted_string <- unlist(strsplit(used_identifier_types[j], ":"))
      identifier <- splitted_string[1]
      type <- splitted_string[2]
      
      if (!(type %in% all_type_names) || !(identifier %in% all_identifier_names) )
        next

      
      
      identifier_type_sim[i,j] <- type_sim[current_type, type] * identifier_sim[current_identifier, identifier]
      }
  }
  identifier_type_sim <- fill_lower_diagonal(identifier_type_sim)
  diag(identifier_type_sim) <- 1
  
  #Create results directory, if it doesn't exist
  dir.create(file.path(getwd(), paste(prname, "Results/ContextModel", sep="/")), showWarnings = FALSE)
  
  #compute_semantic_similarity_clustering
  semantic <- diag(dim(myBoF)[2])
  r <- compute_semantic_similarity_clustering(semantic, myBoF, priori.decomp)
  print_clustering_results(prname, r, txt.file = "Results/ContextModel/Identifier_Type_BoF_Sim.txt")
  
  semantic <- identifier_type_sim
  r <- compute_semantic_similarity_clustering(semantic, myBoF, priori.decomp)
  print_clustering_results(prname, r, txt.file = "Results/ContextModel/Enriched_Identifier_Type_BoF_Sim.txt")
  
  return(r)
}