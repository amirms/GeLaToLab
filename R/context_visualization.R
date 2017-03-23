
#eval.fun=compute_conceptual_density
visualize_document_type_matrix <- function(prname, eval.fun = NULL){
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
  myBoF <- myBoF[get_sample_docs(prname, priori.decomp, 0.5),]
  
  if (!is.null(eval.fun)){
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
  }
  
  # if (is.null(eval.fun))
    semantic <- diag(dim(myBoF)[2])
  # else
  #   semantic <- type_sim
  
  USUt <- svd(semantic)
  S <- USUt$u %*% diag(sqrt(USUt$d))
  
  Phi_d <- apply_tf_idf(myBoF) %*% S
  
  dimnames(Phi_d) <- dimnames(myBoF)
  Phi_d <- Phi_d[order(rownames(Phi_d)),]
  
  print("Dimensions of Phi_d before Cleansing")
  print(dim(Phi_d))
  
  #Remove empty rows
  Phi_d <- Phi_d[ apply(Phi_d!=0, 1, any), , drop=FALSE]
  Phi_d <- Phi_d[!duplicated(Phi_d),]
  Phi_d <- Phi_d[ , apply(Phi_d!=0, 2, any), drop=FALSE] 
  
  print("Dimensions of Phi_d after Cleansing")
  print(dim(Phi_d))
  
  {
  indices <- which(unlist(lapply(colnames(Phi_d), function(name) grep("jedit",name))) > 0)
  
    Y <- Phi_d
    X <-  svd(Phi_d)
    D <- X$d
    D[11:length(D)] <- 0
    Phi_d <- X$u %*% diag(D) %*% t(X$v)
    dimnames(Phi_d) <- dimnames(Y)
    
    packages <- colnames(Y)
    
    # TODO For Document-Context matrix, context vector representation, 
    # make a document vs content heatmap, with distribution of contexts in documents
    # group together types based on packages (probably means eliminating primitive types)
    # also group together modules based on packages
    # explain the usage and the relation between packages
  }
  
  library(dendextend)
  library(vegan)
  
  Rowv  <- Phi_d %>% (function(x) decostand(x, method = "normalize")) %>% dist %>% 
    hclust %>% as.dendrogram %>%
    set("branches_k_color", k = 3) %>% set("branches_lwd", 4) %>%
    ladderize
  #    rotate_DendSer(ser_weight = dist(x))
  Colv  <- Phi_d %>% t %>% (function(x) decostand(x, method = "normalize")) %>% dist %>%
    hclust %>% as.dendrogram %>%
    set("branches_k_color", k = 2) %>% set("branches_lwd", 4) %>%
    ladderize
    
  
  
  
    library(d3heatmap)
  heatmap(Phi_d, Rowv = Rowv, Colv = Colv)
  
}

cosine_distance <- function(x){
  x <- as.matrix(x)
  N <- nrow(x)
  c <- compute_cosine_kernel(x)
  
  1-c
}