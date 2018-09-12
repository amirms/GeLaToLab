#eval.fun=compute_conceptual_density
visualize_document_type_matrix <- function(prname, eval.fun = NULL){
  library(igraph)
  library(GeLaToLab)
  
  setwd("~/workspace")
  
  # Read the authoritative decomposition
  # decomposition <- read.csv(paste("benchmark", prname ,"decomposition.csv", sep="/"), sep=",",  header = TRUE)
  # priori.decomp <- decomposition$x
  # names(priori.decomp) <- decomposition$X
  # priori.decomp <- normalizeVector(priori.decomp)
  
  #Bag of Features
  myBoF <- load_BoF(prname, c(F,T)) 
  myBoF <- merge_names_by_lower_case(myBoF, 2)
  
  #Remove unknown type
  unknownIdx <- which(colnames(myBoF) == "Unknown")
  myBoF <- myBoF[,-unknownIdx]
  
  
  print("The dimension of myBoF before eliminating small packages")
  print(dim(myBoF))
  myBoF <- myBoF[eliminate_small_packages(rownames(myBoF)),]
  print("The dimension of myBoF after eliminating small packages")
  print(dim(myBoF))
  
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
  
  Phi_d <- Phi_d[order(rownames(Phi_d)), order(colnames(Phi_d))]
  
  Rowv  <- Phi_d %>% (function(x) decostand(x, method = "normalize")) %>% dist %>% 
    hclust %>% as.dendrogram %>%
    # set("branches_k_color", k = 3) %>%
    set("branches_lwd", 2) %>%
    ladderize
  #    rotate_DendSer(ser_weight = dist(x))
  Colv  <- Phi_d %>% t %>% (function(x) decostand(x, method = "normalize")) %>% dist %>%
    hclust %>% as.dendrogram %>%
    # set("branches_k_color", k = 2) %>%
    set("branches_lwd", 2) %>%
    # reorder %>%
    ladderize
  
  library(RColorBrewer)
  # hmcol<-brewer.pal(11,"RdBu")
  
  
  heatmap(Phi_d, Rowv = Rowv, Colv = Colv, col=my_palette, xlab = "Type Packages", ylab = "Source Packages", cexCol = 1, cexRow = 1, margins = c(11, 11) )
  
}

cosine_distance <- function(x){
  x <- as.matrix(x)
  N <- nrow(x)
  c <- compute_cosine_kernel(x)
  
  1-c
}

paths <- list("org/gjt/sp/jedit/search",
              
              "org/gjt/sp/jedit/syntax",
              
              "org/gjt/sp/jedit/textarea",
              
              "org/gjt/sp/jedit/pluginmgr",
              
              c("org/gjt/sp/jedit/bsh", "org/gjt/sp/jedit/bsh/classpath"),
              
              c( "org/gjt/sp/jedit/buffer", "org/gjt/sp/jedit/bufferio", "org/gjt/sp/jedit/datatransfer"),
              
              c( "org/gjt/sp/jedit/gui", "org/gjt/sp/jedit/gui/statusbar", "org/gjt/sp/jedit/gui/tray"),
              
              "org/gjt/sp/jedit/indent",
              "org/gjt/sp/util",
              
              c( "org/gjt/sp/jedit/io", "org/gjt/sp/jedit/input")
)

type_paths <- list("java.awt", 
                   "java.beans", 
                   "java.io" , 
                   "java.lang", 
                   "java.util",
                   "javax.swing",
                   "org.xml.sax", 
                   "org.gjt.sp.jedit.search", 
                   "org.gjt.sp.jedit.syntax", 
                   "org.gjt.sp.jedit.textarea",
                   "org.gjt.sp.jedit.pluginmgr",
                   c( "org.gjt.sp.jedit.bsh", "org.gjt.sp.jedit.bsh.classpath"),
                   c( "org.gjt.sp.jedit.buffer", "org.gjt.sp.jedit.bufferset", "org.gjt.sp.jedit.bufferio", "org.gjt.sp.jedit.datatransfer"),
                   c( "org.gjt.sp.jedit.gui", "org.gjt.sp.jedit.gui.statusbar", "org.gjt.sp.jedit.gui.tray"),
                   "org.gjt.sp.jedit.indent",
                   "org/gjt/sp/util",
                   c(  "org.gjt.sp.jedit.io", "org.gjt.sp.jedit.input"))


#DONT FORGOT TO REMOVE PRIMITIVE TYPES

remove_unpresent_paths <- function(mydata, paths, dimension=1, delimiter="/") {
  
  if (dimension ==1)
    filenames <- rownames(mydata)
  else
    filenames <- colnames(mydata)
  
  paths <- unlist(paths)
  
  indices <- c()
  
  for (i in 1:length(filenames)){
    filename <- filenames[i]
    g <- grep(delimiter, strsplit(filename, "")[[1]])
    lastChar <- g[length(g)]
    path <- substr(filename, 1, lastChar-1)
    
    if (!(path %in% paths))
      indices <- c(indices, i)
  }
  
  if (dimension ==1)
    mydata <- mydata[-indices,]
  else
    mydata <- mydata[,-indices]
  
  return(mydata)
}

merge_names_by_common_prefix <- function(mydata, paths_to_include, dimension=1, delimiter = "/"){
  
  filenames <- dimnames(mydata)[[dimension]]
  
  all_paths <- lapply(filenames, function(filename) {
    g <- grep(delimiter, strsplit(filename, "")[[1]])
    lastChar <- g[length(g)]
    substr(filename, 1, lastChar-1)
  })
  
  
  duplicates = c()
  
  equivalent_indices <- lapply(paths_to_include, function(paths) which(all_paths %in% paths))
  
  for (i in 1:length(equivalent_indices))
    if (length(equivalent_indices[[i]]) > 1)
      if (!(any(equivalent_indices[[i]] %in% duplicates))){
        
        if (dimension==1){
          #add by rowsums
          added <- colSums(mydata[equivalent_indices[[i]], ])
          mydata[equivalent_indices[[i]][1],] <- added
          rownames(mydata)[equivalent_indices[[i]][1]] <- paths_to_include[[i]][1]
        }
        else{
          added <- rowSums(mydata[, equivalent_indices[[i]]])
          mydata[,equivalent_indices[[i]][1]] <- added
          colnames(mydata)[equivalent_indices[[i]][1]] <- paths_to_include[[i]][1]
        }
        
        
        duplicates <- c(duplicates, equivalent_indices[[i]][2:length(equivalent_indices[[i]])])
        
      }
  if (any(duplicates))
    if (dimension==1)
      mydata <- mydata[-duplicates,]
  else
    mydata <- mydata[,-duplicates]
  
  
  return(mydata)
  
}