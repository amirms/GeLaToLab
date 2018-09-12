PROJECTS = list(list("apache-ant-1.9.3", "Apache Ant"),list("hadoop-0.20.2", "Apache Hadoop"),list("apache-log4j-1.2.17", "Apache Log4j"),
                list("eclipse-jdt-core-3.8","Eclipse JDT Core"), list("jdom-2.0.5","JDOM"), list("jedit-5.1.0","JEdit"),
                list("jfreechart-1.2.0","JFreeChart"), list("jhotdraw-7.0.6","JHotDraw"), list("junit-4.12","JUnit"), list("weka-3.6.11","Weka"));

prname = PROJECTS[[1]][[1]]


tabularize_information = function(projects, filename="generalInformation.txt"){
  line=""
  txt =""
  
  callsTotal=0;
  wordsTotal=0
  transTotal=0;
  
  
  for(i in 1:length(PROJECTS)) {
    setwd("~/workspace")
    project = PROJECTS[[i]]
    
    projectDirectory = project[[1]]
    projectName = project[[2]]
    
    info= extractInformation(projectDirectory);
    
    callsTotal = callsTotal +info$noOfCalls;
    wordsTotal = wordsTotal +info$noOfWords;
    transTotal = transTotal +info$noOfTransactions;
    
    line <- paste(projectName, paste("$",info$noOfCalls,"$",sep=""), paste("$",info$noOfWords,"$",sep=""), paste("$",info$noOfTransactions,"$",sep="") , sep="&")
  
    
    line <- paste(line, "\\\\ \\hline \n", sep="")
    txt <- paste(txt, line)
  }
  
  line <- paste("Average", paste("$",(callsTotal*100/length(PROJECTS)),"%$",sep=""), 
                paste("$",(wordsTotal*100/length(PROJECTS)),"%$",sep=""), 
                paste("$",(transTotal*100/length(PROJECTS)),"%$",sep="") , sep="&")
  line <- paste(line, "\\\\ \\hline \n", sep="")
  
  txt <- paste(txt, line)
  
  
  setwd("~/workspace")
  dir.create(file.path(getwd(), paste("benchmark", "Multiview/Results", sep="/")), showWarnings = FALSE)
  write(txt, file = paste("benchmark/Multiview/Results", filename, sep="/"))
  
  return(txt)
}


extractInformation = function(prname){
  setwd("~/workspace")
  #Load the adjacency matrix
  extensions= c("java/", "org/xml/", "javax/")
  cfg <- import.bunch.matrix(paste("benchmark", prname ,"dep_graph.txt", sep="/"), exclude.ext=extensions)
  #cfg <- read.table("benchmark/jedit-5.1.0/cfg.csv", sep=",", row.names = 1, header = TRUE, check.names = FALSE)
  #   cfg <- unweight.adjacency(cfg)
  # cfg <- cfg[which(rownames(cfg) %in% classnames), which(colnames(cfg) %in% classnames)]
  cfg <- cfg[order(rownames(cfg)), order(colnames(cfg))]
  cfg[cfg > 0] <- 1
  
  #Load the transaction frequency
  freq <- read.table(paste("benchmark", prname , "mydata-change-freq-matrix.csv", sep="/"), sep=",", row.names = 1, header = TRUE, check.names = FALSE)  
  freq <- as.matrix(freq)
  #freq <- freq[which(rownames(freq) %in% classnames),]
  freq <- freq[order(rownames(freq)),]
  
  freq[is.na(freq)]= 0
  freq[freq> 0] <- 1 
  
  no_transactions <- colSums(freq)
  freq <- freq[, which(no_transactions > 0)]
  
  # no_transactions <- colSums(freq)
  # freq <- freq[, which(no_transactions <= 30)]
  
  #Process the transaction frequency
  # no_transactions <- colSums(freq)
  # 
  # freq <- freq[, which(no_transactions <= 30)]
  
  
  #Load the bag of words
  BoW <- load_BoW(prname)
  # apply tf-idf mechanism and eliminate features that are lower than some threshold,
  # then remove those features from BoW
  x <- apply_tf_idf(BoW)
  dimnames(x) <- dimnames(BoW)
  BoW <- x
  
  thresh = 1
  BoW[BoW < thresh] = 0
  
  #convert BoW into a membership matrix
  BoW[BoW> 0] <- 1 
  no_words <- colSums(BoW)
  BoW <- BoW[, which(no_words > 0)]
  
  no_words_in_document <- rowSums(BoW)
  BoW <- BoW[which(no_words_in_document > 0),]
  
  
  print("dimensions before intersection")
  print(dim(cfg))
  print(dim(freq))
  print(dim(BoW))
  #INTERSECT
  names <- intersect_all(rownames(cfg), rownames(freq), rownames(BoW))
  
  cfg <- cfg[names, names]
  
  remove_empty_nodes <- function(Adj){
    #Remove nodes with no edges
    empty_rows <- which(apply(Adj,1,FUN = function(x){all(x == 0)}))
    empty_cols <- which(apply(Adj,2,FUN = function(x){all(x == 0)}))
    exclude_empty_elements <- intersect(empty_rows, empty_cols)
    
    if (length(exclude_empty_elements) > 0)
      Adj <- Adj[-exclude_empty_elements, -exclude_empty_elements]
    
    Adj
  }
  
  cfg <- remove_empty_nodes(cfg)
  names <- intersect_all(rownames(cfg), rownames(freq), rownames(BoW))
  
  
  freq <- freq[names,]
  freq <- freq[, colSums(freq) > 0]
  
  BoW <- BoW[names,]
  BoW <- BoW[, colSums(BoW) > 0]
  
  names2 <- intersect_all(rownames(cfg), rownames(freq), rownames(BoW))
  stopifnot(length(names) == length(names2))
  
  print("dimensions before intersection")
  print(dim(cfg))
  print(dim(freq))
  print(dim(BoW))
  
  return(list(noOfTransactions=dim(freq)[2], noOfWords=dim(BoW)[2], noOfCalls=sum(colSums(cfg))))
}