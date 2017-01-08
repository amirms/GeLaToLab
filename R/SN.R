prepare.Bow <- function(prname, rootFolder="org", pattern = "*.java", prog.langs=c("java"), nat.langs=c("english")){
  require(igraph)
  require(GeLaToLab)
  setwd("~/workspace")
  setwd(paste("benchmark", prname, sep="/"))
  mydata <- read.directory(rootFolder, pattern)
  
  ################
  # Preprocessing data #
  ################
  
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
  write.table(mydata.BoW.frame, file=paste("BoW", paste(prname, "BoW.csv", sep="-"), sep="/"),row.names=TRUE, 
              col.names=NA,sep=",", quote=FALSE)
  
  return (mydata.BoW.frame)
}


compute.BoW.kernel <- function(prname, size = 0.25, rootFolder="org", pattern = "*.java", prog.langs=c("java"), nat.langs=c("english")){
  require(igraph)
  require(GeLaToLab)
  setwd("~/workspace")
  
  #Load theauthoritative decomposition
  decomposition <- read.csv(paste("benchmark", prname ,"decomposition.csv", sep="/"), sep=",",  header = TRUE)
  priori.decomp <- decomposition$x
  names(priori.decomp) <- decomposition$X
  priori.decomp <- normalizeVector(priori.decomp)
  
  
  if (file.exists(paste("benchmark", prname , "BoW", paste(prname, "BoW.csv", sep="-"), sep="/")))  
      mydata.BoW.frame <- read.table(file=paste("benchmark", prname , "BoW", paste(prname, "BoW.csv", sep="-"), sep="/"), 
                                     sep=",", row.names = 1, header = TRUE, check.names = FALSE)
  else
     mydata.BoW.frame <- prepare.Bow(prname, rootFolder, pattern, prog.langs, nat.langs)
  
  
    myBoF = read.csv(paste("benchmark", prname , "BoF", paste(prname, "BoF.csv", sep="-"), sep="/"),  sep = ",")
    rownames(myBoF) <- myBoF[,1]
    myBoF <- myBoF[,-1]
    myBoF <- data.matrix(myBoF)
  
  #Get the sample src code units
  src.code.units <- intersect(rownames(myBoF), names(priori.decomp))
  myBoF <- myBoF[src.code.units,]
  priori.decomp <- priori.decomp[src.code.units]
  
  #Get the sample src code units
  if (size < 1)
    myBoF <- myBoF[get_sample_docs(prname, priori.decomp, size),]
  
  #Remove unused identifiernames
  myBoF <- myBoF[,which(!apply(myBoF,2,FUN = function(x){all(x == 0)}))]
  
  #Remove empty source code units
  myBoF <- myBoF[which(!apply(myBoF,1,FUN = function(x){all(x == 0)})),]

  #Intersect BoW and Bof
  mydata.BoW.frame <- mydata.BoW.frame[rownames(myBoF),]
  
  #test the idf.weight function to see if it matches 
  mydata.BoW.idf.frame <- idf.weight(mydata.BoW.frame)  
  
  #compute cosine similariry
  
  #   mydata.BoW.idf.frame.t <- t(na.omit(mydata.BoW.idf.frame))
  
  kernel <- compute_cosine_kernel(mydata.BoW.idf.frame)
  
  #Fix priori decomposition 
  dummy_v <- rep(0, dim(kernel)[1])
  names(dummy_v) <- rownames(kernel)
  
  priori.decomp <- find.intersection(priori.decomp, dummy_v)
  priori.decomp <- normalizeVector(priori.decomp)
  priori.decomp <- priori.decomp[order(names(priori.decomp))]
  
  #K number of clusters
  noc <- max(priori.decomp)
  
  print("printing the numer of groups:")
  print(max(priori.decomp))
  
  
  #find the intersection with the available classes
  
  #kmeans
  clusters <- kmeans(kernel, centers = noc, iter.max = 1500, nstart = 20000)$cluster
#   print(clusters$)
  # result <- spectral.clustering(addition.kernel, k, iter = 15)
  # result <- spectral.clustering(laplacian(addition.kernel, TRUE), k, iter = 40)
  
  clusters <- normalizeVector(clusters)
  
  names(clusters) = rownames(kernel)
  
  precision <- compute.precision(clusters, priori.decomp)
  recall <- compute.recall(clusters, priori.decomp)
  f1.score <- compute.f1(clusters, priori.decomp)
  
  mojosim <- compute.MoJoSim(clusters, priori.decomp)
  
  print(paste(round(mojosim, 3),round(precision, 3) ,round(recall, 3), round(f1.score, 3), sep="&"))
  
  
  return(list(mojosim=mojosim, precision=precision, recall=recall, f1.score=f1.score))
  
  
  
  
}

compute.BoF.kernel <- function(prname, size=0.25){
  require(igraph)
  require(GeLaToLab)
  setwd("~/workspace")
  # Read the authoritative decomposition
  decomposition <- read.csv(paste("benchmark", prname ,"decomposition.csv", sep="/"), sep=",",  header = TRUE)
  priori.decomp <- decomposition$x
  names(priori.decomp) <- decomposition$X
  priori.decomp <- normalizeVector(priori.decomp)
  
  #Make BoW from BoF
  #Load Bag of Features
  myBoF = read.csv(paste("benchmark", prname , "BoF", paste(prname, "BoF.csv", sep="-"), sep="/"),  sep = ",")
  rownames(myBoF) <- myBoF[,1]
  myBoF <- myBoF[,-1]
  myBoF <- data.matrix(myBoF)
  
  #Get the sample src code units
  src.code.units <- intersect(rownames(myBoF), names(priori.decomp))
  myBoF <- myBoF[src.code.units,]
  priori.decomp <- priori.decomp[src.code.units]
  
  #Get the sample src code units
  if (size < 1)
    myBoF <- myBoF[get_sample_docs(prname, priori.decomp, size),]
  
  #Remove unused identifiernames
  myBoF <- myBoF[,which(!apply(myBoF,2,FUN = function(x){all(x == 0)}))]
  
  r <- load_process_All_SN(prname, 0.5)
  
  #Find common identifiernames between BoF and the Semantic Network
  identifierNames <- intersect(colnames(myBoF), colnames(r$kernel))
  
  #Filter out names shorter than 4
  identifierNames <- identifierNames[which(unlist(lapply(identifierNames, nchar))>4)]
  
  
  r$kernel <- r$kernel[identifierNames, identifierNames]
  myBoF <- myBoF[,identifierNames]
  
  #remove classes with no identifiers, when combined with the semantic network    
  myBoF <- myBoF[which(!apply(myBoF,1,FUN = function(x){all(x == 0)})),]
  
  
  kernel <- compute_cosine_kernel(myBoF)
  
  #Fix priori decomposition 
  dummy_v <- rep(0, dim(kernel)[1])
  names(dummy_v) <- rownames(kernel)
  
  priori.decomp <- find.intersection(priori.decomp, dummy_v)
  priori.decomp <- normalizeVector(priori.decomp)
  priori.decomp <- priori.decomp[order(names(priori.decomp))]
  
  #K number of clusters
  noc <- max(priori.decomp)
  
  print("printing the numer of groups:")
  print(max(priori.decomp))
  
  
  #find the intersection with the available classes
  
  #kmeans
  clusters <- kmeans(kernel, centers = noc, iter.max = 1500, nstart = 20000)$cluster
  # result <- spectral.clustering(addition.kernel, k, iter = 15)
  # result <- spectral.clustering(laplacian(addition.kernel, TRUE), k, iter = 40)
  
  clusters <- normalizeVector(clusters)
  
  names(clusters) = rownames(kernel)
  
#   precision <- compute.precision(clusters, priori.decomp)
#   recall <- compute.recall(clusters, priori.decomp)
  f1.score <- compute.f1(clusters, priori.decomp)

adjustedRI <- compute.AdjRI(clusters, priori.decomp)

  
  mojosim <- compute.MoJoSim(clusters, priori.decomp)
  
  print(paste(round(mojosim, 3),round(f1.score, 3) ,round(adjustedRI, 3), sep="&"))
  
  return(list(mojosim=mojosim, f1.score=f1.score, adjustedRI=adjustedRI))
}


compute.BoF.BoW.kernel <- function(prname, size=0.25){

    require(igraph)
    require(GeLaToLab)
    setwd("~/workspace")
    # Read the authoritative decomposition
    decomposition <- read.csv(paste("benchmark", prname ,"decomposition.csv", sep="/"), sep=",",  header = TRUE)
    priori.decomp <- decomposition$x
    names(priori.decomp) <- decomposition$X
    priori.decomp <- normalizeVector(priori.decomp)
  
    #Make BoW from BoF
    #Load Bag of Features
    myBoF = read.csv(paste("benchmark", prname , "BoF", paste(prname, "BoF.csv", sep="-"), sep="/"),  sep = ",")
    rownames(myBoF) <- myBoF[,1]
    myBoF <- myBoF[,-1]
    myBoF <- data.matrix(myBoF)
    
    #Get the sample src code units
    src.code.units <- intersect(rownames(myBoF), names(priori.decomp))
    myBoF <- myBoF[src.code.units,]
    priori.decomp <- priori.decomp[src.code.units]
    
    #Get the sample src code units
    if (size < 1)
      myBoF <- myBoF[get_sample_docs(prname, priori.decomp, size),]
    
    #Remove unused identifiernames
    myBoF <- myBoF[,which(!apply(myBoF,2,FUN = function(x){all(x == 0)}))]
    
    r <- load_process_All_SN(prname, beta)
    
    #Find common identifiernames between BoF and the Semantic Network
    identifierNames <- intersect(colnames(myBoF), colnames(r$kernel))
    
    #Filter out names shorter than 4
    identifierNames <- identifierNames[which(unlist(lapply(identifierNames, nchar))>4)]
    
    
    r$kernel <- r$kernel[identifierNames, identifierNames]
    myBoF <- myBoF[,identifierNames]
    
    #remove classes with no identifiers, when combined with the semantic network    
    myBoF <- myBoF[which(!apply(myBoF,1,FUN = function(x){all(x == 0)})),]

    
#     #Find common identifiernames between BoF and the Semantic Network
    identifierNames <- colnames(myBoF)
    stcUnits <- rownames(myBoF)

  txts <- list()

  for (i in 1:dim(myBoF)[1]){
    txt <- ""
    for (j in 1:dim(myBoF)[2]){
      if (myBoF[i, j] > 0)
        for (z in 1: myBoF[i,j])
          txt <- paste(txt, identifierNames[j], sep= " ")
      
    }
    
    txts[[stcUnits[i]]]<- txt
    
  }
    
  txts <- lapply(txts, strip.java.text, 4)

  txts <- prepare.prog.lang.list(txts, "java")

  txts <- prepare.natural.lang.list(txts, "english")



  mydata.BoW.list <- make.BoW.list(txts)

  mydata.BoW.frame <- make.BoW.frame(mydata.BoW.list, names(mydata.BoW.list))

  print(dim(mydata.BoW.frame))

  write.table(mydata.BoW.frame, file=paste("benchmark", prname , "BoF", paste(prname, "BoW.csv", sep="-"), sep="/"),row.names=TRUE, 
                                           col.names=NA,sep=",", quote=FALSE)

  #test the idf.weight function to see if it matches 
  mydata.BoW.idf.frame <- idf.weight(mydata.BoW.frame)  
  
  #compute cosine similariry
  
#   mydata.BoW.idf.frame.t <- t(na.omit(mydata.BoW.idf.frame))
  
  kernel <- compute_cosine_kernel(mydata.BoW.idf.frame)

  #Fix priori decomposition 
  dummy_v <- rep(0, dim(kernel)[1])
  names(dummy_v) <- rownames(kernel)
  
  priori.decomp <- find.intersection(priori.decomp, dummy_v)
  priori.decomp <- normalizeVector(priori.decomp)
  priori.decomp <- priori.decomp[order(names(priori.decomp))]
  
  #K number of clusters
  noc <- max(priori.decomp)
  
  print("printing the numer of groups:")
  print(max(priori.decomp))
  
  
  #find the intersection with the available classes
  
  #kmeans
  clusters <- kmeans(kernel, centers = noc, iter.max = 1500, nstart = 20000)$cluster
  # result <- spectral.clustering(addition.kernel, k, iter = 15)
  # result <- spectral.clustering(laplacian(addition.kernel, TRUE), k, iter = 40)
  
  clusters <- normalizeVector(clusters)
  
  names(clusters) = rownames(kernel)
  
  precision <- compute.precision(clusters, priori.decomp)
  recall <- compute.recall(clusters, priori.decomp)
  f1.score <- compute.f1(clusters, priori.decomp)
  
  mojosim <- compute.MoJoSim(clusters, priori.decomp)

print(paste(round(mojosim, 3),round(precision, 3) ,round(recall, 3), round(f1.score, 3), sep="&"))

return(list(mojosim=mojosim, precision=precision, recall=recall, f1.score=f1.score))
  
}



compute.semantic.kernel <- function(prname, choice="PR", DIFF = 1, beta=0.5, size= 0.25, ISA_IPO_choice=0, ISA_IPO= c(0.5, 0.5), rootFolder = "org", baseline){
  require(igraph)
  require(GeLaToLab)
setwd("~/workspace")
# Read the authoritative decomposition
decomposition <- read.csv(paste("benchmark", prname ,"decomposition.csv", sep="/"), sep=",",  header = TRUE)
priori.decomp <- decomposition$x
names(priori.decomp) <- decomposition$X
priori.decomp <- normalizeVector(priori.decomp)


if (ISA_IPO_choice > 0){
  print("reached ISA_IPO")
    r <- load_process_ISA_IPO_SN(prname, beta, ISA_IPO, rootFolder)

}
else if (DIFF > 0){
  print("reached ALL DIFF")
  r <- load_process_All_SN(prname, beta)
  
}
# else  
# {
#   print("reached ALL CT")
#   r <- load_process_All_CT_SN(prname)
# 
# }

#Bag of Features
myBoF = read.csv(paste("benchmark", prname , "BoF", paste(prname, "BoF.csv", sep="-"), sep="/"),  sep = ",")
rownames(myBoF) <- myBoF[,1]
myBoF <- myBoF[,-1]
myBoF <- data.matrix(myBoF)

#Get the sample src code units
src.code.units <- intersect(rownames(myBoF), names(priori.decomp))
myBoF <- myBoF[src.code.units,]
priori.decomp <- priori.decomp[src.code.units]

if (size < 1)
  myBoF <- myBoF[get_sample_docs(prname, priori.decomp, size),]

#Remove unused identifiernames

myBoF <- myBoF[,which(!apply(myBoF,2,FUN = function(x){all(x == 0)}))]

#Find common identifiernames between BoF and the Semantic Network
identifierNames <- intersect(colnames(myBoF), colnames(r$kernel))

#Filter out names shorter than 4
identifierNames <- identifierNames[which(unlist(lapply(identifierNames, nchar))>4)]

SK <- constant.string.kernel(identifierNames)
dimnames(SK) <- list(identifierNames, identifierNames) 


r$kernel <- r$kernel[identifierNames, identifierNames]
r$dist <- r$dist[identifierNames, identifierNames]
myBoF <- myBoF[,identifierNames]

#remove classes with no identifiers, when combined with the semantic network

myBoF <- myBoF[which(!apply(myBoF,1,FUN = function(x){all(x == 0)})),]


#order the names
# myBoF <- myBoF[order(rownames(myBoF)), ]
# r$kernel <- r$kernel[order(rownames(r$kernel)), order(colnames(r$kernel))]
# SK <- SK[order(rownames(SK)), order(colnames(SK))]

#sorting
# sort(r$kernel[1,], decreasing = T)[1:5]

semantic <- switch(choice,
          SK = SK,
          SN = r$kernel,
          PR = SK * r$kernel)


#SVD to compute to USU^T
USUt <- svd(semantic)
S <- USUt$u %*% diag(sqrt(USUt$d))

#LSA by reducing concepts
# d <- rep(0, nrow(USUt$u))
# D <- USUt$d[USUt$d > 0.7]
# d[1:10] <- USUt$d[1:10]
# semantic <- USUt$u %*% diag(d) %*% t(USUt$u)

#diagonal matrix for term weighings
#TODO CHECK if this is correct
doc.freq <- colSums(myBoF>0)
doc.freq[doc.freq == 0] <- 1
w <- 1/log(nrow(myBoF)/doc.freq)
R <- diag(w)


#Compute cosine similarity
Phi_d <- myBoF %*% R %*% S

dimnames(Phi_d) <- dimnames(myBoF)



kernel <- compute_cosine_kernel(Phi_d)
kernel <- kernel[order(rownames(kernel)), order(colnames(kernel))]




#Fix priori decomposition 
dummy_v <- rep(0, dim(kernel)[1])
names(dummy_v) <- rownames(kernel)

priori.decomp <- find.intersection(priori.decomp, dummy_v)
priori.decomp <- normalizeVector(priori.decomp)
priori.decomp <- priori.decomp[order(names(priori.decomp))]

#K number of clusters
noc <- max(priori.decomp)

print("printing the numer of groups:")
print(max(priori.decomp))


#find the intersection with the available classes

#kmeans
clusters <- kmeans(kernel, centers = noc, iter.max = 1500, nstart = 20000)$cluster
# result <- spectral.clustering(addition.kernel, k, iter = 15)
# result <- spectral.clustering(laplacian(addition.kernel, TRUE), k, iter = 40)

clusters <- normalizeVector(clusters)

names(clusters) = rownames(kernel)

# precision <- compute.precision(clusters, priori.decomp)
# recall <- compute.recall(clusters, priori.decomp)
f1.score <- compute.f1(clusters, priori.decomp)

adjustedRI <- compute.AdjRI(clusters, priori.decomp)

mojosim <- compute.MoJoSim(clusters, priori.decomp)

print.mojosim <- paste(round(mojosim, 3), " (+", round((mojosim - baseline$mojosim)*100/baseline$mojosim, 1) , "\\%)", sep="")
print.f1.score <- paste(round(f1.score, 3), " (+", round((f1.score - baseline$f1.score)*100/baseline$f1.score, 1) , "\\%)", sep="")
print.adjustedRI <- paste(round(adjustedRI, 3), " (+", round((adjustedRI - baseline$adjustedRI), 3) , ")", sep="")

print(paste(print.mojosim, print.f1.score ,print.adjustedRI, sep="&"))

return(list(mojosim=mojosim, f1.score=f1.score, adjustedRI=adjustedRI))

}

load_process_All_CT_SN <- function(prname){
  
  #checking only if kernel exists, is enough
  if (file.exists(paste("benchmark", prname ,  "CT", "SN", paste(prname, "kernel.csv", sep="-"), sep="/"))) {
    kernel <- read.table(file=paste("benchmark", prname , "CT", "SN", paste(prname, "kernel.csv", sep="-"), sep="/"), 
                         sep=",", row.names = 1, header = TRUE, check.names = FALSE)
    
    return(list(kernel=kernel, dist=NULL))
    
  }
  Adj <- load_SN(prname)
  kernel <- process_All_CT_SN(Adj)
  
  write.table(kernel, file=paste("benchmark", prname , "CT", "SN", paste(prname,beta, "kernel.csv", sep="-"), sep="/"),row.names=TRUE, 
              col.names=NA,sep=",", quote=FALSE)
  
  return(list(kernel=kernel, dist=NULL))

  
}


process_All_CT_SN <- function(Adj){
  
  #   Adj <- load_SN(prname)
  
  d <- apply(abs(Adj),1,sum)
  D <- diag(d)
  
  kernel <- compute.avg.commute.time.kernel(Adj, D)
  
  kernel
}

#Semantic Network
load_process_All_SN <- function(prname, beta){
  #checking only if kernel exists, is enough
  if (file.exists(paste("benchmark", prname , "SN", paste(prname,beta, "kernel.csv", sep="-"), sep="/"))) {
    kernel <- read.table(file=paste("benchmark", prname , "SN", paste(prname,beta, "kernel.csv", sep="-"), sep="/"), 
                         sep=",", row.names = 1, header = TRUE, check.names = FALSE)
    dist <- read.table(file=paste("benchmark", prname , "SN", paste(prname,beta, "dist.csv", sep="-"), sep="/"), 
                       sep=",", row.names = 1, header = TRUE, check.names = FALSE)
    return(list(kernel=kernel, dist=dist))
    
  }
  Adj <- load_SN(prname)
  r <- process_All_SN(Adj, beta)
  
  write.table(r$kernel, file=paste("benchmark", prname , "SN", paste(prname,beta, "kernel.csv", sep="-"), sep="/"),row.names=TRUE, 
              col.names=NA,sep=",", quote=FALSE)
  write.table(r$dist, file=paste("benchmark", prname , "SN", paste(prname,beta, "dist.csv", sep="-"), sep="/"),row.names=TRUE, 
              col.names=NA,sep=",", quote=FALSE)
  
  r
}

load_process_ISA_IPO_SN <- function(prname, beta, ISA_IPO = c(0.5,0.5), rootFolder="org"){
  if ((file.exists(paste("benchmark", prname , "SN", paste(prname, beta, "ISA", "kernel.csv", sep="-"), sep="/"))) 
      && (file.exists(paste("benchmark", prname , "SN", paste(prname, beta, "IPO", "kernel.csv", sep="-"), sep="/")))){
    
    ISA_kernel <- read.table(file=paste("benchmark", prname , "SN", paste(prname,beta,"ISA", "kernel.csv", sep="-"), sep="/"), 
                             sep=",", row.names = 1, header = TRUE, check.names = FALSE)
    ISA_dist <- read.table(file=paste("benchmark", prname , "SN", paste(prname,beta,"ISA", "dist.csv", sep="-"), sep="/"), 
                           sep=",", row.names = 1, header = TRUE, check.names = FALSE)
    
    IPO_kernel <- read.table(file=paste("benchmark", prname , "SN", paste(prname,beta, "IPO", "kernel.csv", sep="-"), sep="/"), 
                             sep=",", row.names = 1, header = TRUE, check.names = FALSE)
    IPO_dist <- read.table(file=paste("benchmark", prname , "SN", paste(prname,beta, "IPO", "dist.csv", sep="-"), sep="/"), 
                           sep=",", row.names = 1, header = TRUE, check.names = FALSE) 
    
    kernel = ISA_IPO[1] * ISA_kernel + ISA_IPO[2] * IPO_kernel
    dist = ISA_IPO[1] * ISA_dist + ISA_IPO[2] * IPO_dist
    
    
    return(list(kernel=kernel, dist=dist))
    
  }
  
  Adj <- load_SN(prname)
  r1 <- process_ISA_SN(Adj, beta, rootFolder)
  r2 <- process_IPO_SN(Adj, beta, rootFolder)
  
  write.table(r1$kernel, file=paste("benchmark", prname , "SN", paste(prname, beta, "ISA", "kernel.csv", sep="-"), sep="/"),row.names=TRUE, 
              col.names=NA,sep=",", quote=FALSE)
  write.table(r1$dist, file=paste("benchmark", prname , "SN", paste(prname,beta, "ISA", "dist.csv", sep="-"), sep="/"),row.names=TRUE, 
              col.names=NA,sep=",", quote=FALSE)
  
  write.table(r2$kernel, file=paste("benchmark", prname , "SN", paste(prname, beta, "IPO", "kernel.csv", sep="-"), sep="/"),row.names=TRUE, 
              col.names=NA,sep=",", quote=FALSE)
  write.table(r2$dist, file=paste("benchmark", prname , "SN", paste(prname,beta, "IPO", "dist.csv", sep="-"), sep="/"),row.names=TRUE, 
              col.names=NA,sep=",", quote=FALSE)
  
  kernel = ISA_IPO[1] * r1$kernel + ISA_IPO[2] * r2$kernel
  dist = ISA_IPO[1] * r1$dist + ISA_IPO[2] * r2$dist
  
  return(list(kernel=kernel, dist=dist))
  
}

#TESTED
merge_names_by_lower_case <- function(mydata, dimension=1){
  
#   stopifnot(dim(adj)[1] == dim(adj)[2])
#   stopifnot(rownames(adj) == colnames(adj))
  
  lower_names <- tolower(dimnames(mydata)[[dimension]])
  
  equiv_names = list()
  
  names <- dimnames(mydata)[[dimension]]
  
  duplicates = c()
  
  
  equivalent_indices <- lapply(names, function(name) which(tolower(name)==lower_names))
  
  for (i in 1:length(equivalent_indices))
    if (length(equivalent_indices[[i]]) > 1)
      if (!(any(equivalent_indices[[i]] %in% duplicates))){
    
        if (dimension==1){
          #add by rowsums
          added <- colSums(mydata[equivalent_indices[[i]], ])
          mydata[equivalent_indices[[i]][1],] <- added
        }
        else{
          added <- rowSums(mydata[, equivalent_indices[[i]]])
          mydata[,equivalent_indices[[i]][1]] <- added
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

#TESTED
load_SN <- function(prname, make_symmetric=T, makeTopNode=T, eliminate_local_variables=T, rootNode="java.lang.Object", identifiers=c()){
  setwd("~/workspace")
  mySN = read.csv(paste("benchmark", prname , "SN", paste(prname, "SN.csv", sep="-"), sep="/"), sep = ",", quote = "\"", dec= ".")
  # mySN[which(is.na(mySN[,1])),1] <- NA
  indices <- which(is.na(mySN[,1]))
  
  if(length(indices) > 0){
    mySN <- mySN[-indices, - (indices + 1)]
  }
  
  rownames(mySN) <- mySN[,1]
  mySN <- mySN[,-1]
  gc()
  mySN <- data.matrix(mySN)
  
  gc()
  colnames(mySN) <- rownames(mySN)

  startIndex <- get.start.index.of.types(colnames(mySN))
  
  #Siz of dictionary
  D <-startIndex-1
  
  #Some memory optimization #1
#   x <- merge_names_by_lower_case(mySN[1:startIndex-1,], 1)
#   y <- mySN[startIndex:nrow(mySN),]
#   mySN <- 0
#   gc()
#   
#   mySN <- rbind(x, y)
#   x<-0
#   y <- 0
#   gc()
#   
#   
#   #Some memory optimization #2
#   x <- merge_names_by_lower_case(mySN[,1:startIndex-1], 2)
#   y <- mySN[,startIndex:ncol(mySN)]
#   mySN <- 0
#   gc()
#   mySN <- cbind(x,y)
#   x<-0
#   y <- 0
#   gc()

#Filter on the identifiers; otherwise use all the identifiers
#It is safe to intersect on uppercase (camel-case for Java)

# NOW all info
# identifiers <- rownames(mySN[1:D,])
  if (length(identifiers) > 0){
    lower_case_sensitive_identifiers <- tolower(rownames(mySN)[1:D])
    
    case_insensitive_identifiers <- intersect(lower_case_sensitive_identifiers, tolower(identifiers))
    
    identifier_indices <- which(lower_case_sensitive_identifiers %in% case_insensitive_identifiers)

    mySN <- rbind(merge_names_by_lower_case(mySN[identifier_indices,], 1), mySN[startIndex:nrow(mySN),])  
    mySN <- cbind(merge_names_by_lower_case(mySN[,identifier_indices], 2), mySN[,startIndex:ncol(mySN)])
    
  }else{
    mySN <- rbind(merge_names_by_lower_case(mySN[1:D,], 1), mySN[startIndex:nrow(mySN),])  
    mySN <- cbind(merge_names_by_lower_case(mySN[,1:D], 2), mySN[,startIndex:ncol(mySN)])
  }

  setTopNode <- function(Adj, rootNode){
    require(igraph)
    stopifnot(all(rownames(Adj)== colnames(Adj)))
    names <- colnames(Adj)
    startIndex <- get.start.index.of.types(names)

    print("First Type Node found at:")
    print(startIndex)
    
    #size of dictionary
    D <- startIndex -1
    #size of types
    V <- dim(Adj)[2] - startIndex +1

    gc()
    dims <- dim(Adj)
    
    
    S <- Adj[startIndex:dims[1], startIndex:dims[2]]
    Adj_dimnames <- dimnames(Adj)

    
    rownames(S) <- Adj_dimnames[[1]][startIndex:dims[1]]
    colnames(S) <- colnames(Adj)[startIndex:dims[2]]
    
    topNode <- which(rownames(S) == rootNode)  
    
    g <- graph.adjacency(S, mode='directed', diag=FALSE)
    
    count <-0
    
    for (i in 1:V)
      if (i != topNode){
        s <- shortest.paths(g, v=topNode, to=i, mode="in")
        
        if (s==Inf){
          Adj[colnames(s), rownames(s)] <- 1
          
          count <- count+1
        }
      }
    
    Adj
  }
  

  remove_unknown_node <- function(Adj) {
    
    allNames <- colnames(Adj)
    
    startIndex <- get.start.index.of.types(allNames)
    
    #Remove unknown type
    unknownIdx <- which(allNames[startIndex:length(allNames)] == "Unknown")
    
    if(length(unknownIdx) > 0)
      Adj <- Adj[-unknownIdx,-unknownIdx]
      
    Adj
  }
  

  remove_empty_nodes <- function(Adj){
    #Remove nodes with no edges
    empty_rows <- which(apply(Adj,1,FUN = function(x){all(x == 0)}))
    empty_cols <- which(apply(Adj,2,FUN = function(x){all(x == 0)}))
    exclude_empty_elements <- intersect(empty_rows, empty_cols)
    
    if (length(exclude_empty_elements) > 0)
      Adj <- Adj[-exclude_empty_elements, -exclude_empty_elements]
    
    Adj
  }

  mySN <- remove_unknown_node(mySN)
  mySN <- remove_empty_nodes(mySN)

  #Make a virtual top node, i.e. java.lang.Object
  if (makeTopNode)
    mySN <- setTopNode(mySN, rootNode)
  

  #Eliminate Local variables
  eliminateLocalVariables <- function(Adj){
    #find the start Index for the types in SN
    names <- colnames(Adj)
    startIndex <- get.start.index.of.types(names)

    print("starting Index")
    print(startIndex)
    
    #size of Dictionary
    D <- startIndex - 1
    
    #Remove contaninmet dependencies from method names to local variables
    Adj[1:D,1:D] <- 0

    # Get the type contents
    C <- Adj[startIndex:dim(Adj)[1], 1:D] 
    dimnames(C) <- dimnames(Adj[startIndex:dim(Adj)[1], 1:D])
    
    exposed_identifiers <- colSums(C)
    excluded_names_indices <- which(exposed_identifiers == 0)
    
    Adj <- Adj[-excluded_names_indices, -excluded_names_indices]
    
  
    Adj <- remove_empty_nodes(Adj)
    
    Adj
  }


  if (eliminate_local_variables)
    mySN <- eliminateLocalVariables(mySN)


  if (make_symmetric){
    mySN <- 0.5 * (mySN + t(mySN))
    mySN <- mySN[which(!apply(mySN,1,FUN = function(x){all(x == 0)})),which(!apply(mySN,2,FUN = function(x){all(x == 0)}))]
  }
  #   names <- rownames(Adj)
  
  # Adj <- Adj[which(unlist(lapply(names, nchar))>4),which(unlist(lapply(names, nchar))>4)]
  #TODO get rid of the identifier names not used.
  print(dim(mySN))
  
  
  mySN
}    

process_All_SN <- function(Adj, beta){
  require(igraph)
  #   Adj <- load_SN(prname)
  
  g <- build.graph(Adj)
  r <- graph.diffusion(g, beta= beta, v=V(g))
  dimnames(r$kernel) <- dimnames(Adj)
  dimnames(r$dist) <- dimnames(Adj)
  
  r
}

get.start.index.of.types <- function(mySN_names){
  
  #FIND THE STARTING INDEX OF TYPE NAMES  
  classTypeIndex <- which(unlist(gregexpr(pattern ="\\.",mySN_names)) > 0)[1]
  primitiveTypeIndices <- which(mySN_names %in% c("float", "int", "char", "byte", "void", "double", "boolean"))
  startIndex <- min(c(classTypeIndex, primitiveTypeIndices))
  startIndex
}

process_ISA_SN <- function(Adj, beta, rootFolder="org"){
  #remove Adj of types to identifier names
  #find the first occurence of a type name in the colnames(ADj), set everything else to 0
  names <- colnames(Adj)
  startIndex <- get.start.index.of.types(names)
  
  #Also remove contaninmet dependencies from method names to local variables
  for (i in 1:dim(Adj)[1])
    for (j in 1:(startIndex-1))
      Adj[i,j] <- 0
  
  
  g <- build.graph(Adj)
  r <- graph.diffusion(g, beta= beta, v=V(g))
  dimnames(r$kernel) <- dimnames(Adj)
  dimnames(r$dist) <- dimnames(Adj)
  
  r
}

process_IPO_SN <- function(Adj, beta, rootFolder="org"){
  names <- rownames(Adj)
  startIndex <- get.start.index.of.types(names)
  
  for (i in 1:(startIndex-1))
    for (j in startIndex:dim(Adj)[2])
      Adj[i,j] <- 0
  
  
  g <- build.graph(Adj)
  r <- graph.diffusion(g, beta= beta, v=V(g))
  dimnames(r$kernel) <- dimnames(Adj)
  dimnames(r$dist) <- dimnames(Adj)
  
  r
}

#Input: rows are documents, columns are identifiernames
compute_cosine_kernel <- function(Phi_d){
  cos.sim <- function(ix) 
  {
    A = Phi_d[ix[1],]
    B = Phi_d[ix[2],]
    return( sum(A*B)/sqrt(sum(A^2)*sum(B^2)) )
  }   
  n <- nrow(Phi_d) 
  cmb <- expand.grid(i=1:n, j=1:n) 
  kernel <- matrix(apply(cmb,1,cos.sim),n,n)
  
  
  rownames(kernel) <- rownames(Phi_d)
  colnames(kernel) <- rownames(kernel)
  
  kernel
}


get_top_sample_docs <- function(prname, decomp, size = 30, top = 5){
  
  cluster_sizes <- unlist(lapply(1:max(decomp), function(x) length(decomp[decomp==x])))
  
  top_cluster_sizes = sort(cluster_sizes, decreasing=T)[1:top]
  
  top_clusters <- c()
  
  for (i in 1:top)
    top_clusters <- c(top_clusters, which(cluster_sizes == top_cluster_sizes[i]))
  
  size_sample_per_cluster <- ceiling(size/top)
  
  names <- unlist(lapply(top_clusters, function(x) {
    cls = decomp[decomp==x]
    if (length(cls) > size_sample_per_cluster)
      names(cls[sample(1:length(cls), size_sample_per_cluster)])
    else
      c()
  }))
  
  names
  
}

get_sample_docs <- function(prname, decomp, size=0.25){
  require(igraph)
  require(GeLaToLab)
  setwd("~/workspace")
  
  LOWER_LIMIT = 4
  
  #number of clusters
  noc <- max(decomp)
  
  names <- unlist(lapply(1:noc, function(x) {
    cls = decomp[decomp==x]
    if (length(cls) >= LOWER_LIMIT)
#       names(cls[sample(1:length(cls), ceiling(length(cls) * size))])
      names(cls)
    else
      c()
  }))
  
  names
}