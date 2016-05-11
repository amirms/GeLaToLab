add.kernels = function (Ks) {
  #check if the sizes match, more than one
  
  added_K <- Ks[[1]]
  for (i in 2:length(Ks))
    added_K <- added_K + Ks[[i]]
  dimnames(added_K) <- dimnames(Ks[[1]])
  added_K
}


product.kernels = function(Ks) {
  
  added_K <- Ks[[1]]
  for (i in 2:length(Ks))
    added_K <- added_K * Ks[[i]]
  dimnames(added_K) <- dimnames(Ks[[1]])
  added_K
}

row.normalize <- function(df) {
  df / rowSums(sqrt(df^2))
  
}
  
make.mean.zero <- function(df) {
  rowMeans <- rowSums(df) / dim(df)[2]
  
  for (i in 1:dim(df)[1])
    
    df[i,] <- df[i,] - rowMeans[i]
  
  df
}  
  
#Input: X = [p x n]  
estimate.error <- function(X) {
  
  p <- dim(X)[1]
  N <- dim(X)[2]
  
  phi <- rep(0, p)
  
  for ( j in 1:p) {
    sum = 0
    for (r in 1:N)
      sum<- sum + X[j,r]^2
    
    
    print(sum)
    phi[j] <- 1 - sum
    
  }
  diag(phi)
  
}  

#Input:
#prname = project name
#k= number of clusters
compute.multiview <- function(prname, choice, k=0, rootFolder="org", pattern = "*.java") {
  require(proxy)
  
  setwd("~/workspace")
  
  #Read the set of classnames for running the experiment
  classnames <- unlist(read.table(paste("benchmark", prname , paste("MULTIVIEW", "classnames.txt" ,sep="/") , sep="/")) )
  
  
  #Load the priori decomposition
  decomposition <- read.csv(paste("benchmark", prname ,"decomposition.csv", sep="/"), sep=",",  header = TRUE)
  priori.decomp <- decomposition$x
  names(priori.decomp) <- decomposition$X
  priori.decomp <- priori.decomp[which(names(priori.decomp) %in% classnames)] # find the ones in classnames
  priori.decomp <- priori.decomp[order(names(priori.decomp))]
  priori.decomp <- normalizeVector(priori.decomp)

  #Load the adjacency matrix
  extensions= c("java/", "org/xml/", "javax/")
  cfg <- import.bunch.matrix(paste("benchmark", prname ,"dep_graph.txt", sep="/"), exclude.ext=extensions)
  #cfg <- read.table("benchmark/jedit-5.1.0/cfg.csv", sep=",", row.names = 1, header = TRUE, check.names = FALSE)
  #   cfg <- unweight.adjacency(cfg)
  cfg <- cfg[which(rownames(cfg) %in% classnames), which(colnames(cfg) %in% classnames)]
  cfg <- cfg[order(rownames(cfg)), order(colnames(cfg))]

  #Load the bag of words
#   if (choice)
#   bow <- read.table(paste("benchmark", prname , "mydata-BoW-matrix.csv", sep="/"), sep=",", row.names = 1, header = TRUE, check.names = FALSE)

  #Load the transaction frequency
  freq <- read.table(paste("benchmark", prname , "mydata-change-freq-matrix.csv", sep="/"), sep=",", row.names = 1, header = TRUE, check.names = FALSE)  
  freq <- as.matrix(freq)
  freq <- freq[which(rownames(freq) %in% classnames),]
  freq <- freq[order(rownames(freq)),]

  #Initialize the text of source files
  setwd(paste("benchmark", prname, sep="/"))
  txts <- read.text.directory(rootFolder, pattern)
  txts <- txts[which(names(txts) %in% classnames)] #Need to check if this work FIXME
  txts <- txts[order(names(txts))]

  lexsim.kernel <- spectrum.string.kernel(txts, 6)

  # no self-references
  diag(cfg) <- 0 
  cfg <- make.symmetric(cfg)
  d <- apply(abs(cfg),1,sum)
    
  indices <- which(d==0)
    
  while(length(indices) > 0) {
    cfg <- cfg[-indices, -indices]
    d <- apply(abs(cfg),1,sum)
    indices <- which(d==0)
  }
    
  d <- apply(abs(cfg),1,sum)
  D <- diag(d)

  cfg.laplacian = laplacian(cfg, TRUE)
  
  cfg.kernel <- calc.diffusion.kernel(cfg.laplacian, beta = 5)
  

  N = dim(cfg)[1]
  
  print(length(txts))
  print(dim(cfg.kernel))


  #TEST see if the intersection of all actually works
  for (i in 1:length(names(txts))) {
    print(names(txts)[i] )
    print(rownames(cfg)[i])
    
    if (!(names(txts)[i] == rownames(cfg)[i]))
      stop('names do not match!')
    
  }
  
  #bow: a document-term matrix: documents in columns and terms in rows
#   bow <- t(bow)
#   
#   bow <- lw_logtf(bow) * gw_idf(bow)
#   
#   bow <- bow[complete.cases(bow),]
#   
#   ndims = compute.rank(dim(bow)[2], dim(bow)[1])
#   print("No of Dimensions:")
#   print(ndims)
#   
#   space1 = lsa(bow, ndims)
#   
#   lsa.bow = diag(sqrt(space1$sk)) %*% t(space1$dk)
#   
#   #make lsa bow, to documents in rows and terms in columns
#   lsa.bow <-  t(lsa.bow)
  
  #Check to ensure  lsa.bow is a documentTerm matrix
#   lexsim.kernel = gaussian.kernel(lsa.bow, 0.01) #ndims/N
  

  #Process the transaction frequency
  no_transactions <- colSums(freq)
  
  freq <- freq[, which(no_transactions <= 30)]

  freq.kernel = gaussian.kernel(freq, 0.01)


  print(dim(freq.kernel))
  print(dim(cfg.kernel))
  print(dim(lexsim.kernel))


#   dummy_v <- rep(0, dim(cfg)[1])
#   names(dummy_v) <- rownames(cfg)
#   
#   #Fix priori decomposition 
#   priori.decomp <- find.intersection(priori.decomp, dummy_v)

  #Normalize the priori decomposition
  priori.decomp <- normalizeVector(priori.decomp)
  priori.decomp <- priori.decomp[order(names(priori.decomp))]


  print("printing the numer of groups:")
  print(max(priori.decomp))

  k <- max(priori.decomp)

# return(cfg.kernel)

# return(list(cfg.kernel = cfg.kernel, lexsim.kernel = lexsim.kernel))

#   print(dim(cfg.kernel))
#   print(dim(lexsim.kernel))
  # stop("naoad")

  kernels = list()
  if (choice == "all")
  {
    kernels[[1]] <- cfg.kernel
    kernels[[2]] <- lexsim.kernel
    kernels[[3]] <- freq.kernel
    kernels[[4]] <- add.kernels(list(cfg.kernel, lexsim.kernel, freq.kernel))
    kernels[[5]] <- product.kernels(list(cfg.kernel, lexsim.kernel, freq.kernel))
    kernels[[6]] <- compute.cotraining(list(cfg.kernel, lexsim.kernel, freq.kernel), k, iter=50, priori.decomp, prname)
    kernels[[7]] <- compute_generalized_pareto_multiview(list(laplacian(cfg.kernel, TRUE), laplacian(lexsim.kernel, TRUE), laplacian(freq.kernel, TRUE)))    
  }
  else
    kernels <- list(switch(choice,
           MQ = cfg.kernel,
           CQ = lexsim.kernel,
           FQ = freq.kernel,
           MKL.Add = add.kernels(list(cfg.kernel, lexsim.kernel, freq.kernel)),
           MKL.Prod = product.kernels(list(cfg.kernel, lexsim.kernel, freq.kernel)),
           CT = compute.cotraining(list(cfg.kernel, lexsim.kernel, freq.kernel), k, iter=50, priori.decomp, prname),
           MO = compute_generalized_pareto_multiview(list(laplacian(cfg.kernel, TRUE), laplacian(lexsim.kernel, TRUE), laplacian(freq.kernel, TRUE))))
            )

  
  lapply(kernels, function(kernel) compute_cluster_kernel(kernel, k, priori.decomp) )
  
}

#Input: kernel matrix, no of clusters, and priori.decomp
compute_cluster_kernel <- function(kernel, k, priori.decomp){
  
  result <- kmeans(kernel, centers = k, iter.max = 500, nstart = 20000)$cluster
  # result <- spectral.clustering(addition.kernel, k, iter = 15)
  # result <- spectral.clustering(laplacian(addition.kernel, TRUE), k, iter = 40)
  
  result <- normalizeVector(result)
  
  #THis is necessary for kernel produced from co-training
  names(result) <- names(priori.decomp)
  
  precision <- compute.precision(result, priori.decomp)
  recall <- compute.recall(result, priori.decomp)
  f1.score <- compute.f1(result, priori.decomp)
  
  nmi <- compute.NMI(result, priori.decomp)
  
  mojosim <- compute.MoJoSim(result, priori.decomp)
  
  ri <- compute.RI(result, priori.decomp)
  
  adjri <- compute.AdjRI(result, priori.decomp)
  
  return(list(mojosim=mojosim, precision=precision, recall=recall, f1.score=f1.score, nmi=nmi, ri=ri, adjri=adjri))
}

perform.all <- function(mode="all", indices=0){
#     projects <- list("apache-ant-1.9.3", "hadoop-0.20.2", "apache-log4j-1.2.17", "eclipse-jdt-core-3.8")#,
  projects <- list("jdom-2.0.5", "jedit-5.1.0", "jfreechart-1.2.0", "jhotdraw-7.0.6", "junit-4.12")# ,"weka-3.6.11")
#       
#     names(projects) <- c("Apache Ant", "Apache Hadoop", "Apache Log4j", "Eclipse JDT Core")#, 
  names(projects) <- c("JDOM", "JEdit", "JFreeChart", "JHotDraw", "JUnit")#, "Weka")
# 
    projects <- list("weka-3.6.11")
    names(projects) <- c( "Weka")
  
#     projects <- list("apache-log4j-1.2.17")#, "apache-log4j-1.2.17", "eclipse-jdt-core-3.8")#,
#     names(projects) <- c("Apache Log4j") #,  "Apache Log4j", "Eclipse JDT Core")#, 
#     projects <- list("junit-4.12")
#     names(projects) <- c("JUnit")


  #Project weka has a different root folder than others!

    results <- lapply(projects, function(p) {
                                    rootFolder <- "org"
                                    if (p=="weka-3.6.11")
                                      rootFolder <- "weka"
                                    compute.multiview(p, mode, rootFolder)
                        })
    
    names(results) <- names(projects)

# return(results)
    
    # plot(c(min(decay),max(decay)),c(0,1),type='n',xlab="lambda",ylab='MoJoSim')

    
    #\cellcolor[gray]{0.8}
    
    if (indices == 0)
      indices <- 1:length(projects)
    
    str2tex = ""


  for (i in 1:length(results)) {
    
    line <- paste(names(results[i]), round(results[[i]][[1]]$mojosim, 3), round(results[[i]][[2]]$mojosim, 3), 
                  round(results[[i]][[3]]$mojosim, 3), round(results[[i]][[4]]$mojosim, 3), round(results[[i]][[5]]$mojosim, 3),
                  round(results[[i]][[6]]$mojosim, 3), round(results[[i]][[7]]$mojosim, 3), "\\", sep="&" )
    str2tex <- paste(str2tex, line , "\n" )
    
  }

    
#     for (i in 1:length(results)) {
#       
#       line <- paste(names(projects[i]), round(results[[i]]$precision, 3), round(results[[i]]$recall, 3), round(results[[i]]$f1.score, 3),
#                      round(results[[i]]$mojosim, 3), round(results[[i]]$nmi, 3), round(results[[i]]$adjri, 3), "\\\\", sep="&" )
#       str2tex <- paste(str2tex, line , "\n" )
#       
#     }
    
    str2tex
  
}



