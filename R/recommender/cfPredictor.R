cfPredict <- function(prname, classnames, single_view_results, noOfTopics=7){
  
  require(proxy)
  require(NMF)
  
  setwd("~/workspace")
  
  #Read the set of classnames for running the experiment
  #classnames <- unlist(read.table(paste("benchmark", prname , paste("MULTIVIEW", "classnames.txt" ,sep="/") , sep="/")) )
  
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
  
  freq[is.na(freq)] <- 0
  freq[freq> 0] <- 1 
  
  no_transactions <- colSums(freq)
  freq <- freq[, which(no_transactions > 0)]
  
  #Load the bag of words
  BoW <- load_BoW(prname)
  # apply tf-idf mechanism and eliminate features that are lower than some threshold,
  # then remove those features from BoW
  x <- apply_tf_idf(BoW)
  dimnames(x) <- dimnames(BoW)
  BoW <- x

  result <- nmf(BoW, noOfTopics, "Frobenius")
  w <- basis(result)
  mms <- apply(w, 1, mean)
  w[w<mms] = 0
  w[w>0] =1
  
  BoW <- w
  
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
  
  print("dimensions after intersection")
  
  print(dim(cfg))
  print(dim(freq))
  print(dim(BoW))
  #lexsim kernels
  #cosine normalized linear kernel
  cosine_kernel_func <- function(x) { cos.sim(t(x))}
  
  # polnomial degree: 1,2,3,4,5
  polynomial_params <- c(1,2,3,4,5)
  polynomial_kernel_func <- lapply(polynomial_params, function(p) {function(x) polynomial.kernel(x, p)})
  
  #gaussian parameter: 10^(-5),10^(-4),10^(-3),10^(-2),10^(-1),10^(0),10^(1),10^(2)
  gaussian_params <- c(10^(-5),10^(-4),10^(-3),10^(-2),10^(-1),10^(0),10^(1),10^(2))
  gaussian_kernel_func <- lapply(gaussian_params, function(p) {function(x) gaussian.kernel(x, p)})
  
  #exponential diffusion parameter: 10^(-5),10^(-4),10^(-3),10^(-2),10^(-1),10^(0),10^(1),10^(2)
  exponential_diffusion_params <- c(10^(-5),10^(-4),10^(-3),10^(-2),10^(-1),10^(0),10^(1),10^(2))
  exponential_diffusion_kernel_func <- lapply(exponential_diffusion_params, function(p) {function(x) calc.diffusion.kernel(x, p, TRUE)})
  
  #laplacian exponential diffusion parameter: 10^(-5),10^(-4),10^(-3),10^(-2),10^(-1),10^(0),10^(1),10^(2)
  # what is this about?!?
  laplacian_exponential_diffusion_params <- c(1,2,5,10,15,20)
  laplacian_exponential_diffusion_kernel_func <- lapply(laplacian_exponential_diffusion_params, 
                                                        function(p) {function(x) compute.exponential.diffusion.kernel(x, p)})
  
  
  lex_kernel_parameters_list <- list(0, polynomial_params)
  lexsim_kernel_funcs <- unlist(list(cosine_kernel_func, polynomial_kernel_func))
  
  freq_kernel_parameters_list <- list(polynomial_params, gaussian_params)
  freqsim_kernel_funcs <- unlist(list(polynomial_kernel_func, gaussian_kernel_func))
  
  cfg_kernel_parameters_list <- list(exponential_diffusion_params, laplacian_exponential_diffusion_params)
  cfgsim_kernel_funcs <- unlist(list(exponential_diffusion_kernel_func, laplacian_exponential_diffusion_kernel_func))
  
  
  add_kernel_func <- function(Ks) {
    r <- add.kernels(Ks)
    list(r, r, r)
  }
  
  
  best_cfg_eval_func_idx <- single_view_results$best_eval_func_indices[[1]]
  best_freq_eval_func_idx <- single_view_results$best_eval_func_indices[[2]]
  best_lex_eval_func_idx <- single_view_results$best_eval_func_indices[[3]]
  
  best_cfg_eval_func <- cfgsim_kernel_funcs[[best_cfg_eval_func_idx]]
  best_freq_eval_func <- freqsim_kernel_funcs[[best_freq_eval_func_idx]]
  best_BoW_eval_func <- lexsim_kernel_funcs[[best_lex_eval_func_idx]]
  
  
  cfgsim_kernel <- best_cfg_eval_func(cfg)
  freqsim_kernel <- best_freq_eval_func(freq)
  lexsim_kernel <- best_BoW_eval_func(BoW)
  
  # Kernel addition is used for data fusion
  fused_Ks <- add_kernel_func(list(cfgsim_kernel, freqsim_kernel, lexsim_kernel))
  
  for (c in 1:length(classnames)){
    classname = classnames[c]
    
    #CFG
    classIndex = which(rownames(cfg) == classname)
    featureIndices = which(cfg[classname,] == 0)
    cfgPrediction <- predictMembership(fused_Ks[[1]], cfg, classIndex, featureIndices)
    sortedCFGPrediction <- sort(cfgPrediction, decreasing = T, index.return=TRUE)
    
    print("sortedCFGPrediction")
    print(sortedCFGPrediction$x[1:3])
    
    ns <- colnames(cfg)[featureIndices[sortedCFGPrediction$ix[1:3]]]
    
    print(paste(classname, " should have a calling relationships with ", ns, sep=" "))
  }
  
  
  #Freq
  for (c in 1:length(classnames)){  
    classname = classnames[c]
    
    newFreq <- cbind(freq, rep(0, dim(freq)[1] ))
    
    classIndex = which(rownames(newFreq) == classname)
    newFreq[classIndex,dim(newFreq)[2]] <- 1
    
    allPred <- c()
    otherClasses = c(1:(classIndex-1), (classIndex+1):dim(newFreq)[1])
    
    for (i in otherClasses){
      freqPrediction <- predictMembership(fused_Ks[[2]], newFreq, i, dim(newFreq)[2])
      allPred <- c(allPred, freqPrediction)
    }
    sortedFreqPrediction <- sort(allPred, decreasing = T, index.return=TRUE)
    
    print("sortedFreqPrediction")
    print(sortedFreqPrediction$x[1:3])
    
    ns <- rownames(newFreq)[otherClasses[sortedFreqPrediction$ix[1:2]]]
    print(paste(classname, "should be committed with ", ns, sep=" "))
  }
  
  
  for (c in 1:length(classnames)){
    classname = classnames[c]
    
    #Topic
    classIndex = which(rownames(BoW) == classname)
    featureIndices = which(BoW[classname,] == 0)
    bowPrediction <- predictMembership(fused_Ks[[3]], BoW, classIndex, featureIndices)
    sortedBoWPrediction <- sort(bowPrediction, decreasing = T, index.return=TRUE)
    featureIndices[sortedBoWPrediction$ix[1]]
    
    print("sortedBoWPrediction")
    print(sortedBoWPrediction$x[1:3])
    
    ns <- featureIndices[sortedBoWPrediction$ix[1:3]]
    
    print(paste(classname, " should cover concepts ", ns, sep=" "))
    
  }
}

predictMembership <- function(K, data, classIndex, featureIndices) {
  allPred = c();
  for (i in 1:length(featureIndices)) {
    
    featureIndex = featureIndices[i];
    
    prediction <- compute_indirect_similarity(classIndex, featureIndex, K, data)
    allPred = c(allPred, prediction)
  }
  
  allPred
} 