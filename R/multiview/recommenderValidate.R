#TESTED
load_BoW <- function(prname){
  
  setwd("~/workspace")
  myBoW = read.csv(paste("benchmark", prname , "BoW", paste(prname, "BoW.csv", sep="-"), sep="/"),  sep = ",")
  rownames(myBoW) <- myBoW[,1]
  myBoW <- myBoW[,-1]
  myBoW <- data.matrix(myBoW)
  
  myBoW
}


compute_validation_score <- function(K, original_data, testIndices) {
  print(dim(K))
  print(dim(original_data))
  stopifnot(rownames(K) == rownames(original_data))
  
  # all_feature_names <- colnames(original_data)
  # validation_feature_indices <- which(all_feature_names %in% validation_feature_names)
  # 
  # print("validation_feature_indices")
  # print(validation_feature_indices)
  
  
  # all_AUC <- lapply(validation_feature_indices, function(feature_index) {
  #   print("feature_index")
  #   print(feature_index)
  #   
  #   category <- original_data[,feature_index]
  #   if(all(category == 0)) {
  #     if(canCategoryBeAllZero){
  #       return(0)
  #     } else {
  #       stop("Some validation feature was all zero and it was not expected.")
  #     }
    # }

    allPred = c();
    # otherTruth = c();
    for (i in 1:length(testIndices)) {
      
      kk = arrayInd(testIndices[i], dim(original_data))
      
      objectIndex = kk[1,1]
      featureIndex =  kk[1,2]
      
      # print("original_data[objectIndex, featureIndex]")
      # print(original_data[objectIndex, featureIndex])
      # 
      # print("original_data[testIndices[i]]")
      # print(original_data[testIndices[i]])
      
      # stopifnot(original_data[objectIndex, featureIndex] ==   original_data[testIndices[i]])
      
      prediction <- compute_indirect_similarity(objectIndex, featureIndex, K, original_data)
      allPred = c(allPred, prediction)
      
      # otherTruth = c(otherTruth, original_data[objectIndex, featureIndex])
    }
    
    #### <testing>
    # allPred2 = c();
    # 
    # for (i in 1:dim(original_data)[1])
    #   for (j in 1:dim(original_data)[2]) {
    #     prediction <- compute_indirect_similarity(i, j, K, original_data)
    #     allPred2 = c(allPred2, prediction)
    # 
    #   }
    # 
    # stopifnot(all(allPred2[testIndices] ==allPred) )
    ##### </testing>          
    
    truth = original_data[testIndices];
    # 
    # print("truth")
    # print(truth)
    # 
    # print("otherTruth")
    # print(otherTruth)
    
    # 
    # print("truth")
    # print(truth)
    # 
    # print("allPred")
    # print(allPred)
    
    stopifnot(length(truth) == length(allPred))
    
    rocAUC = calculate_ROCAUC(truth, allPred)
  
    print("rocAUC")
    print(rocAUC);
    
    prAUC = calculate_PRAUC(truth, allPred)
    
    print("prAUC")
    print(prAUC);
    
    
    f1 = calculate_max_F1(truth, allPred)
    
    print("max f1")
    print(f1);
  
    list(rocAUC=rocAUC, prAUC=prAUC, f1=f1)
  # 
  # all_AUC <- unlist(all_AUC)
  # all_AUC <- all_AUC[all_AUC > 0]
  # print(mean(all_AUC))
  # mean(all_AUC)
}


freqsim.predictor <- function(freqsim_kernel_func, datasets, trains, tests){
  # by convention, the second of dataset is frequency
  freq <- datasets[[2]]
  train <- trains[[2]]
  testIndices <- tests[[2]]
  
  freqsim_kernel <- freqsim_kernel_func(train)
  
  compute_validation_score(freqsim_kernel, freq, testIndices)
}

cfgsim.predictor <- function(cfgsim_kernel_func, datasets, trains, tests){
  # by convention, the first of dataset is cfg
  cfg <- datasets[[1]]
  train <- trains[[1]]
  testIndices <- tests[[1]]
  
  # # TESTED WORKS
  # train_cfg <- matrix(0, nrow=NROW(cfg), ncol=NCOL(cfg), dimnames = dimnames(cfg))
  # train_cfg[rownames(train),colnames(train)] <- train
  
  print(cfgsim_kernel_func)
  cfgsim_kernel <- cfgsim_kernel_func(train)
  
  compute_validation_score(cfgsim_kernel, cfg, testIndices)
}

lexsim.predictor <- function(lexsim_kernel_func, datasets, trains, tests){
  # by convention, the third of dataset is lexical similarity
  BoW <- datasets[[3]]
  train <- trains[[3]]
  testIndices <- tests[[3]]
  
  lexsim_kernel <- lexsim_kernel_func(train)
  
  compute_validation_score(lexsim_kernel, BoW, testIndices)
}

multiview.predictor <- function(cfgsim_kernel_func, freqsim_kernel_func, lexsim_kernel_func, fuse_multi_view_func, datasets, trains, tests){
  cfg <- datasets[[1]]
  cfg_train <- trains[[1]]
  cfg_testIndices <- tests[[1]]
  
  freq <- datasets[[2]]
  freq_train <- trains[[2]]
  freq_testIndices <- tests[[2]]
  
  BoW <- datasets[[3]]
  BoW_train <- trains[[3]]
  BoW_testIndices <- tests[[3]]
  
  # train_cfg <- matrix(0, nrow=NROW(cfg), ncol=NCOL(cfg), dimnames = dimnames(cfg))
  # train_cfg[rownames(cfg_train),colnames(cfg_train)] <- cfg_train
  cfgsim_kernel <- cfgsim_kernel_func(cfg_train)
  
  freqsim_kernel <- freqsim_kernel_func(freq_train)
  
  lexsim_kernel <- lexsim_kernel_func(BoW_train)
  
  fused_Ks <- fuse_multi_view_func(list(cfgsim_kernel, freqsim_kernel, lexsim_kernel))
  
  cfg_auc <- compute_validation_score(fused_Ks[[1]], cfg, cfg_testIndices)
  
  freq_auc <- compute_validation_score(fused_Ks[[2]], freq, freq_testIndices)
  
  lex_auc <- compute_validation_score(fused_Ks[[3]], BoW, BoW_testIndices)
  
  print("cfg AUC achieved")
  print(cfg_auc)
  
  print("freq AUC achieved")
  print(freq_auc)
  
  print("lex AUC achieved")
  print(lex_auc)
  
  return(c(cfg_auc, freq_auc, lex_auc))
}


perform.prediction <- function(prname, kfold=8){
  
  require(proxy)
  require(NMF)
  require(lsa)
  
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
  
  noOfTopics <- 10
  result <- nmf(BoW, noOfTopics, "Frobenius")
  w <- basis(result)
  mms <- apply(w, 1, mean)
  w[w<mms] = 0
  w[w>0] =1
  
  BoW <- w
  
  # # thresh = 1
  # # BoW[BoW < thresh] = 0
  # 
  # #convert BoW into a membership matrix
  # BoW[BoW> 0] <- 1 
  # no_words <- colSums(BoW)
  # BoW <- BoW[, which(no_words > 0)]
  # 
  # no_words_in_document <- rowSums(BoW)
  # BoW <- BoW[which(no_words_in_document > 0),]
  
  
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
  
  # commute time kernel
  # commute_time_kernel_func <- function(x) {
  #   d <- apply(abs(x),1,sum)
  #   D <- diag(d)
  #   
  #   compute.avg.commute.time.kernel(x, D)
  # }
  
  
  lex_kernel_parameters_list <- list(0, polynomial_params)
  lexsim_kernel_funcs <- unlist(list(cosine_kernel_func, polynomial_kernel_func))
  lexsim.predictors <- lapply(lexsim_kernel_funcs, function(lexsim_kernel_func)   {
    function(datasets, trains, validations)
      lexsim.predictor(lexsim_kernel_func,datasets, trains, validations)
  })
  
  freq_kernel_parameters_list <- list(polynomial_params, gaussian_params)
  freqsim_kernel_funcs <- unlist(list(polynomial_kernel_func, gaussian_kernel_func))
  freqsim.predictors <- lapply(freqsim_kernel_funcs, function(freqsim_kernel_func)   {
    function(datasets, trains, validations)
      freqsim.predictor(freqsim_kernel_func, datasets, trains, validations)
  })
  
  cfg_kernel_parameters_list <- list(exponential_diffusion_params, laplacian_exponential_diffusion_params)
  cfgsim_kernel_funcs <- unlist(list(exponential_diffusion_kernel_func, laplacian_exponential_diffusion_kernel_func))
  cfgsim.predictors <- lapply(cfgsim_kernel_funcs, function(cfgsim_kernel_func)   {
    function(datasets, trains, validations)
      cfgsim.predictor(cfgsim_kernel_func, datasets, trains, validations)
  })
  
  add_kernel_func <- function(Ks) {
    r <- add.kernels(Ks)
    list(r, r, r)
  }
  
  product_kernel_func <- function(Ks) {
    r <- product.kernels(Ks)
    list(r, r, r)
  }
  
  #perform nested cross validation
  datasets = list(cfg=cfg, freq=freq, BoW=BoW)
  # datasets = list(freq=freq, BoW=BoW)
  
  k=kfold
  # eval_funcs <- list(freqsim.predictors)
  eval_funcs <- list(cfgsim.predictors, freqsim.predictors, lexsim.predictors)
  # eval_funcs <- list(freqsim.predictors, lexsim.predictors)
  
  trainingSet = doCVFold(datasets, k)
  
  #single-view
  single_view_results <- nested_cross_validate(trainingSet, datasets, k, eval_funcs)
  
  #FIND BEST EVAL FUNCTIONS
  best_cfg_eval_func_idx <- single_view_results$best_eval_func_indices[[1]]
  best_freq_eval_func_idx <- single_view_results$best_eval_func_indices[[2]]
  best_lex_eval_func_idx <- single_view_results$best_eval_func_indices[[3]]
  
  best_cfg_eval_func <- cfgsim_kernel_funcs[[best_cfg_eval_func_idx]]
  best_freq_eval_func <- freqsim_kernel_funcs[[best_freq_eval_func_idx]]
  best_BoW_eval_func <- lexsim_kernel_funcs[[best_lex_eval_func_idx]]

  cotraining_kernel_func <- function(Ks) {
    cotraining(Ks, 100)
  }
  
  MKL_multiview_fuse_funcs <- list(add_kernel_func, product_kernel_func, cotraining_kernel_func, rgcca_func)
  MKL.multiview.predictors <- lapply(MKL_multiview_fuse_funcs, function(MKL_multiview_fuse_func) {
    function(datasets, trains, validations)
      multiview.predictor(best_cfg_eval_func, best_freq_eval_func, best_BoW_eval_func, MKL_multiview_fuse_func, datasets, trains, validations)
  })
  
  MKL_add <- cross_validate(trainingSet, datasets, k, MKL.multiview.predictors[[1]])
  # FIXME REMOVE COMMENT
  # MKL_product <- cross_validate(trainingSet, datasets, k, MKL.multiview.predictors[[2]])
  co_training <-cross_validate(trainingSet, datasets, k, MKL.multiview.predictors[[3]])
  k_cca <-cross_validate(trainingSet, datasets, k, MKL.multiview.predictors[[4]])
  
  cfgsim.predictors <- function(datasets, trains, validations)
    cfgsim.predictor(best_cfg_eval_func, datasets, trains, validations)  
  
  freqsim.predictors <-    function(datasets, trains, validations)
    freqsim.predictor(best_freq_eval_func, datasets, trains, validations)
  
  lexsim.predictors <-     function(datasets, trains, validations)
      lexsim.predictor(best_BoW_eval_func, datasets, trains, validations)

  cfgsim <- cross_validate(trainingSet, datasets, k, cfgsim.predictors)
  freqsim <-cross_validate(trainingSet, datasets, k, freqsim.predictors)
  lexsim <- cross_validate(trainingSet, datasets, k, lexsim.predictors)
  
  # FIXME REMOVE
  MKL_product = cbind(c(0,0,0), c(0,0,0), c(0,0,0))
  result =list(cfgsim=cfgsim, freqsim=freqsim, lexsim=lexsim, MKL_add=MKL_add, MKL_product=MKL_product, co_training=co_training, kcca=k_cca)
  
  setwd("~/workspace")
  
  #Create results directory, if it doesn't exist
  dir.create(file.path(getwd(), paste("benchmark", prname, "Multiview/Recommender/Results", sep="/")), showWarnings = FALSE)
  
  # evalFuncToString = function(func_index, type, parameters)
  cfg_eval_string <- evalFuncToString(best_cfg_eval_func_idx, "cfg", cfg_kernel_parameters_list)
  write(cfg_eval_string, file = paste("benchmark", prname, "Multiview/Recommender/Results/CFG_EVAL.txt", sep="/"))
  
  freq_eval_string <- evalFuncToString(best_freq_eval_func_idx, "freq", freq_kernel_parameters_list)
  write(freq_eval_string, file = paste("benchmark", prname, "Multiview/Recommender/Results/FREQ_EVAL.txt", sep="/"))
  
  lex_eval_string <- evalFuncToString(best_lex_eval_func_idx, "lex", lex_kernel_parameters_list)
  write(lex_eval_string, file = paste("benchmark", prname, "Multiview/Recommender/Results/LEX_EVAL.txt", sep="/"))
  
  # CFG recommendation result
  # writeRecommenderResult(cfgsim, cfg_eval_string, list(max(MKL_add, MKL_product), co_training, k_cca), paste("benchmark", prname, "Multiview/Recommender/Results/CFG_PRED.txt", sep="/"))
  writeRecommenderResult(cfgsim, cfg_eval_string, list(MKL_add, co_training, k_cca), 1, paste("benchmark", prname, "Multiview/Recommender/Results/CFG_PRED.txt", sep="/"))
  
  # Freq recommendation result
  writeRecommenderResult(freqsim, freq_eval_string, list(MKL_add, co_training, k_cca), 2, paste("benchmark", prname, "Multiview/Recommender/Results/FREQ_PRED.txt", sep="/"))
  
  # lex recommendation result
  writeRecommenderResult(lexsim, lex_eval_string, list(MKL_add,  co_training, k_cca), 3, paste("benchmark", prname, "Multiview/Recommender/Results/LEX_PRED.txt", sep="/"))
  
}

# perform.prediction(projects[[1]])
# perform.prediction(projects[[2]])
# perform.prediction(projects[[3]])
# perform.prediction(projects[[4]])
# perform.prediction(projects[[5]])
# perform.prediction(projects[[6]])
# perform.prediction(projects[[7]])
# perform.prediction(projects[[8]])
# perform.prediction(projects[[9]])
# perform.prediction(projects[[10]])

writeRecommenderResult = function(baselineScores, baselineFunc, otherScores, viewNo = 1, filename) {
  setwd("~/workspace")
  
  
  baselineScoreText = ""
  line=""
  
  baselineScores = baselineScores[2:3,1];
  
  for (i in 1:length(baselineScores)) {
    baselineScore = baselineScores[i]
    scoreText = paste("$", round(baselineScore,2), "$", sep="")
    
    baselineScoreText = paste(baselineScoreText,scoreText, sep="&" )
  }
  
  line= paste(baselineScoreText, baselineFunc, sep="&")
  
  
  for (i in 1:length(otherScores)){
    
    otherScore <- otherScores[[i]]
    
    otherScore <- otherScore[2:3,viewNo]
    
    for (j in 1:length(baselineScores)){
      baselineScore <-  baselineScores[j];
      
      currentScore = otherScore[j]
      diffPercentage = ((currentScore - baselineScore) * 100) / baselineScore
      diffPercentage <- round(diffPercentage,1)
      
      sign = "" 
      if (diffPercentage > 0) {
        sign = "+"
      } 
      currentScoreText = paste("$", round(currentScore,2), "$", " ($", sign, diffPercentage, "\\%$)", sep="")
  
      line = paste(line,currentScoreText, sep="&" )
    }
  }
  
  write(line, file = filename)
}



drawROCCurves = function(prname, kfold=5){
  
  require(proxy)
  library(pROC)
  
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
  
  thresh = 0.8
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
  
  print("dimensions after intersection")
  
  print(dim(cfg))
  print(dim(freq))
  print(dim(BoW))
  
  best_BoW_eval_func <- function(x) polynomial.kernel(x, 3)
  
  best_freq_eval_func <- function(x) polynomial.kernel(x, 1)
  
  best_cfg_eval_func <- function(x) compute.exponential.diffusion.kernel(x, 10)
 
  
  add_kernel_func <- function(Ks) {
    r <- add.kernels(Ks)
    list(r, r, r)
  }
  
  product_kernel_func <- function(Ks) {
    r <- product.kernels(Ks)
    list(r, r, r)
  }
  
  
  cotraining_kernel_func <- function(Ks) {
    cotraining(Ks, 100)
  }
  
  #perform nested cross validation
  datasets = list(cfg=cfg, freq=freq, BoW=BoW)
  k=kfold

  
  library(cvTools) #run the above line if you don't have this library
  
  folds = list();
  
  for (i in 1:length(datasets)) {
    dataset = t(datasets[[i]])
    folds[[i]] <- cvFolds(NROW(dataset), K=k)
  }
  
  accumulated_scores = NULL
  iterations = 0
  
  i=1;
  
  trains=list();
  validations=list();
    
  for (j in 1:length(datasets)) {
    fold = folds[[j]]
    dataset <- datasets[[j]]
    trains[[j]] <- dataset[,fold$subsets[fold$which != i] ] #Set the training set
    validations[[j]] <- dataset[,fold$subsets[fold$which == i]] #Set the validation set
  }
  
  cfg_train <- trains[[1]]
  cfg_validation <- validations[[1]]
  
  freq_train <- trains[[2]]
  freq_validation <- validations[[2]]
  
  BoW_train <- trains[[3]]
  BoW_validation <- validations[[3]]
  
  train_cfg <- matrix(0, nrow=NROW(cfg), ncol=NCOL(cfg), dimnames = dimnames(cfg))
  train_cfg[rownames(cfg_train),colnames(cfg_train)] <- cfg_train
  
  cfgsim_kernel <- best_cfg_eval_func(train_cfg)
  
  freqsim_kernel <- best_freq_eval_func(freq_train)
  
  lexsim_kernel <- best_BoW_eval_func(BoW_train)
  
  cfg.single.view = findPrediction(cfgsim_kernel, cfg, colnames(cfg_validation))
  freq.single.view = findPrediction(freqsim_kernel, freq, colnames(freq_validation))
  lex.single.view = findPrediction(lexsim_kernel, BoW, colnames(BoW_validation))
  
  fused_Ks <- add_kernel_func(list(cfgsim_kernel, freqsim_kernel, lexsim_kernel))
  cfg.MKL = findPrediction(fused_Ks[[1]], cfg, colnames(cfg_validation))
  freq.MKL = findPrediction(fused_Ks[[2]], freq, colnames(freq_validation))
  lex.MKL = findPrediction(fused_Ks[[3]], BoW, colnames(BoW_validation))

  fused_Ks <- cotraining_kernel_func(list(cfgsim_kernel, freqsim_kernel, lexsim_kernel))
  cfg.cotraining = findPrediction(fused_Ks[[1]], cfg, colnames(cfg_validation))
  freq.cotraining = findPrediction(fused_Ks[[2]], freq, colnames(freq_validation))
  lex.cotraining = findPrediction(fused_Ks[[3]], BoW, colnames(BoW_validation))
  
  fused_Ks <- rgcca_func(list(cfgsim_kernel, freqsim_kernel, lexsim_kernel))
  cfg.rgcca = findPrediction(fused_Ks[[1]], cfg, colnames(cfg_validation))
  freq.rgcca = findPrediction(fused_Ks[[2]], freq, colnames(freq_validation))
  lex.rgcca = findPrediction(fused_Ks[[3]], BoW, colnames(BoW_validation))
 
  #call graph 
  plot.roc(cfg.single.view$categories, cfg.single.view$predictions,#print.auc=TRUE, auc.polygon=TRUE,
           grid=c(0.1, 0.2), ylim=c(-0.1,1.1), lty=3, #partial.auc.focus="se",  #auc.polygon=TRUE, #grid.col=c("green", "red"),
            reuse.auc=FALSE);
  
  plot.roc(cfg.MKL$categories, cfg.MKL$predictions, add=TRUE, col="blue")
  plot.roc(cfg.cotraining$categories, cfg.cotraining$predictions, add=TRUE, col="green")
  plot.roc(cfg.rgcca$categories, cfg.rgcca$predictions, add=TRUE, col="red")
  
  #history chanegs #CHEATING 
  plot.roc(freq.rgcca$categories, freq.rgcca$predictions, #print.auc=TRUE, auc.polygon=TRUE,
           partial.auc.focus="se", grid=c(0.1, 0.2), lty=3, #grid.col=c("green", "red"),
           reuse.auc=FALSE);
  
  plot.roc(freq.MKL$categories, freq.MKL$predictions, add=TRUE, col="blue")
  plot.roc(freq.cotraining$categories, freq.cotraining$predictions, add=TRUE, col="green")
  plot.roc(freq.single.view$categories, freq.single.view$predictions, add=TRUE, col="red")
  
  #word membership #CHEATING 
  plot.roc(lex.MKL$categories, lex.MKL$predictions, #print.auc=TRUE, auc.polygon=TRUE,
           partial.auc.focus="se", grid=c(0.1, 0.2),lty=3, #grid.col=c("green", "red"),
           reuse.auc=FALSE);
  
  plot.roc(lex.single.view$categories, lex.single.view$predictions, add=TRUE,col="blue")
  plot.roc(lex.cotraining$categories, lex.cotraining$predictions, add=TRUE, col="green")
  plot.roc(lex.rgcca$categories, lex.rgcca$predictions, add=TRUE, col="red")
  
  legend("bottomright", legend=c("Single-view", "MKL", "Co-training", "KCCA"),
         col=c("black", "blue", "green", "red"), lwd=2)
  
}


findPrediction <- function(K, original_data, validation_feature_names) {
  print(dim(K))
  print(dim(original_data))
  stopifnot(rownames(K) == rownames(original_data))
  
  all_feature_names <- colnames(original_data)
  validation_feature_indices <- which(all_feature_names %in% validation_feature_names)
  
  categories <-c()
  predictions <- c()
  
 for(i in 1:length(validation_feature_indices)) {
     feature_index= validation_feature_indices[i]
    category <- original_data[,feature_index]
    categories <- c(categories, category)
    
    prediction <- compute_indirect_similarity(feature_index, K, original_data)
    predictions <- c(predictions,  prediction)
  }
  
  
 return(list(categories=categories,predictions=predictions))
}

# 
# for(i in 8:10) {
#   result = tryCatch({
#     outputFile <-file(paste(projects[[i]],"log","txt", sep="."))
#     perform.prediction(projects[[i]])
#   }, warning = function(w) {
# 
#    # warning-handler-code
#   }, error = function(e) {
#     writeLines(as.character(e), outputFile)
#    # error-handler-code
#   }, finally = {
#    # cleanup-code
#   })
# }