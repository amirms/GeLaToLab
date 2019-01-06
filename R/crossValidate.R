#Input
# data
# k: the number of folds
# pred_func
# eval_func is a list of list of evaluation functions


doCVFold = function(datasets, k){
  library(cvTools) #run the above line if you don't have this library
  
  folds = list();
  
  #check it is Document x features
  
  for (i in 1:length(datasets)) {
    dataset = datasets[[i]];
    
    x <- dataset
    
    #hacky. assumes binary dataset
    x[which(rowSums(x) <= 1),] = 0;
    
    linkIndices <- which(x > 0)
    linkFolds = cvFolds(length(linkIndices), K=k)
    
    maxIndex <- dim(dataset)[1] * dim(dataset)[2]
    
    otherIndices <- which(dataset == 0)
    otherFolds <- cvFolds(length(otherIndices), K=k)
    
    folds[[i]] <- list(linkFolds=linkFolds, linkIndices=linkIndices, otherIndices=otherIndices, otherFolds=otherFolds)
  }
  
  # perform nested cross validation
  #exclude set number 1
  
  cvSetsPerDataset = list()
  
  for (i in 1:length(datasets)) {
    dataset = datasets[[i]]
    fold = folds[[i]]
    
    cvSets=list()
    
    for (j in 1:k){
      linkFolds = fold$linkFolds;
      otherFolds = fold$otherFolds;
      
      linkIndices = fold$linkIndices;
      otherIndices = fold$otherIndices;
      
      excludedLinkIndices = linkIndices[linkFolds$subsets[linkFolds$which == j]]
      excludedOtherIndices = otherIndices[otherFolds$subsets[otherFolds$which != 1]]
      
      cvSet =  dataset;
      cvSet[excludedLinkIndices] = 0;
        
      cvSets[[j]] <- list(cvSet=cvSet, indices = union(excludedLinkIndices, excludedOtherIndices))
    }
    
    cvSetsPerDataset[[i]] <- cvSets
  }
  
 return(cvSetsPerDataset)
}

nested_cross_validate = function(trainingSet, datasets, k, eval_funcs) {
  # #install_packages("cvTools")
  # library(cvTools) #run the above line if you don't have this library
  # 
  # folds = list();
  # 
  # for (i in 1:length(datasets)) {
  #   dataset = t(datasets[[i]])
  #   console.log("nrows of dataset")
  #   console.log(NROW(dataset))
  #   folds[[i]] <- cvFolds(NROW(dataset), K=k)
  # }
  # 
  # # perform nested cross validation
  # #exclude set number 1
  # validationsets=list()
  # for (i in 1:length(datasets)) {
  #   dataset = datasets[[i]]
  #   fold = folds[[i]]
  #   validationset =  dataset[,fold$subsets[fold$which != 1]]
  #   validationsets[[i]] <- validationset
  # }
  # 
  # vfolds = list();
  # 
  # for (i in 1:length(datasets)) {
  #   validationset = t(validationsets[[i]])
  #   vfolds[[i]] <- cvFolds(NROW(validationset), K=k)
  # }

  validationSet = doCVFold(lapply(trainingSet, function(s) s[[1]]$cvSet), k)
  
  best_eval_func_indices <- list()
  best_scores <- list()
  
  for (i in 1:length(eval_funcs)) {
    eval_func = eval_funcs[[i]]
    
    best_score <- 0
    best_eval_func_index = -1
    
    for (j in 1:length(eval_func)) {
      func <- eval_func[[j]]
      
      accumulated_score =0
      iterations = 0
      
      for (m in 1:k) {
        trains = lapply(validationSet, function(s) s[[m]]$cvSet)
        validationIndices = lapply(validationSet, function(s) s[[m]]$indices)
        
        # Use ROC AUC to find the best eval function
        scores <- func(datasets, trains, validationIndices)
        print("scores")
        print(scores)
        
        score <- scores$rocAUC
        
        if (score == 0) {
          stop("score returned 0");
        }
        
        iterations <- iterations + 1
        accumulated_score <- accumulated_score + score
      }
      
      avg_score <- accumulated_score / iterations 
      
      if (is.finite(avg_score) && avg_score > best_score ){
        best_score <- avg_score
        best_eval_func_index <- j
      }
    }
    
    best_eval_func_indices[[i]] <- best_eval_func_index
    best_scores[[i]] <- best_score
  }
  
  print(best_eval_func_indices)
  print(best_scores)
  
  test_scores = list()
  
  for (m in 1:length(best_eval_func_indices)) {
    best_eval_func_index <- best_eval_func_indices[[m]]
    func = eval_funcs[[m]][[best_eval_func_index]]
  
    accumulated_score <- 0
    iterations = 0
    
    for(i in 1:k){
      # score <- tryCatch(func(datasets, trains, validations), error = function(e) { print(e); return(0)})
      trains = lapply(trainingSet, function(s) s[[i]]$cvSet)
      testIndices = lapply(trainingSet, function(s) s[[i]]$indices)
      
      # Use ROC AUC to find the best eval function
      scores <- func(datasets, trains, testIndices)
      score <- scores$rocAUC
      
      if (score == 0) {
        stop("score returned 0");
      }

      iterations <- iterations + 1
      accumulated_score <- accumulated_score + score
    }
    
    avg_score <- accumulated_score / iterations
    
    test_scores[[m]] <- avg_score
  }
  
  return(list(score=test_scores, best_eval_func_indices=best_eval_func_indices))
}


cross_validate = function(trainingSet, datasets, k, eval_func) {
  library(cvTools) #run the above line if you don't have this library
  
  accumulated_scores = matrix(0, 3, length(trainingSet))
  iterations = 0
  
  for(i in 1:k){

    trains = lapply(trainingSet, function(s) s[[i]]$cvSet)
    testIndices = lapply(trainingSet, function(s) s[[i]]$indices)
    
    # scores <- tryCatch(eval_func(datasets, trains, validations), error = function(e) { print(e); return(0)})
    scores <- eval_func(datasets, trains, testIndices)
    
    iterations <- iterations + 1
    
    #3 is the number of scores
    mScores <- matrix(as.numeric(scores), 3, length(trains))

    accumulated_scores <- accumulated_scores + mScores
    
  }
  
  avg_scores <- accumulated_scores / iterations
  
  rownames(avg_scores) <- names(unlist(scores))[1:dim(avg_scores)[1]]
  
  avg_scores
}