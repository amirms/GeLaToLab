#Input
# data
# k: the number of folds
# pred_func
# eval_func is a list of list of evaluation functions

nested_cross_validate = function(datasets, k, eval_funcs) {
  #install_packages("cvTools")
  library(cvTools) #run the above line if you don't have this library
  
  folds = list();
  
  for (i in 1:length(datasets)) {
    dataset = t(datasets[[i]])
    folds[[i]] <- cvFolds(NROW(dataset), K=k)
  }
  
  # perform nested cross validation
  #exclude set number 1
  validationsets=list()
  for (i in 1:length(datasets)) {
    dataset = datasets[[i]]
    fold = folds[[i]]
    validationset =  dataset[,fold$subsets[fold$which != 1]]
    validationsets[[i]] <- validationset
  }
  
  vfolds = list();
  
  for (i in 1:length(datasets)) {
    validationset = t(validationsets[[i]])
    vfolds[[i]] <- cvFolds(NROW(validationset), K=k)
  }
  
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
        trains = list()
        validations = list()
        
        for (l in 1:length(validationsets)) {
          vfold = vfolds[[l]]
          
          validationset = validationsets[[l]]
          
          trains[[l]] <- validationset[,vfold$subsets[vfold$which != m] ] #Set the training set
          validations[[l]] <- validationset[,vfold$subsets[vfold$which == m]] #Set the validation set
        }
        
        score <- tryCatch(func(datasets, trains, validations), error = function(e) { print(e); return(0)})
        
        # score <- func(datasets, trains, validations)
        
        if (score == 0) {
          next;
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
      trains = list()
      validations = list()
    
      for (j in 1:length(datasets)) {
        fold = folds[[j]]
        dataset <- datasets[[j]]
        trains[[j]] <- dataset[,fold$subsets[fold$which != i] ] #Set the training set
        validations[[j]] <- dataset[,fold$subsets[fold$which == i]] #Set the validation set
      }
      
      score <- tryCatch(func(datasets, trains, validations), error = function(e) { print(e); return(0)})
      
      if (score == 0) {
        next;
      }
      
      iterations <- iterations + 1
      accumulated_score <- accumulated_score + score
     
    }
    
    avg_score <- accumulated_score / k
    
    test_scores[[m]] <- avg_score
  }
  
  return(list(score=test_scores, best_eval_func_indices=best_eval_func_indices))
}


cross_validate = function(datasets, k, eval_func) {
  library(cvTools) #run the above line if you don't have this library
  
  folds = list();
  
  for (i in 1:length(datasets)) {
    dataset = t(datasets[[i]])
    folds[[i]] <- cvFolds(NROW(dataset), K=k)
  }
  
  accumulated_scores = NULL
  iterations = 0
  
  for(i in 1:k){
    trains = list()
    validations = list()
    
    for (j in 1:length(datasets)) {
      fold = folds[[j]]
      dataset <- datasets[[j]]
      trains[[j]] <- dataset[,fold$subsets[fold$which != i] ] #Set the training set
      validations[[j]] <- dataset[,fold$subsets[fold$which == i]] #Set the validation set
    }
    
    scores <- tryCatch(eval_func(datasets, trains, validations), error = function(e) { print(e); return(0)})
    
    iterations <- iterations + 1
    if (is.null(accumulated_scores)) {
      accumulated_scores = rep(0, length(scores))
    }
    accumulated_scores <- accumulated_scores + scores
    
  }
  
  avg_scores <- accumulated_scores / iterations
  avg_scores
  
}

