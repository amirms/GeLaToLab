# Input
# X is a membership matrix
# K is a similarity matrix
# pIndicex: the object for which we would like to predict

compute_indirect_similarity = function(feature_index, all_sim, X, k_neareset_neighbor=0){
  
  no_of_objects = dim(X)[1]
  no_of_features = dim(X)[2]
  
  if (no_of_features > 300) {
    k_neareset_neighbor = 300
  }
  
  p_predict <- rep(0, no_of_objects)
  
  for (i in 1:no_of_objects) {
    
    p_predict[i] <- cf_predict(i, feature_index, all_sim, X, k_neareset_neighbor)
  }
  
  p_predict
}

cf_predict <- function(obj_index, item_index, sim, data, k=0) {
  # k nearest neighbor
  if (k==0){
    k = dim(sim)[1]
  }
  
  obj_sim_info <- sort(sim[obj_index,], decreasing = T, index.return=TRUE)
  
  obj_sim <- obj_sim_info$x
  obj_sim_Idx <- obj_sim_info$ix
  
  sum_of_sim_item = 0
  sum_of_sim = 0
  
  for (j in 1:k){
    current_sim <- obj_sim[j]
    current_sim_Idx <- obj_sim_Idx[j]
    
    if (current_sim_Idx == obj_index)
      next;
    
    item_value <- data[current_sim_Idx, item_index]
    
    sum_of_sim_item <- sum_of_sim_item + current_sim * item_value
    sum_of_sim <- sum_of_sim + current_sim
  }
  
  pred_value = 0
  
  if (sum_of_sim > 0) {
    pred_value <- sum_of_sim_item / sum_of_sim
  }
  
  return(pred_value)
}

calculate_AUC = function(category, prediction){
  library(pROC)
  
  stopifnot(any(category > 0))
    
  roc_obj <- roc(category, prediction)
  auc(roc_obj)
}