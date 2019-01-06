# Input
# X is a membership matrix
# K is a similarity matrix
# pIndicex: the object for which we would like to predict

compute_indirect_similarity = function(objectIndex, featureIndices, all_sim, X, k_neareset_neighbor=0){
  # print("objectIndex")
  # print(objectIndex)
  noOfObjects = dim(X)[1];
  
  if (noOfObjects > 100) {
    k_neareset_neighbor = 100
  } 
  
  p_predict <- rep(0, length(featureIndices))
  
  for (i in 1:length(featureIndices)) {
    p_predict[i] <- cf_predict(objectIndex, featureIndices[i], all_sim, X, k_neareset_neighbor)
  }
  
  p_predict
}

cf_predict <- function(obj_index, item_index, sim, data, k=0) {
  # k nearest neighbor
  if (k==0){
    k = dim(sim)[1]
  }
  
  if (obj_index > dim(sim)[1]){
    print("obj_index")
    print(obj_index)
    stop("object index too big")
  }

    obj_sim_info <- sort(sim[obj_index,], decreasing = T, index.return=TRUE)
  
  obj_sim <- obj_sim_info$x
  obj_sim_Idx <- obj_sim_info$ix
  
  sum_of_sim_item = 0
  sum_of_sim = 0
  
  for (j in 1:k){
    current_sim <- obj_sim[j]
    current_sim_Idx <- obj_sim_Idx[j]
    
    if (!is.finite(current_sim_Idx) || !is.finite(obj_index)){
      stop("something went wrong")
    }

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

calculate_ROCAUC = function(truth, prediction){
  library(PRROC)
  
  stopifnot(any(truth > 0))
  
  fg <- prediction[truth==1]
  bg <- prediction[truth==0]
  roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  roc$auc
}

calculate_PRAUC = function(truth, prediction){
  library(PRROC)
  
  stopifnot(any(truth > 0))
  
  fg <- prediction[truth==1]
  bg <- prediction[truth==0]
  pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)

  pr$auc.integral
}


calculate_max_F1 = function(truth, prediction) {
  library(Laurae)
  
  scores <- get.max_f1(prediction, truth)
  
  scores[1]
}