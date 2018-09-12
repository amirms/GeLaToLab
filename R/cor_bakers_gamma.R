#' @title Bakers Gamma for two k matrices
#' @description Bakers Gamma for two k matrices
#' @param k_matrix_dend1 a matrix of k cluster groupings from a dendrogram
#' @param k_matrix_dend2 a (second) matrix of k cluster groupings from a dendrogram
#' @param to_plot logical (FALSE). Should a scaterplot be plotted, showing the
#' correlation between the lowest shared branch between two items in the two
#' compared trees.
#' @param ... not used
#' @seealso
#' \link{cor_bakers_gamma}
#' @return 
#' Baker's Gamma coefficient.
bakers_gamma_for_2_k_matrix <- function(k_matrix_dend1, k_matrix_dend2, to_plot = FALSE)
{
  
  if(dim(k_matrix_dend1)[1] != dim(k_matrix_dend2)[1]) stop("The k_matrixes seems to show a different number of items - we can not compare trees in this case!")
  if(!all(sort(rownames(k_matrix_dend1)) == sort(rownames(k_matrix_dend2)))) 
  {		# we are using "sort" since the rownames may be of different order - depending on the way the two trees were constructed.
    print(paste("Item names (rownames) of k_matrix_dend1:" ,rownames(k_matrix_dend1)))
    print(paste("Item names (rownames) of k_matrix_dend2:" ,rownames(k_matrix_dend2)))
    stop("The k_matrixes seems to have different item names - \n we can not compare trees in this case! \n Consider using use_labels_not_values = T (or F) in cutree")
  }
  
  all_combinations_of_items <- t(combn(seq_len(dim(k_matrix_dend1)[1]),2))
  number_of_combinations_of_items <- dim(all_combinations_of_items)[1]
  
  
  cor_mat <- matrix(0, number_of_combinations_of_items, 2)
  
  
  for(i in seq_len(number_of_combinations_of_items))
  {
    item_id_1_name <- rownames(k_matrix_dend1)[all_combinations_of_items[i,1]]
    item_id_2_name <- rownames(k_matrix_dend1)[all_combinations_of_items[i,2]]
    # The names must be identical in both trees - we've made a stopping rule for that already.
    item1 <- k_matrix_dend1[item_id_1_name,]
    item2 <- k_matrix_dend1[item_id_2_name,]
    # print(paste(i, item1, item2))
    cor_mat[i,1] <- lowest_common_branch(item1,  item2)
    
    item1 <- k_matrix_dend2[item_id_1_name,]
    item2 <- k_matrix_dend2[item_id_2_name,]
    cor_mat[i,2] <- lowest_common_branch(item1,  item2)
  }
  
  COR_object <- cor(cor_mat[,1], cor_mat[,2], method = "spearman")
  if(is.na(COR_object)) COR_object <- 1 # because this is NA only if the two vectors are identical
  #  that happens only when the two trees have only leaves. Which is o.k. to define as correlation of 1.
  
  if(to_plot) {
    plot(cor_mat, sub = paste("COR =", round(COR_object,4)))
  }	
  
  return(COR_object)
}

lowest_common_branch <- function(item1, item2,...)
{
  # for two rows of cluster belonging,
  # this finds all the cases where the two items are identical (e.g: belong to the same branch)
  # Then it takes the most extreme value (lowest branch) where this happens
  # and finally - it extracts the level of that lowest brnach (names) and turn it to be numeric
  value <- as.numeric(names(tail(which(item1 == item2),1))	)
  if(length(value) == 0) value <- 0
  return(value)
}

cor_bakers_gamma.dendrogram <- function(dend1, dend2, use_labels_not_values = TRUE, to_plot = FALSE, warn = dendextend_options("warn"), ...)
{
  tree1_heights_per_k <- heights_per_k.dendrogram(dend1)
  tree2_heights_per_k <- heights_per_k.dendrogram(dend2)
  
  if (length(tree1_heights_per_k) > length(tree2_heights_per_k))
    dend_heights_per_k <- tree2_heights_per_k
  else 
    dend_heights_per_k <- tree1_heights_per_k
  
  # k <- everything except 1 and nleaves
  ks_to_use <- as.integer(names(dend_heights_per_k[which(!(names(dend_heights_per_k) %in% c(1,nleaves(dend1))))]))
  
  cutree_tree1 <- lapply(ks_to_use, function(k) cutree.k.dendrogram(dend1, k = k))
  cutree_tree2 <- lapply(ks_to_use, function(k) dendextend::cutree(dend2, k = k))
  
  # k_matrix_dend1 <- cutree(dend1, k = 1:nleaves(dend1), use_labels_not_values=use_labels_not_values,warn=warn,...)
  # k_matrix_dend2 <- cutree(dend2, k = 1:nleaves(dend2), use_labels_not_values=use_labels_not_values,warn=warn,...)
  k_matrix_dend1 <- do.call(cbind, cutree_tree1)
  k_matrix_dend2 <- do.call(cbind, cutree_tree2)
  
  colnames(k_matrix_dend1) <- ks_to_use
  colnames(k_matrix_dend2) <- ks_to_use
  
  bakers_gamma <- bakers_gamma_for_2_k_matrix(k_matrix_dend1, k_matrix_dend2, to_plot = to_plot)
  return(bakers_gamma)
}