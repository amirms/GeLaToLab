cutree.k.dendrogram <- function (dend, k, dend_heights_per_k = NULL, use_labels_not_values = TRUE, 
          order_clusters_as_data = TRUE, warn = dendextend_options("warn"), 
          ...) 
{
  k <- as.integer(k)
  if (!is.natural.number(k)) 
    stop(paste("k must be a natural number!  The k you used (", 
               k, ") is not a natural number"))
  if (k > nleaves(dend)) 
    stop(paste("elements of 'k' must be between 1 and", nleaves(dend)))
  if (k == 1L) {
    h_to_use <- attr(dend, "height") + 1
    cluster_vec <- cutree.h.dendrogram(dend, h = h_to_use, 
                                        use_labels_not_values = use_labels_not_values, order_clusters_as_data = order_clusters_as_data, 
                                        ...)
    return(cluster_vec)
  }
  if (k == nleaves(dend)) {
    labels_dend <- labels(dend)
    cluster_vec <- 1:length(labels_dend)
    names(cluster_vec)[order.dendrogram(dend)] <- labels_dend
    return(cluster_vec)
  }
  if (is.null(dend_heights_per_k)) {
    dend_heights_per_k <- heights_per_k.dendrogram(dend)
  }
  height_for_our_k <- which(names(dend_heights_per_k) == k)
  if (length(height_for_our_k) != 0) {
    h_to_use <- dend_heights_per_k[height_for_our_k]
    cluster_vec <- cutree.h.dendrogram(dend, h = h_to_use, 
                                        use_labels_not_values = use_labels_not_values, order_clusters_as_data = order_clusters_as_data, 
                                        ...)
  }
  else {
    cluster_vec <- rep(NA, nleaves(dend, method = "order"))
    if (warn) {
      warning("Couldn't cut the dend - returning NA.")
      k_s <- as.numeric(names(dend_heights_per_k))
      if (k > max(k_s) || k < min(k_s)) {
        range_for_clusters <- paste("[", paste(range(k_s), 
                                               collapse = "-"), "]", sep = "")
        warning(paste("No cut exists for creating", k, 
                      "clusters.  The possible range for clusters is:", 
                      range_for_clusters))
      }
      else {
        warning(paste("You (probably) have some branches with equal heights so that there exist no height(h) that can create", 
                      k, " clusters"))
      }
    }
  }
  return(cluster_vec)
}

cutree.h.dendrogram <- function (dend, h, order_clusters_as_data = TRUE, use_labels_not_values = TRUE, 
          warn = dendextend_options("warn"), ...) 
{
  if (missing(h)) 
    stop("h is missing")
  if (length(h) > 1) {
    if (warn) 
      warning("h has length > 1 and only the first element will be used")
    h <- h[1]
  }
  if (h < 0) {
    labels_dend <- labels(dend)
    cluster_vec <- 1:length(labels_dend)
    # names(cluster_vec)[order.dendrogram(dend)] <- labels_dend #FIXME may be pass the order as an argument
    names(cluster_vec) <- labels_dend
    return(cluster_vec)
  }
  if (use_labels_not_values) {
    FUN <- labels
  }
  else {
    FUN <- order.dendrogram
  }
  names_in_clusters <- cut_lower_fun(dend, h, FUN)
  number_of_clusters <- length(names_in_clusters)
  number_of_members_in_clusters <- sapply(names_in_clusters, 
                                          length)
  cluster_vec <- rep(rev(seq_len(number_of_clusters)), times = number_of_members_in_clusters)
  if (h > attr(dend, "height")) 
    cluster_vec <- rep(1L, length(cluster_vec))
  names(cluster_vec) <- unlist(names_in_clusters)
  # clusters_order <- order.dendrogram(dend)
  # if (order_clusters_as_data) {
  #   if (!all(clusters_order %in% seq_along(clusters_order))) {
  #     if (warn) {
  #       warning("rank() was used for the leaves order number! \nExplenation: leaves tip number (the order), and the ranks of these numbers - are not equal.\n  The dend was probably subsetted, pruned and/or merged with other dends- and now the order \n labels don't make so much sense (hence, the rank on them was used).")
  #       warning("Here is the cluster order vector (from the dend tips) \n", 
  #               paste(clusters_order, collapse = ", "), "\n")
  #     }
  #     clusters_order <- rank(clusters_order, ties.method = "first")
  #   }
  #   cluster_vec <- cluster_vec[order(clusters_order)]
  # }
  dend_size <- nleaves(dend, method = "order")
  if (number_of_clusters == dend_size) 
    cluster_vec[seq_len(dend_size)] <- seq_len(dend_size)
  return(cluster_vec)
}