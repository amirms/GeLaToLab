#Input: a list of kernel matrices
kgcca = function(Ks, features = 0, th = 1e-4) {
  
  require(RGCCA)
  
  pcvs = list();
  
  p_features <- features
  
  for(i in 1:length(Ks)){
    features <- p_features
    km <- Ks[[i]]
    km <- as.matrix(km)
    m <- nrow(km)
    
    ## center kernel matrix
    kc <- t(t(km - colSums(km)/m) -  rowSums(km)/m) + sum(km)/m^2
    
    ## compute eigenvectors
    res <- eigen(kc/m,symmetric=TRUE)
    
    if(features == 0)
      features <- sum(res$values > th)
    else {
      if(res$values[features] < th)
        warning(paste("eigenvalues of the kernel matrix are below threshold!"))
    }
    
    pcvs[[i]] <- t(t(res$vectors[,1:features])/sqrt(res$values[1:features]))
    
  }
  
  return(rgcca(pcvs))
  
}

rgcca_func <- function(Ks) {
  r <- kgcca(Ks)
  
  results <- list()
  for (i in 1:length(Ks)) {
   Xr <- r$Y[[i]]  %*% t(r$a[[i]])
   #TODO make this a similarity matrix
   Xr <- apply(Xr, 2, normalize_min_zero_unit)
   results[[i]] <- Xr %*% t(Xr)
   dimnames(results[[i]]) <- dimnames(Ks[[i]])
  }
  
  results
}

normalize_min_zero_unit <- function(x) {
  (x- min(x) )/ (max(x) - min(x))
}