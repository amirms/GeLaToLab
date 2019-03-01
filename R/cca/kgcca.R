#Input: a list of kernel matrices
kgcca = function(Ks, ncomp=rep(3,4)) {
  
  require(RGCCA)
  
  pcvs = eigFeatures(Ks)
  
  pcvs[[4]] <- cbind(pcvs[[1]], pcvs[[2]], pcvs[[3]])
  # 
  # #TODO take the average of the components?
  # 
  C = matrix(c(0,0,0,1,0,0,0,1,0,0,0,1,1,1,1,0), 4, 4)
  # 
  # tau = rep(1,4)
  tau = c(rep(1,3), 0) # co-intertia
  # tau = rep(0,4) # Generalizied CCA
  return(rgcca(pcvs, C= C, tau=tau, ncomp=ncomp, scheme="factorial"))
}

eigFeatures <- function(Ks,features = 0, th = 1e-4) {
  pcvs = list();
  
  p_features <- features
  
  for(i in 1:length(Ks)){
    features <- p_features
    km <- Ks[[i]]
    km <- as.matrix(km)
    m <- nrow(km)
    
    #t(t(m/rowSums(m)) + rowSums(t(m)))
    
    ## center kernel matrix - correct
    kc <- t(t(km - colSums(km)/m) -  rowSums(km)/m) + sum(km)/m^2
    
    ## compute eigenvectors
    res <- eigen(kc/m,symmetric=TRUE)
    
    if(features == 0)
      features <- sum(res$values > th)
    else {
      if(res$values[features] < th)
        stop(paste("eigenvalues of the kernel matrix are below threshold!"))
    }
    
    print("features")
    print(features)
    
    pcvs[[i]] <- t(t(res$vectors[,1:features])/sqrt(res$values[1:features]))
  }
  
  return(pcvs)
}

rgcca_func <- function(Ks) {
  r <- kgcca(Ks)
  
  results = list()
  
  for (i in 1:length(Ks)) {
    Xr <- r$Y[[4]] #%*% diag(r$AVE$AVE_X[[4]]) #t(r$a[[4]])
    Xr <- apply(Xr, 2, normalize_min_zero_unit)
    
   D <- as.matrix(dist(Xr))
   sim <- 1 / (D / max(D))
   
   # print("dim(sim)")
   # print(dim(sim))
   
   results[[i]] <- sim
   
   # print("dim(Ks[[i]])")
   # print(dim(Ks[[i]]))
   
   dimnames(results[[i]]) <- dimnames(Ks[[i]])
  }
  
  # avgXr <- matrix(0, nrow = dim(Xrs[[1]])[1], ncol = dim(Xrs[[1]])[2])
  # for (i in 1:length(Xrs)) {
  #   avgXr <- avgXr + Xrs[[i]]
  # }
  #   
  # Xr <- (1/length(Ks)) * avgXr

  
  # D <- dist(Xr)
  # sim <- 1 / (D / max(D))
  # results <- list()
  # 
  # for (i in 1:length(Ks)) {
  #   results[[i]] <- sim
  #   dimnames(results[[i]]) <- dimnames(Ks[[i]])
  # }
  
  results
}

normalize_min_zero_unit <- function(x) {
  (x- min(x) )/ (max(x) - min(x))
}