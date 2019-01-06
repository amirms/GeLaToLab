query<- "The class that handles GUI elements such as search dialog"

findCrossModal <- function(prname, query, single_view_results){
  
  require(proxy)
  require(NMF)
  require(lsa)
  
  setwd("~/workspace")
  
  #Load the adjacency matrix
  extensions= c("java/", "org/xml/", "javax/")
  cfg <- import.bunch.matrix(paste("benchmark", prname ,"dep_graph.txt", sep="/"), exclude.ext=extensions)
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
  original_BoW <- load_BoW(prname)
  # apply tf-idf mechanism and eliminate features that are lower than some threshold,
  # then remove those features from BoW
  x <- apply_tf_idf(original_BoW)
  dimnames(x) <- dimnames(original_BoW)
  BoW <- x

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
  
  
  #FIND BEST EVAL FUNCTIONS
  best_cfg_eval_func_idx <- single_view_results$best_eval_func_indices[[1]]
  best_freq_eval_func_idx <- single_view_results$best_eval_func_indices[[2]]
  best_lex_eval_func_idx <- single_view_results$best_eval_func_indices[[3]]
  
  best_cfg_eval_func <- cfgsim_kernel_funcs[[best_cfg_eval_func_idx]]
  best_freq_eval_func <- freqsim_kernel_funcs[[best_freq_eval_func_idx]]
  best_BoW_eval_func <- lexsim_kernel_funcs[[best_lex_eval_func_idx]]
  
  # return(single_view_results)
  
  original_BoW <- original_BoW[,colnames(BoW)]
  
  query_vector <- decomposeQuery(query, original_BoW)
  
  # 
  cfgsim <- best_cfg_eval_func(cfg)
  freqsim <- best_freq_eval_func(freq)
  lexsim <- best_BoW_eval_func(BoW)
  Ks <- list(cfgsim, freqsim, lexsim)
  
  newBoW <- rbind(BoW, query_vector)
  alllexsim <- best_BoW_eval_func(newBoW)
  
  pcvs <- computeKernalPCA(Ks)
  
  
  # In single view setting: 
  qv <- projectQueryIntoFeatureSpace(as.vector(alllexsim[dim(alllexsim)[1],1:(dim(alllexsim)[1] -1)]), pcv = pcvs[[3]], k = lexsim) 
  
  singleXr <- rbind(pcvs[[3]], qv)
  singleXr <- apply(singleXr, 2, normalize_min_zero_unit)
  d <- as.matrix(dist(singleXr))
  ranks <- sort(d[dim(alllexsim)[1],], decreasing = F)[1:10]
  ranks <- ranks[-1]
  
  
  # In multi view setting: 
  r <- computeCanonicalMatching(pcvs)
  lowqv <- qv %*% r$a[[3]] %*% diag(sqrt(r$AVE$AVE_X[[3]]))

  Xr <- r$Y[[4]]
  Xr <- rbind(Xr, lowqv)
  
  Xr <- apply(Xr, 2, normalize_min_zero_unit)
  d <- as.matrix(dist(Xr))
  ranks <- sort(d[dim(alllexsim)[1],], decreasing = F)[1:10]
  ranks <- ranks[-1]
}

center_scale <- function(x) {
  as.vector(scale(as.vector(x), scale = FALSE))
}


computeKernalPCA = function(Ks, features = 0, th = 1e-4) {
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
        stop(paste("eigenvalues of the kernel matrix are below threshold!"))
    }
    
    print("features")
    print(features)
    
    pcvs[[i]] <- t(t(res$vectors[,1:features])/sqrt(res$values[1:features]))
    rownames(pcvs[[i]]) <- rownames(Ks[[i]])
  }
  
  return(pcvs)
}

# Input:
# x: query similarity vector
# pcv: feature space
# original query

projectQueryIntoFeatureSpace  <- function(x, pcv, k) {
    x  <- if (is.vector(x)) t(t(x)) else if (!is(x,"list")) x <- as.matrix(x)
    
    if (is.vector(x) || is.data.frame(x))
      x <- as.matrix(x)
    if (!is.matrix(x) && !is(x,"list")) stop("x must be a matrix a vector, a data frame, or a list")
    
    if(is(x,"matrix"))
    {
      n <- nrow(x)
      m <- nrow(k)
    }
    else
    {
      n <- length(x)
      m <- length(k)
    }
    
    knc <- x
    ka <- k

    ## center
    ret <- t(t(knc - rowSums(knc)/m) - rowSums(ka)/m) + sum(ka)/(m*n) #[1xn]
    
    print(dim(ret))
    print(dim(pcv))
    
    qk <- t(ret) %*% pcv
    return(qk)
}


computeCanonicalMatching <- function(pcvs, ncomp=rep(2,4)) {
  require(RGCCA)
  
  pcvs[[4]] <- cbind(pcvs[[1]], pcvs[[2]], pcvs[[3]])
  
  C = matrix(c(0,0,0,1,0,0,0,1,0,0,0,1,1,1,1,0), 4, 4)
  
  # generalized CCA
  return(rgcca(pcvs, C= C, tau=rep(0,4), scheme="centroid", ncomp=ncomp))
}



tf_idf <- function(x){
  #diagonal matrix for term weighings
  doc.freq <- colSums(x>0)
  doc.freq[doc.freq == 0] <- 1
  
  w <- log(nrow(x)/doc.freq)
  R <- diag(w)
  
  R
}


decomposeQuery <- function(query, original_BoW){
  txt <- strip.text(query)
  
  mydata <- prepare.natural.lang.list(txt, "english")
  
  mydata <- mydata[-which(unlist(lapply(mydata, function(data) length(data) ==0)))]
  
  bow.list <- unlist(make.BoW.list(mydata))
  
  x <- rep(0, dim(original_BoW)[2])
  names(x) <- colnames(original_BoW)
  
  for (i in 1:length(bow.list)){
    name <- names(bow.list)[i]
    
    if (name %in% names(x)){
      x[name] <- bow.list[name]
    }
  }
  
  R <- tf_idf(original_BoW)
  
  z<- x %*% R
  
  z<- as.vector(z)
  
  names(z) <- names(x)
  z
}
