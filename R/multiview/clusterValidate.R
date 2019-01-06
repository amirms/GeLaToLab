singleview.clusterer <- function(datasets, index, kernel_funcs) {
  mydata = datasets[[index]]
  
  #Build the dendrogram from data set module names
  priori.decomp <- build.dendrogam(rownames(mydata$mtx))
  
  best_kernel_func = NULL
  best_kernel_func_index = 1
  best_score = Inf
  
  for (i in 1:length(kernel_funcs)) {
    kernel_func <- kernel_funcs[[i]]
    
    # print(kernel_func)
    
    K <- tryCatch(kernel_func(mydata), error=function(e) e)
    
    stopifnot(all(rownames(K) == rownames(mydata)))
    
    if(inherits(K, "error")) next
    
    result <- measure_cluster_analysis(K, priori.decomp)
    
    score <- result$diff
    
    # XXX the smaller the score, the better the score
    if (score < best_score){
      best_score = score
      best_kernel_func = kernel_func
      best_kernel_func_index = i
    }
    
  }
  
  return(list(score=best_score, kernel_func = best_kernel_func, kernel_func_index=best_kernel_func_index))
}

measure_cluster_analysis <- function(kernel, priori.decomp){
  #compute distance from kernel
  myDist <- squared.euclidean.distance.of.kernel.matrix(kernel)
  myDist <- as.dist(myDist)
  
  # pinned it to complete linkage
  clusters <- hclust(myDist, method = 'complete')
  
  # compute tree distance 
  treeDistance = compute_tree_edit_distance_for_hc(clusters, priori.decomp$graph)
  
  clusters.tree <- ape::as.phylo(clusters)
  priori.tree <- priori.decomp$tree
  path.difference <- phangorn::path.dist(clusters.tree, priori.tree, check.labels = T)
  
  return(list(baker=0, cophcor=0, Bk=0, diff=path.difference, mojosim = 0, treeDistance = treeDistance))
}


multiview.clusterer <- function(cfgsim_kernel_func, freqsim_kernel_func, lexsim_kernel_func, fuse_multi_view_func, datasets){
  cfg <- datasets[[1]]
  freq <- datasets[[2]]
  lex <- datasets[[3]]
  
  #Build the dendrogram from data set module names
  priori.decomp <- build.dendrogam(rownames(cfg$mtx))
  
  cfgsim_kernel <- cfgsim_kernel_func(cfg)
  freqsim_kernel <- freqsim_kernel_func(freq)
  lexsim_kernel <- lexsim_kernel_func(lex)
  
  fused_Ks <- fuse_multi_view_func(list(cfgsim_kernel, freqsim_kernel, lexsim_kernel))
  
  scores <- lapply(fused_Ks, function(K) measure_cluster_analysis(K, priori.decomp))
  
  return(c(scores[[1]]$diff, scores[[2]]$diff, scores[[3]]$diff))
}

perform.clustering <- function(prname, rootFolder="org", pattern = "*.java"){
  require(GeLaToLab)
  require(proxy)
  
  
  setwd("~/workspace")
  
  #Read the set of classnames for running the experiment
  # classnames <- unlist(read.table(paste("benchmark", prname , paste("MULTIVIEW", "classnames.txt" ,sep="/") , sep="/")) )
  
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
  
  no_transactions <- colSums(freq)
  freq <- freq[, which(no_transactions <= 30)]
  
  #Load the bag of words
  BoW <- load_BoW(prname)
  # apply tf-idf mechanism and eliminate features that are lower than some threshold,
  # then remove those features from BoW
  x <- apply_tf_idf(BoW)
  dimnames(x) <- dimnames(BoW)
  BoW <- x

  # Clean up BoW
  no_words <- colSums(BoW)
  BoW <- BoW[, which(no_words > 0)]
  
  no_words_in_document <- rowSums(BoW)
  BoW <- BoW[which(no_words_in_document > 0),]
  
  #LOAD the text of the source code
  #Initialize the text of source files
  setwd(paste("benchmark", prname, sep="/"))
  txts <- read.text.directory(rootFolder, pattern)
  txts <- txts[order(names(txts))]
  
  
  print("dimensions before intersection")
  print(dim(cfg))
  print(dim(freq))
  print(dim(BoW))
  print(length(txts))
  #INTERSECT
  names <- intersect_all(rownames(cfg), rownames(freq), rownames(BoW), names(txts))
  
  cfg <- cfg[names, names]
  
  freq <- freq[names,]
  freq <- freq[, colSums(freq) > 0]
  
  BoW <- BoW[names,]
  BoW <- BoW[, colSums(BoW) > 0]
  
  txts <- txts[names]
  
  print("dimensions after intersection")
  
  print(dim(cfg))
  print(dim(freq))
  print(dim(BoW))
  print(length(txts))
  
  
  #lexsim kernels
  #cosine normalized linear kernel
  cosine_kernel_func <- function(x) { cos.sim(t(x$mtx))}
  
  # polnomial degree: 1,2,3,4,5
  polynomial_params <- c(1,2,3,4,5)
  polynomial_kernel_func <- lapply(polynomial_params, function(p) {function(x) polynomial.kernel(x$mtx, p)})
  
  #gaussian parameter: 10^(-5),10^(-4),10^(-3),10^(-2),10^(-1),10^(0),10^(1),10^(2)
  gaussian_params <- c(10^(-5),10^(-4),10^(-3),10^(-2),10^(-1),10^(0),10^(1),10^(2))
  gaussian_kernel_func <- lapply(gaussian_params, function(p) {function(x) gaussian.kernel(x$mtx, p)})
  
  # string kernels
  
  # p_spectrum parameters: 1,2,5,10,15,20
  p_spectrum_params <- c(1,2,5,10,15,20)
  p_spectrum_kernel_func <- lapply(p_spectrum_params, function(p) {function(x) spectrum.string.kernel(x$txt, p)})
  
  # constant kernel
  constant_kernel_func <- function(x) { constant.string.kernel(x$txt) }
  
  #exponential decay parameter: 1,2,5,10,15,20
  # what is this about?!?
  exponential_decay_params <- c(1,2,5,10,15,20)
  exponential_decay_kernel_func <- lapply(exponential_decay_params, function(p) {function(x) exponential.decay.kernel(x$txt, p)})
  
  #exponential diffusion parameter: 10^(-5),10^(-4),10^(-3),10^(-2),10^(-1),10^(0),10^(1),10^(2)
  exponential_diffusion_params <- c(10^(-5),10^(-4),10^(-3),10^(-2),10^(-1),10^(0),10^(1),10^(2))
  exponential_diffusion_kernel_func <- lapply(exponential_diffusion_params, function(p) {function(x) calc.diffusion.kernel(x$mtx, p, TRUE)})
  
  #laplacian exponential diffusion parameter: 10^(-5),10^(-4),10^(-3),10^(-2),10^(-1),10^(0),10^(1),10^(2)
  # what is this about?!?
  laplacian_exponential_diffusion_params <- c(1,2,5,10,15,20)
  laplacian_exponential_diffusion_kernel_func <- lapply(laplacian_exponential_diffusion_params, 
                                                        function(p) {function(x) compute.exponential.diffusion.kernel(x$mtx, p)})
  
  # commute time kernel
  commute_time_kernel_func <- function(x) {
    x <- x$mtx
    if (!isSymmetric(x)){
      x <- (x + t(x))/2
    }
    
    d <- apply(abs(x),1,sum)
    D <- diag(d)
    
    compute.avg.commute.time.kernel(x, D)
  }
  
  # exponential_decay_kernel_func
  lex_kernel_parameters_list <- list(0, polynomial_params, p_spectrum_params, 0)
  lexsim_kernel_funcs <- unlist(list(cosine_kernel_func, polynomial_kernel_func, p_spectrum_kernel_func, constant_kernel_func))
  
  freq_kernel_parameters_list <- list(polynomial_params, gaussian_params)
  freqsim_kernel_funcs <- unlist(list(polynomial_kernel_func, gaussian_kernel_func))
  
  #commute_time_kernel_func
  cfg_kernel_parameters_list <- list(exponential_diffusion_params, laplacian_exponential_diffusion_params)
  cfgsim_kernel_funcs <- unlist(list(exponential_diffusion_kernel_func, laplacian_exponential_diffusion_kernel_func))
  
  add_kernel_func <- function(Ks) {
    r <- add.kernels(Ks)
    list(r, r, r)
  }
  
  product_kernel_func <- function(Ks) {
    r <- product.kernels(Ks)
    list(r, r, r)
  }
  
  #perform nested cross validation
  datasets = list(cfg=list(mtx=cfg), freq=list(mtx=freq), lex=list(mtx=BoW, txt=txts))

  cfg_single_view <- singleview.clusterer(datasets, 1, cfgsim_kernel_funcs)
  freq_single_view <- singleview.clusterer(datasets, 2, freqsim_kernel_funcs)
  lex_single_view <- singleview.clusterer(datasets, 3, lexsim_kernel_funcs)
  
  
  cfgsim <- cfg_single_view$score
  freqsim <- freq_single_view$score
  lexsim <- lex_single_view$score
  
  best_cfg_eval_func <- cfg_single_view$kernel_func
  best_freq_eval_func <- freq_single_view$kernel_func
  best_lex_eval_func <- lex_single_view$kernel_func
  
  cotraining_kernel_func <- function(Ks) {
    cotraining(Ks, 50)
  }
  
  
  MKL_multiview_fuse_funcs <- list(add_kernel_func, product_kernel_func, cotraining_kernel_func, rgcca_func)
  MKL.multiview.clusterers <- lapply(MKL_multiview_fuse_funcs, function(MKL_multiview_fuse_func) {
    function(datasets)
      multiview.clusterer(best_cfg_eval_func, best_freq_eval_func, best_lex_eval_func, MKL_multiview_fuse_func, datasets)
  })
  
  MKL_add <- MKL.multiview.clusterers[[1]](datasets)
  MKL_product <-  MKL.multiview.clusterers[[2]](datasets)
  co_training <-  MKL.multiview.clusterers[[3]](datasets)
  k_cca <-  MKL.multiview.clusterers[[4]](datasets)

  result = list(cfgsim=cfgsim, freqsim=freqsim, lexsim=lexsim, MKL_add=min(MKL_add), MKL_product=min(MKL_product), co_training=min(co_training), kcca=min(k_cca))
  
  setwd("~/workspace")
  
  #Create results directory, if it doesn't exist
  dir.create(file.path(getwd(), paste("benchmark", prname, "Multiview/Results", sep="/")), showWarnings = FALSE)
  
  print_clustering_results(prname, result, txt.file = "Multiview/Results/Clustering.txt", rnd=2)
  
  #FIND BEST EVAL FUNCTIONS
  best_cfg_eval_func_idx <- cfg_single_view$kernel_func_index
  best_freq_eval_func_idx <- freq_single_view$kernel_func_index
  best_lex_eval_func_idx <- lex_single_view$kernel_func_index
  
  # evalFuncToString = function(func_index, type, parameters)
  cfg_eval_string <- evalFuncToString(best_cfg_eval_func_idx, "cfg", cfg_kernel_parameters_list)
  write(cfg_eval_string, file = paste("benchmark", prname, "Multiview/Results/CFG_EVAL.txt", sep="/"))
  
  freq_eval_string <- evalFuncToString(best_freq_eval_func_idx, "freq", freq_kernel_parameters_list)
  write(freq_eval_string, file = paste("benchmark", prname, "Multiview/Results/FREQ_EVAL.txt", sep="/"))
  
  lex_eval_string <- evalFuncToString(best_lex_eval_func_idx, "lex", lex_kernel_parameters_list)
  write(lex_eval_string, file = paste("benchmark", prname, "Multiview/Results/FREQ_EVAL.txt", sep="/"))
  
  # result = list(cfgsim=cfgsim, freqsim=freqsim, lexsim=lexsim, MKL_add=min(MKL_add), MKL_product=min(MKL_product), co_training=min(co_training), kcca=min(k_cca))
  all_string <- paste(round(cfgsim,2),cfg_eval_string,round(freqsim,2), freq_eval_string, round(lexsim,2), lex_eval_string, round(min(MKL_add),2), round(min(co_training),2), round(min(k_cca),2), sep="&")
  write(all_string, file = paste("benchmark", prname, "Multiview/Results/All_EVAL.txt", sep="/"))
}