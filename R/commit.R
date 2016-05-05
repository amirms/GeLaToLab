

#Remember: For "eclipse-jdt-core-3.8", "weka-3.6.11", the root position is the 2nd occurring root name
# For "weka-3.6.11" the root name is weka
extract.transactions <- function(prname, root = "org", rootpos = 1) {
  require(XML)
  
  setwd("~/workspace")
  
  #Load the log histoy
  log_data <- xmlParse(paste("benchmark", prname ,"log.xml", sep="/"))
  xml_data <- xmlToList(log_data)
  
  #simple  transaction extraction
  paths <- lapply(xml_data, function(l) l$paths)
  
  texts <- lapply(paths, function(path) unlist(lapply(path, function(p) p$text)))  
  
  unique_texts <- unique(unlist(texts))
  
  change_freq <- matrix(0, nrow = length(unique_texts), ncol = length(texts))
  
  rownames(change_freq) <- unique_texts
  
  for (i in 1:length(texts))
    change_freq[texts[[i]], i] <- 1
  
  java_files <- grep(pattern = "\\.java$", unique_texts, ignore.case = FALSE, perl = FALSE, value = TRUE)
  
 
  
  #Rename to reflect root package structure
  #java_files
  
  #need to fix this - find the first occurrence of root, then get the substring starting from that index
#   java_files <- unlist(lapply(java_files, function(filename) paste(root, strsplit(filename, root)[[1]][2])))
  
  #Find the indices first
  indices <- unlist(lapply(java_files, function(filename) gregexpr(root, filename)[[1]][rootpos]))
  
  java_files <- java_files[which(indices > 0)]
  
  change_freq <- change_freq[java_files, ]
  
  changes <- colSums(change_freq)
  
  change_freq <- change_freq[, changes >= 2]

  
  #sbstring and rename java_files
  java_files <- unlist(lapply(java_files, function(filename) substring(filename, gregexpr(root, filename)[[1]][rootpos] ) ))
  
  #java_files <- gsub("/", ".", java_files, fixed=TRUE)
  
  rownames(change_freq) <- java_files

  unique_java_files <- unique(java_files)
  unique_change_freq <- matrix(0, length(unique_java_files), ncol(change_freq))  
  rownames(unique_change_freq) <- unique_java_files


  # merge duplicated rows in change_freq
  for(i in 1:length(java_files))
    unique_change_freq[java_files[i],] <- unique_change_freq[java_files[i],] + change_freq[java_files[i],]



  write.table(unique_change_freq, file = paste("benchmark", prname ,"mydata-change-freq-matrix.csv", sep="/"),
              row.names=TRUE, col.names=NA,sep=",", quote=FALSE)

unique_change_freq
  
}


find.best.commit.kernel <- function(prname, par){
  
  require(lsa)
  require(proxy)
  
  setwd("~/workspace")
  
  #Load the priori decomposition
  decomposition <- read.csv(paste("benchmark", prname ,"decomposition.csv", sep="/"), sep=",",  header = TRUE)
  priori.decomp <- decomposition$x
  names(priori.decomp) <- decomposition$X
  priori.decomp <- priori.decomp[which(names(priori.decomp) %in% classnames)] # find the ones in classnames
  priori.decomp <- priori.decomp[order(names(priori.decomp))]
  priori.decomp <- normalizeVector(priori.decomp)

  
  freq <- read.table(paste("benchmark", prname , "mydata-change-freq-matrix.csv", sep="/"), sep=",", row.names = 1, header = TRUE, check.names = FALSE)
  
  freq <- as.matrix(freq)

  #Fix freq
  freq <- freq[which(rownames(freq) %in% classnames),] # find the ones in classnames
  freq <- freq[order(rownames(freq)),]
    
  no_transactions <- colSums(freq)

  freq <- freq[, which(no_transactions <= 30)]

  priori.decomp.names <- names(priori.decomp)
  freq.names <- rownames(freq)
  
  for (i in 1:length(names(priori.decomp.names)))
    if (!(priori.decomp.names[i] == freq.names[i]))
      stop('wrong order of names')  
  
  priori.decomp <- normalizeVector(priori.decomp)
  
  print("printing the numer of groups:")
  print(max(priori.decomp))
  
  k <- max(priori.decomp)
  
  commit.kernel <- list()
  
#   commit.kernel[[length(commit.kernel) + 1]] <- as.matrix(dist(freq, method="Euclidean"))
#   commit.kernel[[length(commit.kernel) + 1]] <- as.matrix(dist(freq, method="cosine"))
  
  #Linear Kernel
#   commit.kernel[[length(commit.kernel) + 1]] <- linear.kernel(freq)

  degrees <- 1:7
  
  for(i in 1:length(degrees))
    commit.kernel[[length(commit.kernel) + 1]] <- polynomial.kernel(freq, degrees[i])

  for(i in 1:length(par))
    commit.kernel[[length(commit.kernel) + 1]] <- gaussian.kernel(freq, par[i])

  #compute the euclidean distance, and take the median as the parameter
  d <- dist(freq, method = "Euclidean", diag = FALSE, upper = FALSE)
  
  med <- median(d)
  
  print("median")
  print(med)
  
#   avg <- mean(d)
#   
#   print("mean")
#   print(avg)
  
  commit.kernel[[length(commit.kernel) + 1]] <- gaussian.kernel(freq, med)
#   commit.kernel[[length(commit.kernel) + 1]] <- gaussian.kernel(freq, avg)
  
  laps <- lapply(commit.kernel, function(kern) laplacian(kern, TRUE))
  
  print(length(laps))
  results=list()
  for (i in 1:length(laps)){
    print(paste("computing the spectral clustering for kernel matrix:", i))
    results[[i]] <- spectral.clustering(laps[[i]], k)
    
  }
  #   results <- lapply(laps, function(lap) {
  #     spectral.clustering(lap, k)})
  
  for (i in 1:length(results)) {
    names(results[[i]]) <- rownames(laps[[i]])
    results[[i]] <- normalizeVector(results[[i]])
    
  }
  
  N <- length(priori.decomp)
  
  mojoSims <- unlist(lapply(results, function(r) 1 - (compute.MoJo(r, priori.decomp)/N)))
  
  #   plot(x=par, y=mojoSims[-1])
  
  best = list()

  best.mjsim <- max(mojoSims[1:length(degrees)])

    index <- which(mojoSims[1:length(degrees)]==best.mjsim)
  
    variance <- var(mojoSims[1:length(degrees)])
    
    best[[length(best) + 1]] <- list(par=degrees[index], value = best.mjsim, variance = variance)

  

  #Gaussian
  best.mjsim <- max(mojoSims[(length(degrees)+1):(length(degrees)+length(par))])
  
  variance <- var(mojoSims[(length(degrees)+1):(length(mojoSims))])

  if (best.mjsim >= mojoSims[length(mojoSims)])
  {
  
    index <- which(mojoSims[(length(degrees)+1):(length(degrees)+length(par))]==best.mjsim)   

    best[[length(best) + 1]] <- list(par=par[index], value = best.mjsim, variance = variance)
  
  }
  else
  {
    
    best.mjsim <-  mojoSims[length(mojoSims)]
    
    best[[length(best) + 1]] <- list(par=med, value = best.mjsim, variance = variance)
    
  }

print(mojoSims)
  return(best)


  best.mjsim = max(mojoSims)

  
  #   return(mojoSims)
  #   return(which(mojoSims==best.mjsim))
  return(list(mojosim = mojoSims, indices= which(mojoSims==best.mjsim)))
  
}

test.best.commit.kernel <- function(indices = 0) {
# 
  projects <- list("apache-ant-1.9.3", "hadoop-0.20.2", "apache-log4j-1.2.17", "eclipse-jdt-core-3.8",
                   "jdom-2.0.5", "jedit-5.1.0", "jfreechart-1.2.0", "jhotdraw-7.0.6", "junit-4.12")#,"weka-3.6.11")
  
  names(projects) <- c("Apache Ant", "Apache Hadoop", "Apache Log4j", "Eclipse JDT Core", "JDOM", 
                    "JEdit", "JFreeChart", "JHotDraw", "JUnit")#, "Weka")
  
#   projects <- list("apache-ant-1.9.3", "apache-log4j-1.2.17", "eclipse-jdt-core-3.8")#,
#                   # "jdom-2.0.5", "jfreechart-1.2.0", "jhotdraw-7.0.6", "junit-4.12", "weka-3.6.11")
#   
#   names(projects) <- c("Apache Ant", "Apache Log4j", "Eclipse JDT Core")#, "JDOM", 
                       # "JFreeChart", "JHotDraw", "JUnit", "Weka")
  
#   projects <- list("weka-3.6.11")
#   names(projects) <- c("Weka")
 
  pars <- c(10^(-2), 10^(-1), 10^(0), 5, 10^(1), 10^(2))
  
  results <- lapply(projects, function(p) find.best.commit.kernel(p, pars))
  
  names(results) <- names(projects)
  
  # plot(c(min(decay),max(decay)),c(0,1),type='n',xlab="lambda",ylab='MoJoSim')
  
  
  
  #\cellcolor[gray]{0.8}
  
  if (indices == 0)
    indices <- 1:length(projects)
  
  str2tex = ""
  
  for (i in 1:length(results))
    str2tex <- paste(str2tex, print_latex(results[[i]], names(results[i])), "\n" )
  
  str2tex
  
  #Results
#   $`jhotdraw-7.0.6`
#   $`jhotdraw-7.0.6`$mojosim
#   [1] 0.2647059 0.2794118 0.5147059 0.3970588 0.3676471 0.3088235 0.3529412 0.3382353 0.3382353 0.3676471 0.3382353
#   
#   $`jhotdraw-7.0.6`$indices
#   [1] 3
#   
#   
#   $`junit-4.12`
#   $`junit-4.12`$mojosim
#   [1] 0.3275862 0.3448276 0.5517241 0.4827586 0.4310345 0.3793103 0.2758621 0.2931034 0.3620690 0.2931034 0.2758621
#   
#   $`junit-4.12`$indices
#   [1] 3
#   
#   
#   $`jfreechart-1.2.0`
#   $`jfreechart-1.2.0`$mojosim
#   [1] 0.2057416 0.2057416 0.7129187 0.4114833 0.3779904 0.3779904 0.3014354 0.3014354 0.2775120 0.3349282 0.2822967
#   
#   $`jfreechart-1.2.0`$indices
#   [1] 3
#   
#   
#   $`jdom-2.0.5`
#   $`jdom-2.0.5`$mojosim
#   [1] 0.4193548 0.4354839 0.4193548 0.4193548 0.3709677 0.4032258 0.3709677 0.4354839 0.4193548 0.3548387 0.4354839
#   
#   $`jdom-2.0.5`$indices
#   [1]  2  8 11
#   
#   
#   $`apache-log4j-1.2.17`
#   $`apache-log4j-1.2.17`$mojosim
#   [1] 0.3333333 0.3684211 0.4035088 0.4210526 0.3684211 0.3684211 0.3333333 0.3684211 0.3684211 0.3333333 0.3333333
#   
#   $`apache-log4j-1.2.17`$indices
#   [1] 4
#   
#   
#   $`jedit-5.1.0`
#   $`jedit-5.1.0`$mojosim
#   [1] 0.2012012 0.2492492 0.5645646 0.4864865 0.3483483 0.3903904 0.4384384 0.3093093 0.3333333 0.3933934 0.4204204
#   
#   $`jedit-5.1.0`$indices
#   [1] 3
#   
#   
#   $`hadoop-0.20.2`
#   $`hadoop-0.20.2`$mojosim
#   [1] 0.2872340 0.2588652 0.2836879 0.2801418 0.2801418 0.2907801 0.2907801 0.2801418 0.2872340 0.2801418 0.2907801
#   
#   $`hadoop-0.20.2`$indices
#   [1]  6  7 11
# 
# $`eclipse-jdt-core-3.8`
# $`eclipse-jdt-core-3.8`$mojosim
# [1] 0.4109589 0.4452055 0.5000000 0.4794521 0.4178082 0.4041096 0.4520548 0.4863014 0.5273973 0.4383562 0.5068493
# 
# $`eclipse-jdt-core-3.8`$indices
# [1] 9
# 
# 
# $`weka-3.6.11`
# $`weka-3.6.11`$mojosim
# [1] 0.1696429 0.1763393 0.3638393 0.3348214 0.3013393 0.1607143 0.2991071 0.1629464 0.2991071 0.1629464 0.1629464
# 
# $`weka-3.6.11`$indices
# [1] 3
  
}

print_latex <- function(result, prname) {
  
  all_values <- unlist(lapply(result, function(r) r$value))
  max_value <- max(all_values)
  
  str <- prname
  
  for(i in 1:length(result)) {
    
    merged_pars = result[[i]]$par[1]
    
    if (length(result[[i]]$par) > 1)
      for (j in 2:length(result[[i]]$par))
        merged_pars <- paste(merged_pars, result[[i]]$par[j], sep=",")
    
    
#     val <- paste(round(result[[i]]$value, 3), "(", round(result[[i]]$variance, 4) , ")", sep="")
    val <- round(result[[i]]$value, 3)
    
    if (result[[i]]$value == max_value) {
      merged_pars <- paste("\\cellcolor[gray]{0.8}{", merged_pars, "}")
      val <- paste("\\cellcolor[gray]{0.8}{", val, "}")
      
    }
    
    if (length(result[[i]]$par) > 0)
      l <- paste(merged_pars, val , sep="&")
    else
      l <- val
    
    str <- paste(str, l, sep="&")
    
  }
  
  paste(str, "\\\\")
  
  
}