#Lexsim test


find.best.lexsim.kernel <- function(prname, par, strs, decays, rootFolder="org") {
  
  require(lsa)
  require(proxy)
  
  setwd("~/workspace")
  
  #Read the set of classnames for running the experiment
  classnames <- unlist(read.table(paste("benchmark", prname , paste("MULTIVIEW", "classnames.txt" ,sep="/") , sep="/")) )
  
  #Load the priori decomposition
  decomposition <- read.csv(paste("benchmark", prname ,"decomposition.csv", sep="/"), sep=",",  header = TRUE)
  priori.decomp <- decomposition$x
  names(priori.decomp) <- decomposition$X
  priori.decomp <- priori.decomp[which(names(priori.decomp) %in% classnames)] # find the ones in classnames
  priori.decomp <- priori.decomp[order(names(priori.decomp))]
  priori.decomp <- normalizeVector(priori.decomp)
  
  
  #Initialize the text of source files
  
  setwd(paste("benchmark", prname, sep="/"))
  txts <- read.text.directory(rootFolder, pattern = "*.java")
  
  #Fix bow and txts
  bow <- bow[intersect(rownames(bow),rownames(cfg)),]
  bow <- bow[intersect(rownames(bow),names(txts)),]
  bow <- bow[intersect(rownames(bow),names(priori.decomp)),]
  
  bow <- bow[order(rownames(bow)),]
  
  txts <- txts[rownames(bow)]
  txts <- txts[order(names(txts))]
  
  
  #bow: a document-term matrix: documents in columns and terms in rows  
  bow <- t(bow)
  
  bow <- lw_logtf(bow) * gw_idf(bow)
  
  bow <- bow[complete.cases(bow),]
  
  ndims = compute.rank(dim(bow)[2], dim(bow)[1])
  print("No of Dimensions:")
  print(ndims)
  
  space1 = lsa(bow, ndims)
  #   lsa.bow = space1$tk %*% diag(space1$sk) %*% t(space1$dk)
  
  lsa.bow = diag(sqrt(space1$sk)) %*%  t(space1$dk)
  
  #   lsa.bow = t(space1$dk)
  
  #make lsa bow, to documents in rows and terms in columns
  lsa.bow <-  t(lsa.bow)
  
  #Check to ensure  lsa.bow is a documentTerm matrix
  #   lexsim.kernel = gaussian.kernel(lsa.bow, 0.00000001) #ndims/N
  
  #   lexsim.kernel = cosine(t(lsa.bow))
  
  #Fix priori decomposition
  dummy_v <- rep(0, dim(lsa.bow)[1])
  names(dummy_v) <- rownames(lsa.bow)
  
  priori.decomp <- find.intersection(priori.decomp, dummy_v)
  priori.decomp <- normalizeVector(priori.decomp)
  priori.decomp <- priori.decomp[order(names(priori.decomp))]
  
  for (i in 1:length(names(dummy_v)))
    if (!(names(dummy_v)[i] == names(priori.decomp)[i]))
      stop('shsdf')  
  
  priori.decomp <- normalizeVector(priori.decomp)
  
  print("printing the numer of groups:")
  print(max(priori.decomp))
  
  k <- max(priori.decomp)
  
  lexsim.kernel <- list()
  
  lexsim.kernel[[length(lexsim.kernel) + 1]] <- as.matrix(dist(lsa.bow, method="cosine"))
  #   lexsim.kernel[[length(lexsim.kernel) + 1]] <- as.matrix(dist(lsa.bow, method="Euclidean"))
  
  degrees <- 1:7
  for(i in 1:length(degrees))
    lexsim.kernel[[length(lexsim.kernel) + 1]] <- polynomial.kernel(lsa.bow, degrees[i])
  
  
  for(i in 1:length(par))
    lexsim.kernel[[length(lexsim.kernel) + 1]] <- gaussian.kernel(lsa.bow, par[i])
  
  #compute the euclidean distance, and take the median as the parameter
  
  d <- dist(lsa.bow, method = "Euclidean", diag = FALSE, upper = FALSE)
  
  m <- median(d)
  
  print("median")
  print(m)
  
  lexsim.kernel[[length(lexsim.kernel) + 1]] <- gaussian.kernel(lsa.bow, m)
  lexsim.kernel[[length(lexsim.kernel) + 1]] <- gaussian.kernel(lsa.bow, ndims)
  
  
  
  #String kernels
  
  #constant kernel
  lexsim.kernel[[length(lexsim.kernel) + 1]] <- constant.string.kernel(txts) 

  #p-spectrum Kernel
  for(i in 1:length(strs))
        lexsim.kernel[[length(lexsim.kernel) + 1]] <- spectrum.string.kernel(txts, strs[i])

  #exponential kernel
  for(i in 1:length(decays))
      lexsim.kernel[[length(lexsim.kernel) + 1]] <- exponential.decay.kernel(txts, decays[i]) 
  
  #   return(lexsim.kernel[[length(lexsim.kernel)]])
  
  laps <- lapply(lexsim.kernel, function(kern) laplacian(kern, TRUE))
  
  print(length(laps))
  results=list()
  for (i in 1:length(laps)){
    print(i)
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
  print(mojoSims)
  
  #Separating the results
  
  best = list()
  
  #initialize current position of the window
  cur_window = 0
  
  #cosine window
  window <- 1
  
  best.mjsim <- max(mojoSims[(cur_window + 1):(cur_window + window)])
  
  index <- which(mojoSims[(cur_window + 1):(cur_window + window)]==best.mjsim)
  
  variance <- var(mojoSims[(cur_window + 1):(cur_window + window)])
  
  best[[length(best) + 1]] <- list(par=NULL, value = best.mjsim, variance = variance)
  
  #set the new current window
  cur_window <- cur_window + window
  
  #For polynomials
  window <- length(degrees)

  best.mjsim <- max(mojoSims[(cur_window + 1):(cur_window + window)])
  
  index <- which(mojoSims[(cur_window + 1):(cur_window + window)]==best.mjsim)
  
  variance <- var(mojoSims[(cur_window + 1):(cur_window + window)])
  
  best[[length(best) + 1]] <- list(par=degrees[index], value = best.mjsim, variance = variance)
  
  #set the new current window
  cur_window <- cur_window + window
  
  #For gaussian
  window <- length(par) + 2
  
  #Gaussian
  best.mjsim <- max(mojoSims[(cur_window + 1):(cur_window + window)])
  
  index <- which(mojoSims[(cur_window + 1):(cur_window + window)]==best.mjsim)
  
  variance <- var(mojoSims[(cur_window + 1):(cur_window + window)])
  
  if (index == (length(par) + 1))
    best[[length(best) + 1]] <- list(par=m, value = best.mjsim, variance = variance)
  else if (index == (length(par) + 2))
    best[[length(best) + 1]] <- list(par=ndims, value = best.mjsim, variance = variance)
  else
    best[[length(best) + 1]] <- list(par=par[index], value = best.mjsim, variance = variance)
  
  #set the new current window
  cur_window <- cur_window + window
  
  #For constant
  window <- 1
  
  #Constant
  best.mjsim <- max(mojoSims[(cur_window + 1):(cur_window + window)])
  
  index <- which(mojoSims[(cur_window + 1):(cur_window + window)]==best.mjsim)
  
  variance <- var(mojoSims[(cur_window + 1):(cur_window + window)])
  
  best[[length(best) + 1]] <- list(par=NULL, value = best.mjsim, variance = variance)
  
  #set the new current window
  cur_window <- cur_window + window
  
  #For p-spectrum
  window <- length(strs)
  
  best.mjsim <- max(mojoSims[(cur_window + 1):(cur_window + window)])
  
  index <- which(mojoSims[(cur_window + 1):(cur_window + window)]==best.mjsim)
  
  variance <- var(mojoSims[(cur_window + 1):(cur_window + window)])
  
  best[[length(best) + 1]] <- list(par=strs[index], value = best.mjsim, variance = variance)
  
  #set the new current window
  cur_window <- cur_window + window
  
  #For exponential
  window <- length(decays)
  
  best.mjsim <- max(mojoSims[(cur_window + 1):(cur_window + window)])
  
  index <- which(mojoSims[(cur_window + 1):(cur_window + window)]==best.mjsim)
  
  variance <- var(mojoSims[(cur_window + 1):(cur_window + window)])
  
  best[[length(best) + 1]] <- list(par=decays[index], value = best.mjsim, variance = variance)
  
  
  return(best)
  
  best.mjsim = max(mojoSims)
  #   return(mojoSims)
  #   return(which(mojoSims==best.mjsim))
  return(list(mojosim = mojoSims, indices= which(mojoSims==best.mjsim)))
  
}

test.best.lexsim.kernel <- function(indices = 0) {
#   projects <- list("apache-ant-1.9.3", "hadoop-0.20.2", "apache-log4j-1.2.17", "eclipse-jdt-core-3.8",
#                    "jdom-2.0.5", "jedit-5.1.0", "jfreechart-1.2.0", "jhotdraw-7.0.6", "junit-4.12" ,"weka-3.6.11")
#   
#   names(projects) <- c("Apache Ant", "Apache Hadoop", "Apache Log4j", "Eclipse JDT Core", "JDOM", 
#                        "JEdit", "JFreeChart", "JHotDraw", "JUnit", "Weka")
  
#   projects <- list("hadoop-0.20.2", "apache-log4j-1.2.17", "jdom-2.0.5")#,
#   projects <- list("jedit-5.1.0", "jfreechart-1.2.0", "jhotdraw-7.0.6")
#   projects <- list("weka-3.6.11")
#   projects <- list("apache-ant-1.9.3")
#   projects <- list("junit-4.12")
  projects <- list("eclipse-jdt-core-3.8")
  
#   names(projects) <- c("Apache Hadoop", "Apache Log4j", "JDOM")#, 
#   names(projects) <-  c("JEdit", "JFreeChart", "JHotDraw")
#   names(projects) <-  c("Weka")
# names(projects) <-  c("Apache Ant")
# names(projects) <-  c("JUnit")
names(projects) <-  c("Eclipse JDT Core")

  
  #r <-find.best.kernel(prname, c(10^(-6), 10^(-5), 10^(-4), 10^(-3), 10^(-2), 10^(-1), 10^(0), 5, 10^(1), 5*10^(1), 10^(2), 10^(3) ))
  
  pars <- c(10^(-2), 10^(-1), 10^(0), 5, 10^(1), 10^(2) )
  
  lens <- c(2,4,6, 8, 10,12, 15, 20)

# lens <- c(2,4)
  
  #   decays <- c(0.1, 1, 2, 5, 10)
  
  decays <- 1.1^seq(6)
# decays <- 1.1^seq(2)
  
  results <- lapply(projects, function(p) find.best.lexsim.kernel(p, pars, lens, decays))
  
  # plot(c(min(decay),max(decay)),c(0,1),type='n',xlab="lambda",ylab='MoJoSim')
  
  
  names(results) <- names(projects)
  
  # plot(c(min(decay),max(decay)),c(0,1),type='n',xlab="lambda",ylab='MoJoSim')
  
  
  
  #\cellcolor[gray]{0.8}
  
  if (indices == 0)
    indices <- 1:length(projects)
  
  str2tex = ""
  
  for (i in 1:length(results))
    str2tex <- paste(str2tex, print_latex(results[[i]], names(results[i])), "\n" )
  
  str2tex
  
}

get.best.mojosims <- function(results) {
  mojosims <- lapply(results, function(r) r$mojosim)
  best.mojosims <- lapply(mojosims, function(m) max(m))
  unlist(best.mojosims)
  
}

prepare.string.kernel <- function(prname, dirname="org", pattern = "*.java") {
  require("tm")
  require("kernlab")
  
  setwd(paste("~/workspace/benchmark/", prname, sep=""))
  corpus <- Corpus(DirSource(dirname, recursive = T, pattern = pattern ))
  
  (f <- content_transformer(function(x, pattern) gsub(pattern, "", x)))
  corpus <- tm_map(corpus, f, "[\t\n]+")
  corpus <- tm_map(corpus, stemDocument)
  
  corpus <- tm_map(corpus, removeWords, stopwords("english"))
  
  corpus
}