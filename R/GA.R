

#Input: mydata
# eval.funcs: a list of functions
two.archive.search <- function(mydata, names, eval.funclist, popsize, 
                               iters, mutchance, crossoverprob, cdprob, nochanges, limit ){
  
  if (length(mydata) != length(eval.funclist))
    stop("The size of data and functions do not match")
  
  dimensions <- sapply(mydata, function(d){dim(d)})
  
  #is it one or two
  if (!all(apply(dimensions, 1, function(r) unique(r))))
    stop("the matrices have different dimensions")
  
  if (!all(apply(dimensions, 2, function(c) unique(c))))
    stop("the matrices are not square")
  
  #a trick
  len <- dimensions[1,1]
  if (len == 0) 
    stop("Empty dataset")
  
  nobjs <- length(mydata)
  
  archives = list(converging = list(),  diversity = list())

  make_population <- function(k) {
    partition = list(
      group = NULL,
      result = NULL
    )

    partition$group <- gelato::normalizeVector(sample(c(1:k), len, replace = TRUE))
    
    names(partition$group) = names

    partition$result <- lapply(eval.funclist, function(eval) eval(mydata, partition$group))
    
    return(partition) 
  }
  
  
  noc <- sample(c(1:len), popsize, replace=TRUE)
  
  #populate the population
  population <- lapply(noc, make_population)  
  
  
#Precondition: 0 <= mutation chance <= 1 
mutate <- function(group, mutationchance, len) {
  
  
  
}  
  
  
#Input: 
mutate<- function(group, mutationchance, len) {
  
  k <- max(group)
  
  nmutations = ceiling(mutationchance * len)
  
  positions = sample(1:len, nmutations, replace=FALSE)
  values = sample(1:(k+1), nmutations, replace=TRUE)
  
  group[positions] = values
  
  return(gelato::normalizeVector(group))
   
}

cross.over <- function(parents, crossoverpoint, len) {
  
  children=list()
  
  children[[1]] <- c(head(parents[[1]], crossoverpoint), 
                     tail(parents[[2]], len - crossoverpoint))
  children[[2]] <- c(head(parents[[2]], crossoverpoint), 
                     tail(parents[[1]], len - crossoverpoint))

  return(children)
}


find.nondominated <- function(pop) {

  results <- sapply(pop, function(x) x$result)
  
  for (i in length(pop):1)
    for (j in 1:length(pop)) {
      if (i != j)
        if (dominates(pop[[j]]$result, pop[[i]]$result)) {
          pop <- pop[-i]
          break
        }
    }  
  return(pop)
}

#Input: The sizes match
dominates <- function(r1, r2) {

  if (all(r1%in%r2))
    return(FALSE)
  
  for (i in 1:length(r1))
    if (r1[[i]] < r2[[i]])
      return(FALSE)
  
  return(TRUE)
}


check.whether.dominates <- function(r1s, r2) {
  
  if (length(r1s) == 0)
    return(FALSE)

  for (i in 1:length(r1s)) 
    if (dominates(r1s[[i]], r2))
      return(TRUE)

  
  return(FALSE)
} 

is.dominated <- function(r1, r2) {
  
 if (dominates(r2, r1))
   return(TRUE)
 
  return(FALSE)
}



find.unique <- function(archive) {
  
  groups <- lapply(archive, function(x) x$group)
  
  return(archive[which(!duplicated(groups))])
  
}


collect_nondominated <- function(pop, archives) {
  
  collection = list(archives=NULL, changed=FALSE)
  
  nondominated.pop <- find.nondominated(pop)
  
  for(i in 1:length(nondominated.pop)) {
    
    #if no member in both archives can dominate the indivisual(i)
    
    ca.results <- lapply(archives$converging, function(x) x$result)
    
    da.results <- lapply(archives$diversity, function(x) x$result)
    
    result <- nondominated.pop[[i]]$result
    
    if(!check.whether.dominates(ca.results, result))
      if (!check.whether.dominates(da.results, result)){        
        collection$changed=TRUE
        
      #if individual(i) dominates any member in both archives then 
      dflag = FALSE
      
      ca.dominated <- sapply(ca.results, is.dominated, result)
      da.dominated <- sapply(da.results, is.dominated, result)
      
      if (any(ca.dominated)) {
        dflag = TRUE
        
        #Delete elements
        archives$converging <- archives$converging[-which(ca.dominated)]
        
      }
      
      if (any(da.dominated)) {
        dflag = TRUE
        
        #Delete elements
        archives$diversity <- archives$diversity[-which(da.dominated)]
      }
      
      if (dflag)
        archives$converging <- append(archives$converging, nondominated.pop[i])
      else 
        archives$diversity <- append(archives$diversity, nondominated.pop[i])

      
    }
  }
  
  archives$converging <- find.unique(archives$converging)
  archives$diversity <- find.unique(archives$diversity)
  
  
  
  
  # Removing strategy to shrink the size of the archives below the limit
  #FIXME based on the euclidean distance
  while((length(archives$converging) + length(archives$diversity) > limit)
        && (length(archives$diversity) != 0)) {
    
    index <- sample(1:length(archives$diversity), 1)    
    archives$diversity <-  archives$diversity[-index]
    
  }

  collection$archives = archives
  
  return(collection)
  
}

#Converging.archive must be nonempty
select_parents <- function(cdprob, nparents =2) {
  parents = list()
  np = 1
  
  while(np <= nparents) {
    
    cd <- sample(1:2, 1, prob= c(cdprob, 1-cdprob))

    if ((cd==1) && (length(archives$diversity)>0)) {
      index <- sample(1:length(archives$diversity), 1)
      
      parents[[np]] <- archives$diversity[[index]]$group
      
      np <- np+1
      
    }else  #The converging archive can be empty
      if (length(archives$converging)>0){ 
        index <- sample(1:length(archives$converging), 1)
        
        parents[[np]] <- archives$converging[[index]]$group        
        
        np <- np+1
        
      }
  }
  
  return(parents)
  
}

iter = 1

nochanges.seen = 0

while(iter <= iters) {
            
  collection <- collect_nondominated(population, archives)
  
  archives = collection$archives
  
  if (!collection$changed)
    nochanges.seen <- nochanges.seen +1
  else
    nochanges.seen = 0

  print(nochanges.seen)
  
  if (nochanges.seen > nochanges)
    return(archives)

  #print(system.time(archives <- collectNondominated(population, archives)))

  print(length(archives$converging))
  print(length(archives$diversity))
  
  newpopulation = list(group=NULL, result=NULL)
  
  newpopsize = 1
  
  while (newpopsize <= popsize) {
    
    parents <- select_parents(cdprob)
  

    crossover = sample(c(0,1), 1, prob=c(1-crossoverprob, crossoverprob))
    
    if (crossover==1) {
      #one-point cross-over
      crossoverpoint <- sample(1:len, 1)
      children <- cross.over(parents, crossoverpoint, len) 
    }
    else
        children = parents
    
    
    
    children <- lapply(children, mutate, mutchance, len)
    
    for (i in 1:length(children))  
    {
      result <- lapply(eval.funclist, function(eval) eval(mydata, children[[i]]))
      
      #populate new population
      newpopulation[[newpopsize]] <- list(group=children[[i]], result = result)      
      
      newpopsize <- newpopsize +1
    } 
    
  }
  
  population = newpopulation
  
  iter <- iter+1
}

return(archives)
}

dist_fun <- function(x){
  require(gelato)
  print(x[1])
  apply.two.archive(x[1], x[2])
  
}


run.GA.experiments <- function() {

prnames = c("apache-log4j-1.2.17", "eclipse-jdt-core-3.8", "jedit-5.1.0", "jfreechart-1.2.0")
iters = c(200, 50, 50, 100)

df = data.frame(prnames, iters)

sfInit(parallel=TRUE, cpus=4)

sfExport("df", "dist_fun", "remove.diagonal", "two.archive.search", "apply.two.archive", "Cpp.MQ.evaluator", "Cpp.LQ.evaluator", "import.bunch.matrix", "exclude.scope.ext" , "extensions", "remove.documents", "compute.rank", "make.compatible")

sfClusterApplyLB(df, dist_fun)
sfStop()

}

apply.two.archive <- function(prname, iters.coefficient) {

  setwd("~/workspace")
  
  require(igraph)
  require(Rcpp)
  require(lsa)
  
  #Load the priori decomposition
  decomposition <- read.csv(paste("benchmark", prname ,"decomposition.csv", sep="/"), sep=",",  header = TRUE)
  priori.decomp <- decomposition$x
  names(priori.decomp) <- decomposition$X
  priori.decomp <- gelato::normalizeVector(priori.decomp)
  
  #Load the adjacency matrix
  extensions= c("java/", "org/xml/", "javax/")
  cfg <- import.bunch.matrix(paste("benchmark", prname ,"dep_graph.txt", sep="/"), exclude.ext=extensions)
  #cfg <- read.table("benchmark/jedit-5.1.0/cfg.csv", sep=",", row.names = 1, header = TRUE, check.names = FALSE)
  cfg <- cfg[intersect(rownames(cfg),names(priori.decomp)), intersect(colnames(cfg),names(priori.decomp))]
  
  
  #Load and compute the lexical sim matrix
  bow <- read.table(paste("benchmark", prname , "mydata-idf-BoW-matrix.csv", sep="/"), sep=",", row.names = 1, header = TRUE, check.names = FALSE)
  bow <- remove.documents(bow, rownames(cfg))
  #write.table(bow, file =paste("benchmark", prname , "mydata-compatible-idf-BoW-matrix.csv", sep="/"),row.names=TRUE, col.names=NA,sep=",", quote=FALSE)
  
  #Apply LSA
  ndims = compute.rank(dim(bow)[1], dim(bow)[2])
  print("No of Dimensions:")
  print(ndims)
  
  space1 = lsa(bow, ndims)
  lsa.bow = as.textmatrix(space1)
  lexsim = cosine(t(lsa.bow))
  
  
  #Build mydata
  mydata = list(cfg=as.matrix(cfg), lexsim=as.matrix(lexsim))
  mydata <- make.compatible(mydata)
  mydata$cfg <- remove.diagonal(mydata$cfg)
  mydata$lexsim <- remove.diagonal(mydata$lexsim)
   
  N <- dim(mydata$cfg)[1]
  
  popsize <- 10 * N
  
  if (N >= 100)
    coprob = 1
  else  
    coprob = 0.8
    
  
  iters <- iters.coefficient * N
  
  mutchance <- 0.004 * log2(N)
  
  nochanges = 300

  cdprob <- 0.8
  
  limit = 3 * popsize
  
  archives <- two.archive.search(mydata, rownames(mydata$cfg), list(mq=Cpp.MQ.evaluator, lq=Cpp.LQ.evaluator), 
                     popsize, iters, mutchance, coprob, cdprob, nochanges, limit)

  converging_groups<- lapply(archives$converging, function(x) x$group)

  diversity_groups <- lapply(archives$diversity, function(x) x$group)
  
  archive_groups <- append(converging_groups, diversity_groups)
    
  
  #Compute distance  
  priori.decomp <- find.intersection(priori.decomp, archive_groups[[1]])

  priori.decomp <- gelato::normalizeVector(priori.decomp)
  
  noc <- lapply(archive_groups, function(g) max(g) )
  
  mojo <- lapply(archive_groups, function(c) compute.MoJo(c, priori.decomp))
  
  mojosim <- lapply(mojo, function(m) 1 - (m/N))
  
  #Save Results
  results = matrix(nrow = length(archive_groups), ncol = 3)
  
  for (i in 1:length(archive_groups)) {
    
    results[i,1] <- noc[[i]]
    results[i,2] <- mojo[[i]]
    results[i,3] <- mojosim[[i]]
    
  }
  
  colnames(results) <- c("NOC", "MoJo", "MoJoSim")
  
  write.table(results, file =paste("benchmark", prname ,"multiobj-results.csv", sep="/"), row.names=TRUE, col.names=NA,sep=",", quote=FALSE)
  
  #Store groups
  groups = matrix(nrow = 0, ncol = dim(mydata$cfg)[1])
  for (i in 1:length(archive_groups))
    groups <- rbind(groups, archive_groups[[i]])
  
  colnames(groups) <- rownames(mydata$cfg)
  #print(groups)
  write.table(groups, file =paste("benchmark", prname ,"multiobj-groups.csv", sep="/"),row.names=TRUE, col.names=NA,sep=",", quote=FALSE)
  
  
  converging_results <- lapply(archives$converging, function(x) x$result)
  
  diversity_results <- lapply(archives$diversity, function(x) x$result)
  
  archive_results <- append(converging_results, diversity_results)
  
  
  #Save Objective Scores
  scores = matrix(nrow = length(archive_results), ncol = 2)
  
  for (i in 1:length(archive_results)) {
    
    scores[i,1] <- archive_results[[i]][[1]]
    scores[i,2] <- archive_results[[i]][[2]]
    
  }
  
  colnames(scores) <- c("MQ", "LQ")
  
  print(scores)
  write.table(scores, file =paste("benchmark", prname ,"multiobj-scores.csv", sep="/"),row.names=TRUE, col.names=NA,sep=",", quote=FALSE)
  
}

analyze.multiobj.results <- function(prname) {
  
  require(scatterplot3d)
  require(ggplot2)
  
  #setwd("~/workspace")
  setwd("~/workspace")
    
  #Load the results
  results <- read.csv(paste("benchmark", prname ,"multiobj-results.csv", sep="/"), sep=",",  header = TRUE)

  results <- results[,-1]
  
  scores <- read.csv(paste("benchmark", prname ,"multiobj-scores.csv", sep="/"), sep=",",  header = TRUE)
  
  scores <- scores[,-1]
  
  if(dim(results)[1] != dim(scores)[1])
    stop("incompatible dataset")
  
  results <- cbind(results, scores)

  
  print(dim(results))
  results <- unique( results[ , 1:5 ] )
  
  
  print(dim(results))
  
  print(results[1,5])
  
  
  print("size of the pareto front:")
  print(dim(results)[1])
  
  best_case <- find.element(results, max)
  worst_case <- find.element(results, min)
  median_case <- find.element(results, median)
  
  average_case <- avg.elements(results)
  
  best_output <- paste(best_case$noc, round(best_case$mojosim, digits = 2), sep = "&")
  worst_output <- paste(worst_case$noc, round(worst_case$mojosim, digits = 2), sep = "&")
  median_output <- paste(median_case$noc, round(median_case$mojosim, digits = 2), sep = "&")
  average_output <- paste(round(average_case$noc, digits = 2), round(average_case$mojosim, digits = 2), sep = "&")
  
  
  output <- paste(dim(results)[1], best_output, worst_output, median_output, average_output, sep="&")
  print(output)
  
  #png(paste("benchmark", prname ,"multiobj_m_scatterplot.png", sep="/"), pointsize = 15, width = 1000, height = 1000)
 # with(results, {
 #attach(results) 
 #s3d <- scatterplot3d(MQ, LQ, MoJoSim,        # x y and z axis
 #                        color="blue", pch=15,        # filled blue circles
 #                        type="h",                    # vertical lines to the x-y plane
  #                       main=paste("3-D Scatterplot for", prname, sep = " "),
  #                       xlab="MQ",
  #                       ylab="LQ",
  #                       zlab="MoJoSim")
  #s3d.coords <- s3d$xyz.convert(MQ, LQ, MoJoSim) # convert 3D coords to 2D projection
 #fit <- lm(MoJoSim ~ MQ+LQ)
 #s3d$plane3d(fit)
 
  #  text(s3d.coords$x, s3d.coords$y,             # x and y coordinates
        # labels=row.names(mtcars),               # text to plot
  #       cex=.5, pos=4)           # shrink text 50% and place to right of points)
 # })
 
 results <- data.frame(results)
 
 colnames(results)[5] = "CQ"
 
 mjsim = results$MoJoSim
 
 best.MoJoSim = max(mjsim)
 snd.best.MoJoSim = max(mjsim[mjsim != best.MoJoSim])
 
 print(snd.best.MoJoSim)
 best <- subset(results, MoJoSim == best.MoJoSim)
 snd.best <- subset(results, MoJoSim == snd.best.MoJoSim)
 
 print(colnames(results))
 
 #format mojosim
#   results$MoJoSim <- round(results$MoJoSim, 2)
 
 print(min(results$MoJoSim))
 print(max(results$MoJoSim))
 
 t1<-theme(                              
   plot.background = element_blank(), 
   panel.grid.major = element_blank(), 
   panel.grid.minor = element_blank(), 
   panel.border = element_blank(), 
   panel.background = element_blank(),
   axis.line = element_line(size=.4),
   axis.text=element_text(size=12),
   axis.title=element_text(size=14,face="bold")
#    legend.title=element_text(size=10, vjust=-12)
#     guide_colourbar.title = element_text(draw.ulim = FALSE, draw.llim = FALSE)
#    legend.key.height=unit(3,"line"),
#    legend.key.width=unit(3,"line")
 )
 
  so.scores <- read.csv(paste("benchmark", prname ,"singleobj-scores.csv", sep="/"), sep=",",  header = TRUE)
  so.scores <- so.scores[,-1]

  so.scores <- data.frame(so.scores)
  colnames(so.scores)[2]<- "CQ"


 g1 <- ggplot(results, aes(x=MQ, y=CQ)) + t1 + geom_point(size=6, aes(color=MoJoSim))
 g1 <- g1 + geom_point(data=best, aes(x=MQ, y=CQ), shape=22, fill="red", size=10)


# g1 <- g1 + geom_smooth(se=FALSE)
 g1 <- g1 + scale_color_continuous(breaks = with(results, 
            round(seq(min(MoJoSim), max(MoJoSim),length.out=5), 3)), 
             guide = guide_colourbar(label.theme = element_text(size=12, angle=0),
                                    title.theme = element_text(size=14, face="bold", angle=0),
                                     draw.ulim = F,
                                     draw.llim = F,
                                    nbin=100),
              low = "black", high = "lightgray")  

  
# g1 <- g1 + geom_point(data=so.scores, aes(x=MQ, y=CQ), shape=24, fill="blue", size=10)
  
 ggsave(filename= paste("benchmark", prname ,"multiobj_mv_scatterplot.png", sep="/"), plot=g1, pointsize = 15, width = 10, height = 10)
 
 print(g1)
 
dev.off() 
 # return(results)  



}


multi.plot <- function(prnames, cols=c("red", "green", "blue")) {
  require(ggplot2)
  require(splines)
  
  t1<-theme(                              
    plot.background = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    panel.background = element_blank(),
    axis.line = element_line(size=.4),
    axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold")
    #     guide_colourbar.title = element_text(draw.ulim = FALSE, draw.llim = FALSE)
    #    legend.key.height=unit(3,"line"),
    #    legend.key.width=unit(3,"line")
  )
  

  
  #setwd("~/workspace")
  setwd("~/workspace/sci.hage0101.GELATO")
  
  #Load the results
  results <- lapply(prnames, function(prname) read.csv(paste("benchmark", prname ,"multiobj-results.csv", sep="/"), sep=",",  header = TRUE))
  
  results <- lapply(results, function(x) x[,-1])
  
  scores <- lapply(prnames, function(prname) read.csv(paste("benchmark", prname ,"multiobj-scores.csv", sep="/"), sep=",",  header = TRUE))
  
  scores <- lapply(scores, function(x) x[,-1])
  
  #if(dim(results)[1] != dim(scores)[1])
  #  stop("incompatible dataset")
  
  for(i in 1:length(results))
    results[[i]] <- cbind(results[[i]], scores[[i]])
  
  results <- lapply(results, function(r) unique( r[ , 1:5 ] ))
  
  
  d <- lapply(results, function(r) dist(r[,4:5])) # euclidean distances between the COLS
  fits <- lapply(d, function(x) cmdscale(x,eig=TRUE, k=1)) # k is the number of dim
  xs <- lapply(fits, function(fit) fit$points[,1])
  
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  
  xs <- lapply(xs, function(x) range01(x))
  
  MoJoSims = lapply(results, function(r) r[,3])
  
  dfs = list()
  for (i in 1:length(prnames)) {
    x = xs[[i]]
    MoJoSim = MoJoSims[[i]]
   
    dfs[[i]] = data.frame(x, MoJoSim)
  }
  

  plot <- ggplot(NULL, aes(x, MoJoSim)) + t1 + ggtitle("Plot") + t1
  
  for (i in 1:length(dfs)) {
    
    plot <- plot + 
      stat_smooth(data = dfs[[i]], method = "lm", se= FALSE, formula = y ~ ns(x,10), colour = cols[i])
    
  }
     print(plot)
  
}

avg.elements <- function(mydata) {
  element = list(noc= 0, mojosim = 0)

  element$mojosim = mean(mydata$MoJoSim)
  element$noc = mean(mydata$NOC)
  
  return(element)
  
}

find.element <- function(mydata, fun) {
  element = list(noc= 0, mojosim = 0)
  
  element$mojosim = fun(mydata$MoJoSim)
  
  for (i in 1:dim(mydata)[1])
    if(mydata[i, ]$MoJoSim == element$mojosim) {
      element$noc <- mydata[i, ]$NOC
      break()
    }
  
  return(element)
  
}
