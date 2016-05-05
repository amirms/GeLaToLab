

find.best.graph.kernel <- function(prname, par) {
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
  
  
  #Load the adjacency matrix
  extensions= c("java/", "org/xml/", "javax/")
  cfg <- import.bunch.matrix(paste("benchmark", prname ,"dep_graph.txt", sep="/"), exclude.ext=extensions)
  #   cfg <- unweight.adjacency(cfg)
  cfg <- cfg[which(rownames(cfg) %in% classnames), which(colnames(cfg) %in% classnames)]
  cfg <- cfg[order(rownames(cfg)), order(colnames(cfg))]
  
  
  #   return(cfg)
  if (!all(rownames(cfg) == names(priori.decomp)))
    stop('names do not match!')

  
  print("printing the numer of groups:")
  print(max(priori.decomp))
  
  k <- max(priori.decomp)
  
  # no self-references
  diag(cfg) <- 0 
  cfg <- make.symmetric(cfg)
  d <- apply(abs(cfg),1,sum)
  
  indices <- which(d==0)
  
  while(length(indices) > 0) {
    cfg <- cfg[-indices, -indices]
    d <- apply(abs(cfg),1,sum)
    indices <- which(d==0)
    
  }
  
  cfg <- make.symmetric(cfg)
  d <- apply(abs(cfg),1,sum)
  D <- diag(d)
  
  
  cfg.kernel <- list()
  
  #Calculate exponential diffusion kernel
  for (i in 1:(length(par)-4))
  cfg.kernel[[length(cfg.kernel) + 1]] <- compute.exponential.diffusion.kernel(cfg, alpha = par[i])

  
  #Compute the normalized graph laplacian of symmetric CFG
  cfg.laplacian = laplacian(cfg, TRUE)
  
  #Calculate the laplacian exponential diffusion kernel
  for (i in 1:length(par))
    cfg.kernel[[length(cfg.kernel) + 1]] <- calc.diffusion.kernel(cfg.laplacian, beta = par[i])
  
  #Calculate the commute-time kernel
  cfg.kernel[[length(cfg.kernel) + 1]] <- compute.avg.commute.time.kernel(cfg, D)
  
  norm_vec <- function(x) sqrt(sum(x^2))
  #Calculate the von neumann diffusion kernel
  params <- c(10^(-5), 10^(-4), 10^(-3), 10^(-2), 10^(-1), 0.99) * norm_vec(cfg)^-1
  for (i in 1:length(params))
    cfg.kernel[[length(cfg.kernel) + 1]] <- compute.von.neumann.diffusion.kernel(cfg, params[i])
  
  #Calculate the sigmoid commute-time kernel
#   sig <- c(1,2,4,5,10)
#   for (i in 1:length(sig))
#     cfg.kernel[[length(cfg.kernel) + 1]] <- compute.sigmoid.commute.time.kernel(cfg, D, sig[i])
  
#   plot.influence(cfg.kernel[[length(par) - 3 + 5]], cfg.kernel[[length(cfg.kernel)]])
  
  
  #   for (i in 1:length(cfg.diffusion.kernel2)){
  #     dimnames(cfg.diffusion.kernel2[[i]]) <- dimnames(cfg)
  # #     cfg.diffusion.kernel2[[i]] <- normalize.influence(cfg.diffusion.kernel2[[i]])
  #     cfg.kernel[[length(cfg.kernel) + 1]] <- cfg.diffusion.kernel2[[i]]
  #     
  #   }
  
  #   return(cfg.diffusion.kernel[[1]])  
  laps <- lapply(cfg.kernel, function(kern) laplacian(kern, TRUE))                     
  
  for (i in 1:length(laps)){
    print(i)
    results[[i]] <- spectral.clustering(laps[[i]], k)
    
  }
  
  for (i in 1:length(results)) {
    names(results[[i]]) <- rownames(laps[[i]])
    results[[i]] <- normalizeVector(results[[i]])
    
  }
  
  N <- length(priori.decomp)
  
  mojoSims <- unlist(lapply(results, function(r) 1 - (compute.MoJo(r, priori.decomp)/N)))
  
  print(mojoSims)
  
  #Separating the results  
  best = list()
  
  #initialize current position of the window
  cur_window = 0
  
  #For Exponential Diffusion
  window <- length(par) - 4
  
  best.mjsim <- max(mojoSims[(cur_window + 1):(cur_window + window)])
  
  index <- which(mojoSims[(cur_window + 1):(cur_window + window)]==best.mjsim)
  
  variance <- var(mojoSims[(cur_window + 1):(cur_window + window)])
  
  best[[length(best) + 1]] <- list(par=par[index], value = best.mjsim, variance = variance)
  
  #set the new current window
  cur_window <- cur_window + window
  
  #For Laplacian Exponential Diffusion
  window <- length(par)
  
  best.mjsim <- max(mojoSims[(cur_window + 1):(cur_window + window)])
  
  index <- which(mojoSims[(cur_window + 1):(cur_window + window)]==best.mjsim)
  
  variance <- var(mojoSims[(cur_window + 1):(cur_window + window)])
  
  best[[length(best) + 1]] <- list(par=par[index], value = best.mjsim, variance = variance)
  
  #set the new current window
  cur_window <- cur_window + window
  
  #For Commute-Time Kernel
  window <- 1
  
  best.mjsim <- max(mojoSims[(cur_window + 1):(cur_window + window)])
  
  index <- which(mojoSims[(cur_window + 1):(cur_window + window)]==best.mjsim)
  
  variance <- var(mojoSims[(cur_window + 1):(cur_window + window)])
  
  best[[length(best) + 1]] <- list(par=NULL, value = best.mjsim, variance = variance)
  
  return(best)
  
  best.mjsim = max(mojoSims)
  #   return(mojoSims)
  #   return(which(mojoSims==best.mjsim))
  return(list(mojosim = mojoSims, indices= which(mojoSims==best.mjsim)))
  
}
test.best.graph.kernel <- function(indices=0) {
#     projects <- list("apache-ant-1.9.3", "hadoop-0.20.2", "apache-log4j-1.2.17", "eclipse-jdt-core-3.8",
#                      "jdom-2.0.5", "jedit-5.1.0", "jfreechart-1.2.0", "jhotdraw-7.0.6", "junit-4.12" ,"weka-3.6.11")
#     
#     names(projects) <- c("Apache Ant", "Apache Hadoop", "Apache Log4j", "Eclipse JDT Core", 
#                          "JDOM", "JEdit", "JFreeChart", "JHotDraw", "JUnit", "Weka")
    
#     projects <- list("jedit-5.1.0", "jfreechart-1.2.0", "jhotdraw-7.0.6", "junit-4.12" ,"weka-3.6.11")
#     
#     names(projects) <- c("JEdit", "JFreeChart", "JHotDraw", "JUnit", "Weka")
    
    projects <- list("jdom-2.0.5")
    
    names(projects) <- c("JDOM")
  
  #   pars <- c(10^(-2), 10^(-1), 10^(0), 5, 10^(1), 25, 50, 75, 10^(2), 10^(3))
  
  pars <- c(10^(-3), 10^(-2), 10^(-1), 10^(0), 5, 10^(1), 50, 10^(2))
  
  results <- lapply(projects, function(p) find.best.graph.kernel(p, pars))
  
  
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


plot.influence <- function(influence_matrix_diffusion, influence_matrix_randwalk) {
  
  # TotalInfluenceMatrix ------------------------------------------------
  
  sum_diffusion_matrix <- rowSums(influence_matrix_diffusion)
  
  sum_randwalk_matrix <- rowSums(influence_matrix_randwalk)
  
  
  # CalculateNormalizedDiffusionInfluence -----------------------------------
  
  influence_matrix_diffusion.norm <- normalize.influence(influence_matrix_diffusion)
  
  plotDir = ""
  # PlotDiffusionResults ----------------------------------------------------
  
  # Plot parameters: diffusion
  xlabel='Ranked node'
  ylabel='Influence'
  xmin <- 0
  xmax <- nrow(influence_matrix_diffusion)
  ymin <- 0
  ymax <- max(influence_matrix_diffusion)*1.2
  title <- 'Diffusion Influence Graph'
  
#   plotInfluence(influence_matrix_diffusion, xlabel, ylabel, xmin, xmax, ymin, 
#                 ymax, title, saveDir=plotDir, 
#                 saveFileName='diffusion_influence_graph.pdf', 
#                 sortdecreasing=TRUE)
  
  
  # PlotRandomWalkInfluenceGraph --------------------------------------------
  
  # Plot parameters: random walk
  xlabel='Ranked node'
  ylabel='Average Passage Time'
  xmin <- 0
  xmax <- nrow(influence_matrix_randwalk)
  ymin <- 0
  ymax <- max(influence_matrix_randwalk)*1.2
  title <- 'Random Walk Influence Graph'
  
#   plotInfluence(influence_matrix_randwalk, xlabel, ylabel, xmin, xmax, ymin, 
#                 ymax, title, 
#                 saveDir=plotDir,
#                 saveFileName='randomwalk_influence_graph.pdf',
#                 sortdecreasing=FALSE)
  
  
  # PlotSumsofInfluence -----------------------------------------------------
  
  x11()
  barplot(sum_diffusion_matrix, names.arg=rownames(influence_matrix_diffusion), xlab='node', 
          ylab='total_influence', 
          main='Total influence of each node: diffusion model')
#   dev.copy2pdf(file=file.path(plotDir, 'total_influence_diffusion_model.pdf'))
  
  x11()
  barplot(sum_randwalk_matrix, names.arg=rownames(influence_matrix_randwalk), xlab='node', 
          ylab='total_influence', 
          main='Total influence of each node: random walk model')
#   dev.copy2pdf(file=file.path(plotDir, 'total_influence_randomwalk_model.pdf'))
  
  
  
  # PlotNormalizedInfluence -------------------------------------------------
  
  xlabel='Ranked node'
  ylabel='Normalized Influence'
  xmin <- 0
  xmax <- nrow(influence_matrix_diffusion.norm)
  ymin <- 0
  ymax <- max(influence_matrix_diffusion.norm)*1.2
  title <- 'Diffusion Influence Graph: Normalized by Total Influence'
  
#   plotInfluence(influence_matrix_diffusion.norm, xlabel, ylabel, xmin, xmax, ymin, 
#                 ymax, title, 
#                 saveDir=plotDir,
#                 saveFileName='diffusion_normalized_influence_graph.pdf',
#                 sortdecreasing=TRUE)
  
  x11()
  barplot(rowSums(influence_matrix_diffusion.norm), names.arg=rownames(L), xlab='node', 
          ylab='total normalized influence', main='Total normalized influence of each node: diffusion model')
#   dev.copy2pdf(file=file.path(plotDir, 'total_normalized_influence_diffusion.pdf'))
  
  
  
}

plotInfluence <- function(influence_matrix, xlabel, ylabel, xmin, xmax, ymin, 
                          ymax, title, saveDir, saveFileName, sortdecreasing=TRUE){
  x11()
  
  # Plot blank
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)    # Create extra space on right
  plot(1, type="n", xlab=xlabel, ylab=ylabel, xlim=c(xmin, xmax), 
       ylim=c(ymin, ymax), main=title)
  
  # Set plot vectors
  xvector <- seq(1, nrow(influence_matrix))
  colvector <- palette()[2:length(palette())]
  plotcolorvector = rep(NA, length(xvector))
  
  # Plot in loop
  for(i in 1:nrow(influence_matrix)){
    paletteindex <- i - (i%/%length(colvector)) * (length(colvector)-1)
    shadeindex <- length(colvector) - (i %/% length(colvector) + 1)
    plotcolor <- colors()[grep(colvector[paletteindex], colors())][shadeindex]
    points(xvector, sort(influence_matrix[i, ], decreasing=sortdecreasing),
           type='b', pch=i, 
           col=plotcolor)
    plotcolorvector[i] = plotcolor
  }
  legend('topright', inset=c(-0.28, 0), legend=xvector, 
         pch=1:length(colvector), col=plotcolorvector, title='Influence of Node')
  # legend(0.9*xmax, ymax, xvector, pch=1:length(colvector), col=plotcolorvector)
  
  dev.copy2pdf(file=file.path(saveDir, saveFileName))
}
