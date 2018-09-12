prname = "ant";
prname = "hadoop";
prname = "jdom";
prname = "eclipse.jdt.core";
prname = "jfreechart";
prname ="jhotdraw";
prname ="junit4";
prname = "log4j";
prname = "weka";

parseLogHistory <- function(prname, rootFolder = "org/"){
  filename = paste("c:\\src\\git\\benchmarks\\",prname,"\\log.txt", sep="" );
  #x <- read.csv(filename);
  
  
  oracle <- readLines(filename)
  nvec <- length(oracle)
  breaks <- which(! nzchar(oracle))
  nbreaks <- length(breaks)
  if (breaks[nbreaks] < nvec) {
    breaks <- c(breaks, nvec + 1L)
    nbreaks <- nbreaks + 1L
  }
  
  changesPerCommit = list();
  
  changes =c();
  index = 1;
  for(i in 1:nvec){
    if (i %in% breaks) {
      changesPerCommit[[index]] = changes;
      changes = c();
      index = index + 1;
      next;
    }
    changes = c(changes, oracle[i]);
  }
  
  x <- lapply(changesPerCommit, function(cs) {
    z <- cs[grepl(rootFolder, cs)]
    
    if (length(z) > 0){
      return(z);
    }
    
    return(NULL);
  })
  
  x[sapply(x, is.null)] <- NULL
  
  names <- unlist(x);
  names <- lapply(names, function(strs) substring(strs, gregexpr(pattern =rootFolder, strs)[[1]][length(gregexpr(pattern =rootFolder, strs)[[1]])]))
  
  names <- unlist(names)
  names <- names[!duplicated(names)]
  
  m <- matrix(0, nrow=length(names), ncol=length(x), dimnames= list(names, paste("T", seq(1:length(x)), sep="")))
  
  for (i in 1:length(x)){
    z = unlist(x[[i]]);
    z = unlist(lapply(z, function(cs) substring(cs, gregexpr(pattern =rootFolder, cs)[[1]][length(gregexpr(pattern =rootFolder, cs)[[1]])])))
    indices = which(names %in% z)  
  
    m[indices, i] <- 1
  }
  
  m <- m[,colSums(m)>1]
  m <- m[rowSums(m)>0,]
  
  m[which(m==0)] = NaN
  
  
  write.table(m, file =paste("c:\\src\\git\\benchmarks\\",prname,"\\commit-story.csv", sep=""),row.names=TRUE, col.names=NA,sep=",", quote=TRUE,na="")
}