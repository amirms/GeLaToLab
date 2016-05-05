#Input



#lex.eval.fun could be constant.string.kernel or normalized_LCU_kernel or normalized_LCS_kernel



compute_semantic_kernel <- function(prname, eval.fun, lex.eval.fun, normalized=F, weights= c(0.5, 0.5), size=0.25, rootFolder = "org", baseline){
  
  require(igraph)
  require(gelato)
  setwd("~/workspace")
  # Read the authoritative decomposition
  decomposition <- read.csv(paste("benchmark", prname ,"decomposition.csv", sep="/"), sep=",",  header = TRUE)
  priori.decomp <- decomposition$x
  names(priori.decomp) <- decomposition$X
  priori.decomp <- normalizeVector(priori.decomp)
  
  
  Adj <- load_SN(prname,make_symmetric = F)
  
  r <- compute_semantic_similarity(prname, Adj, eval.fun, weights)
  
  if (normalized)
    sim_kernel <- r$normalized_sim
  else
    sim_kernel <- r$sim
  
  

#   myBoF = read.csv(paste("benchmark", prname , "BoF", paste(prname, "BoF.csv", sep="-"), sep="/"),  sep = ",")
#   rownames(myBoF) <- myBoF[,1]
#   myBoF <- myBoF[,-1]
#   myBoF <- data.matrix(myBoF)

  #Bag of Features
  myBoF <- load_BoF(prname, c(T,F))

  myBoF <- merge_names_by_lower_case(myBoF, 2)

  
#RUN WHEN PERFORMING EXPERIMENT WITH JEDIT
  noc <- max(priori.decomp)
  
  priori.decomp.names <- unlist(lapply(1:noc, function(x) {
    cls = priori.decomp[priori.decomp==x]
    if (length(cls) > 4)
      names(cls)
    else
      c()
  }))
  
#FIX priori.decomp
#priori.decomp <- priori.decomp[which(names(priori.decomp) %in% priori.decomp.names)]
priori.decomp <- priori.decomp[priori.decomp.names]
priori.decomp <- priori.decomp[order(names(priori.decomp))]

  #Get the sample src code units
  src.code.units <- intersect(rownames(myBoF), names(priori.decomp))
  myBoF <- myBoF[src.code.units,]
  priori.decomp <- priori.decomp[src.code.units]
  
  if (size < 1)
    myBoF <- myBoF[get_sample_docs(prname, priori.decomp, size),]
  
  #Remove unused identifiernames
  
  myBoF <- myBoF[,which(!apply(myBoF,2,FUN = function(x){all(x == 0)}))]
  
  #Find common identifiernames between BoF and the Semantic Network
  identifierNames <- intersect(colnames(myBoF), colnames(sim_kernel))
  
  #Filter out names shorter than 4
  identifierNames <- identifierNames[which(unlist(lapply(identifierNames, nchar))>4)]
  
  #String kernel
#   SK <- lex.eval.fun(identifierNames)
#   dimnames(SK) <- list(identifierNames, identifierNames) 
  
  
  myBoF <- myBoF[,identifierNames]
  sim_kernel <- sim_kernel[identifierNames, identifierNames]
  #remove classes with no identifiers, when combined with the semantic network
  
  myBoF <- myBoF[which(!apply(myBoF,1,FUN = function(x){all(x == 0)})),]
  
  #FIX priori.decomp again  
priori.decomp <- priori.decomp[rownames(myBoF)]


if (!all(rownames(myBoF) == names(priori.decomp)))
  stop('names do not match!')

if (!all(colnames(myBoF) == rownames(sim_kernel)))
  stop('identifier names do not match!')

  #order the names
  # myBoF <- myBoF[order(rownames(myBoF)), ]
  # r$kernel <- r$kernel[order(rownames(r$kernel)), order(colnames(r$kernel))]
  # SK <- SK[order(rownames(SK)), order(colnames(SK))]
  
  #sorting
  # sort(r$kernel[1,], decreasing = T)[1:5]
  
#   semantic <- switch(choice,
#                      SK = SK,
#                      SN = r$kernel,
#                      PR = SK * r$kernel)

#   semantic <- sim_kernel
# S <- sim_kernel

semantic <- sim_kernel
  
  #SVD to compute to USU^T
  USUt <- svd(semantic)
  S <- USUt$u %*% diag(sqrt(USUt$d))
  
  #LSA by reducing concepts
  # d <- rep(0, nrow(USUt$u))
  # D <- USUt$d[USUt$d > 0.7]
  # d[1:10] <- USUt$d[1:10]
  # semantic <- USUt$u %*% diag(d) %*% t(USUt$u)
  
  #diagonal matrix for term weighings
  #TODO CHECK if this is correct
  doc.freq <- colSums(myBoF>0)
  doc.freq[doc.freq == 0] <- 1
  w <- 1/log(nrow(myBoF)/doc.freq)
  R <- diag(w)
  
  
  #Compute cosine similarity
  Phi_d <- myBoF %*% R %*% S
#   dimnames(Phi_d)


  colnames(Phi_d) <- colnames(myBoF)
  
# return(Phi_d)
  
  kernel <- compute_cosine_kernel(Phi_d)
  kernel <- kernel[order(rownames(kernel)), order(colnames(kernel))]
  
  
  
  
  #Fix priori decomposition 
  dummy_v <- rep(0, dim(kernel)[1])
  names(dummy_v) <- rownames(kernel)
  
  priori.decomp <- find.intersection(priori.decomp, dummy_v)
  priori.decomp <- normalizeVector(priori.decomp)
  priori.decomp <- priori.decomp[order(names(priori.decomp))]
  

if (!all(rownames(kernel) == names(priori.decomp)))
  stop('names do not match!')

  #K number of clusters
  noc <- max(priori.decomp)
  
  print("printing the numer of groups:")
  print(max(priori.decomp))
  
  
  #find the intersection with the available classes
  
  #kmeans
  clusters <- kmeans(kernel, centers = noc, iter.max = 1500, nstart = 20000)$cluster
  # result <- spectral.clustering(addition.kernel, k, iter = 15)
  # result <- spectral.clustering(laplacian(addition.kernel, TRUE), k, iter = 40)
  
  clusters <- normalizeVector(clusters)
  
#   names(clusters) = rownames(kernel)
  
  # precision <- compute.precision(clusters, priori.decomp)
  # recall <- compute.recall(clusters, priori.decomp)
  f1.score <- compute.f1(clusters, priori.decomp)
  
  adjustedRI <- compute.AdjRI(clusters, priori.decomp)
  
  mojosim <- compute.MoJoSim(clusters, priori.decomp)
  
  print.mojosim <- paste(round(mojosim, 3), " (+", round((mojosim - baseline$mojosim)*100/baseline$mojosim, 1) , "\\%)", sep="")
  print.f1.score <- paste(round(f1.score, 3), " (+", round((f1.score - baseline$f1.score)*100/baseline$f1.score, 1) , "\\%)", sep="")
  print.adjustedRI <- paste(round(adjustedRI, 3), " (+", round((adjustedRI - baseline$adjustedRI), 3) , ")", sep="")
  
  print(paste(print.mojosim, print.f1.score ,print.adjustedRI, sep="&"))
  
  return(list(mojosim=mojosim, f1.score=f1.score, adjustedRI=adjustedRI))
  
}

output_document_plots <- function(fit, kernel, noc){
  require(ggfortify)
  require(data.table)
  
  require(cluster)
  
  names <- rownames(kernel)
  
  names <- lapply(names, function(name) substr(name, regexpr("/[^/]*$", name)[1] + 1, regexpr("\\.[^\\.]*$", name)[1] - 1))
  
  dimnames(kernel) <- list(names, names)
  
  Phi_d_BoF <- myBoF %*% R
  dimnames(Phi_d) <- dimnames(myBoF)
  dimnames(Phi_d_BoF) <- dimnames(myBoF)
  kernel_BoF <- compute_cosine_kernel(Phi_d_BoF)
  kernel_BoF <- kernel_BoF[order(rownames(kernel_BoF)), order(colnames(kernel_BoF))]
  clusters_BoF <- kmeans(kernel_BoF, centers = noc, iter.max = 1500, nstart = 20000)$cluster
  
  positions <- prcomp(kernel)
  
  dt <- data.table(PC1=positions$x[,1], PC2=positions$x[,2], level=priori.decomp, key="level")
  
  hulls <- dt[, .SD[chull(PC1, PC2)], by = level]
  
  
  #autoplot(pam(kernel, noc), frame = TRUE, label = TRUE, label.size = 3) +  geom_polygon(data = hulls,aes(fill=level,alpha = 0.5))# +  geom_polygon (data=positions, group= clusters_BoF, fill=red) 


  autoplot(pam(kernel, noc), frame = TRUE, label = TRUE, label.size = 3) +
    geom_polygon(data = hulls,fill=NA, colour = "black", alpha = 0.5)
}


#Input: not: no of topics
output_topic_plots <- function(fit, Phi_d, not){
  require(ggfortify)
  require(data.table)
  
  require(cluster)
  
  setwd("~/workspace")
  Phi_d <- read.table(paste("benchmark", prname , "Phi_d.csv", sep="/"), sep=",", row.names = 1, header = TRUE, check.names = FALSE)
  
  names <- colnames(Phi_d)
  
  indices <- c()
  patterns <- c("Time", "Variable",  "error", "Directory",  "jj", "element", "array", "getOut", "macOS", "found",
                "flag", "param", "message", "numLine", "nsName", "complete", "finish",
                "value", "throw", "iterator", "oldstr", "float", "unknown", "instance", "namespace", "usage",
                "provider", "readonly", "tostring", "toclass", "method", "str11", "str21", "short", "double", "warning",
                "synch", "member", "class", "register", "bool", "import", "classpath", "obj", "debug", "statement", "static",
                "package", "type", "dummy", "wrapp", "args", "empty", "isOS", "have", "child", "char", "windows_", 
                "table", "this", "modifier", "resource", "iconst", "default", "starti", "command", "version", "localv",
                "local_", "exit", "prn", "while", "resolve",  "backup1", "backup2", "backups", "replaceWith", "replaceAll", "with", "URI", "temp", "public", "private", "read", "write",
                "connect","state1", "state2", "idx" , "transient", "iswindows" , "notice", "string", "filepath", "nextX",
  
  "binary", "unary", "operation", "linelist", "letter", "primitive", "currentbar", "properties", "invoke", "lineno", "load_data",
  "endline", "toremove", "removeall", "keycode", "setlist", "bracket", "counter", "paths", "freturn", "other", "filename",
  "varname", "islocal", "global", "getIn", "stream", "oldContext", "baseURL", "Literal", "sourceFile", "exact", "operator", "force", "COMMENT1",
  "COMMENT2", "COMMENT3", "COMMENT4", "stub", "addMe", "frameworks", "newParent", "newText", "getField", "IMPLEMENTS", "getDelegated", "isDirty",
  "CATCH", "consArgNames", "tmp", "newDir", "dryRun", "userFile", "subvector", "constructPath", "sourceIn", "newCurrent", "lowerX",
  "TITLE_CASE","initializ", "retVal", "finally", "getsuper", "repeat", "TRUE", "FALSE", "NULL",
  "Inf", "NA_integer_", "NA_real_", "NA_complex_", "NA_character_",   "action", "oldContext", "format", "color1", "color3", "color2", 
  "getColor", "operand", "symbol", "oldCount", "getCount", "setCount", "condExp", "assign","actual", "nextfile",
  "ireturn", "DRETURN", "ARETURN", "LRETURN", "branch", "currLine", "delegate",  "createmodel", "obsolete", "native",
  
  "my", "LPAREN", "RPAREN", "parent", "ID_count", "newfile", "nested", "upper", "lower", "magic1", "magic2", "entry",
  "group", "col_", "row_", "estimate", "which", "getmode", "result", "newline", "labels", "aName","rename",
  "nameSet", "varNames", "newname", "argument", "slash", "private", "protected", "public",
  "LBRACE", "wordLen", "EA_SIZE", "lineNum", "optional", "NameList", "inputHandler", "interf", "STRICTFP",
"newTop", "exception", "MAXLINE", "vfsDst", "getSize", "iteration", "ileft", "iright", "chunkList", "listfiles",
"ForUpdate", "fromRow", "maxItems", "newcount", "support", "migration", "textN", "nothing", "newpath",
"body", "compare", "old_nrows", "newStr", "EA_STATUS", "LAST_ALT", "setLimit", "fileOut", "maximumSize", "CONTINUE",
"path2", "varPattern2", "M_OPEN", "FINAL","getFiles", "iteratee", 
"ncols", "nrows", "old_d", "ntabs", "targetname", "valid_", "outFull", "fully", "parens","mixed",
"continue", "compare", "migration",

"panel2", "getIcon", "source_length", "rowcount", "println", "varPattern", "forward", "toTitleCase",
"memory", "noRecord", "currentSize", "divider", "getrow", "lastLine", "running" , "getCurrent", "settingsDir",
"getStatus", "loadColors", "DOCKING_OPTIONS", "showIcons", "oldIndent", "removeBufferListener", "where", "isLoaded", "stdio", "rootNode",
"COMMA", "uninit", "rhsNode", "end_index", "lname", "isGzipped", "rowNumber",


"replaced", "SUBPIXEL_", "hiByte", "getProperty", "getStyle", "visible_", "ViewCount", "rcfile", "getRules",
"insertion_", "popNode", "fieldVal", "owner", "addMode", "cellText", "loadMode", "addToBus", "getToken", "mkdir",
"UNTITLED",


"SERVICE_NAME", "TRAILING_EOL", "resizing", "AUTORELOAD_DIALOG", "tocRoot", "mapLength", "constructor",  "super",
"getState", "baseLength", "usedBefore", "longTitle", "JAVACODE", "realErr", "newToken", "selRows", "ConditionalExpression",
"indices", "getSource", "locPanel",

"trailingEOL", "PrimaryExpression", "getScope", "snapshot", "dirIcon", "doSuffix", "evaluateCondition", "getTitle", 
"clock", "leftWidth", "country", "lookAndFeel", "keyword", "dispatch", "maximumUpdated", "evaluate", "SCROLL_HORIZ", 
"tokenHandler", "rectSelect", "searchFailed", "encoding", "fieldPanel", "propertyLock", "java", "active1", "identifier",
"worker", "desire", "oldPath", "path1", "topDir", "mkdir", "BOM", "tttext", "mouse", "lhs", "rhs", "piece", "abstract",
"remain", "getBelowPosition", "beginUndo", "fireEndUndo", "requestFocus",  "dockableName", "MINUS", "plus",
"dimension", "normal", "aload", "m_val", "fileExists", "ALOAD", "ILOAD","DLOAD", "LLOAD", "FLOAD", "lfOld",  
"lfNew", "sLfNew", "M_INSERT", "M_OPEN", "invalid", "Offset", "used", "api","teststr", "tester", "DESTROYED",
"hasNext", "hasPrevious", "same", "getLines", "toMerge", "total_", "argNum",
"intfs", "isUnix", "print", "intfs",  "currentText", "byte", "PRIORITY", "newWord", "FILESYSTEM", "_LAYOUT", 
"_LAYER","curPos", "flavor", "install", "new", "next", "root", "prompt", "loaded", "number", "notabs",
"SHIFT", "free_", "startPos", "getLine", "context", "caret_", "oldCaret", "expander", "addSeparator",
"CLOSING", "block","getFromMap", "fromCol","depFrom", "build", "isMac", "curTok", "mnode", "state",
"only", "whole", "fixed", "total", "prefix", "skip", "vector", "TILDE", "task", "runnable", "model",
"modified", "startCaret", "setBounds", "subst", "view")


#TRYING ALOT OF STUFF
#"set", "get", "new", "old", "add", "delete", "open", "close", "begin", "end", "size", "state", "stop", "start",
#"remove", "current", "next", "prev", "old" )


  
  for(i in 1:length(patterns)){
   
    temp <- which(unlist(lapply(names, function(x) grepl(pattern= patterns[i] ,x, ignore.case = T ) )))
    
    print(length(temp))
    
    indices <- c(indices, temp )
    
    
  }
  
  indices <- unique(indices)
#   names[which(unlist(lapply(names, function(x) grepl(pattern= "selecti" ,x, ignore.case = T ) )))]
# "BSH_PACKAGE", "textArea", "regex", "buffer" , "getXml", "submenu", "parse"
  
  Phi_d <- Phi_d[, -indices]



Phi_d <- Phi_d[,which(!apply(Phi_d,2,FUN = function(x){length(which(x > 0)) <= 1}))]

Phi_d <- Phi_d[which(!apply(Phi_d,1,FUN = function(x){all(x == 0)})),]

  not <- 7  

  USVt <- svd(Phi_d)

  D <- diag(sqrt(USVt$d))

  topics <- USVt$v %*% D[,1:not]

  rownames(topics) <- colnames(Phi_d)

  top_identifiers = list()
  
  id_names <- colnames(Phi_d)

  for (i in 1:not){
    
    sorted <- sort(topics[,i], decreasing = T)
    
    top_identifiers[[i]] <- id_names[which(sorted==topics[,i])]
    
    
    
    
  }
  
#MDS
fit <- cmdscale(topics,eig=TRUE, k=2) # k is the number of dim
fit # view results

# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric  MDS",  type="n")
text(x, y, labels = row.names(mydata), cex=.7)


  #Phi(d)
  
  normalized_Phi_d <- div.by.euc.length(t(Phi_d))
  
  d <- distances(normalized_Phi_d)
  
  positions <- prcomp(d)
  
  #plain Phi(d)
  
  plain_Phi_d <- myBoF %*% R
  
  

  dimnames(plain_Phi_d) <- dimnames(myBoF)
  
  normalized_plain_Phi_d <- div.by.euc.length(t(plain_Phi_d))
  
  plain_d <- distances(normalized_plain_Phi_d)
  
  cluster <- kmeans(plain_d, not)$cluster
  
  dt <- data.table(PC1=positions$x[,1], PC2=positions$x[,2], level=cluster, key="level")
  
  hulls <- dt[, .SD[chull(PC1, PC2)], by = level]
  

  compute_ranks <- function(top){
    
    not <- max(top$clustering)
    
    ranks = list()
    
    for (i in 1:not){
      indices <- which(top$silinfo$widths[,1]==i)
      
      sil_widths <- sort(top$silinfo$width[indices, 3], decreasing = T)
      
      ranks[[i]] <- sil_widths
    }
    
    ranks
    
  }

#FIND the ones under 0
wrong_clustered_words <- lapply(ranks, function(r) names(r[which(r < 0)]))
not <- 7

medoids <- which(colnames(d) %in% c("doBsh", "textArea", "regex", "buffer" , "getXml", "submenu", "parse"))

top <- pam(d, not, medoids = medoids)


autoplot(top, frame = TRUE, label = FALSE, label.size = 3)
  
  autoplot(pam(d, not), frame = TRUE, label = FALSE, label.size = 3) + 
    geom_polygon(data = hulls,fill="black", colour = "black", alpha = 0.5)
}


find_exclusive_terms <- function(clusters )

#For each cluster (topics), ranks the terms based on how close they are to the centroid of the cluster
#Input: distance matrix of terms
# not: number of topics
rank_terms <- function(term_dist, not){
  
  euc_dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
  
#   fit <- kmeans(term_dist, centers = not, iter.max = 100, nstart = 1000)
  

fit <- kmeans(term_dist, centers = not, iter.max = 50, nstart = 10)
  ranks = list()
  
  for (i in 1:not){
    
    members_i <- which(fit$cluster==i)
    
    centroid_i <- fit$centers[i,]
    
    ranked_members_i <- unlist(lapply(members_i, function(elem) euc_dist(term_dist[elem,], centroid_i) ))
    
    
    ranks[[i]] <- sort(ranked_members_i, decreasing = F)
    
  }

ranks
  
}



shiny_prep <- function(myBoF, Phi_d, k){
  vocab = colnames(myBoF)
  
  term.frequency <- colSums(myBoF)
  doc.length <- rowSums(myBoF)
  N <- sum(doc.length)
  
  USUt <- svd(Phi_d)
  
  D <- diag(USUt$d)[,1:k]
  
  theta <- USUt$u %*% D
  phi <- t(D) %*% t(USUt$v)
  
#   phi <- t(apply(phi, 1, function(x) scale(x)) ) + 1/dim(myBoF)[2]

phi <- t(apply(phi, 1, function(x)(x-min(x))/(max(x)-min(x)))) + 0.01

phi <- phi / rowSums(phi)
  
theta <- t(apply(theta, 1, function(x)(x-min(x))/(max(x)-min(x)))) + 0.01

theta <- theta / rowSums(theta)
  
  
  json <- createJSON(phi = phi, 
                     theta = theta, 
                     doc.length = doc.length, 
                     vocab = vocab, 
                     term.frequency = term.frequency)


  
}

load_BoF <- function(prname, id2t = c(T, F) ){
  
  
  stopifnot(id2t[1] || id2t[2])
  
#   if (!id2t[1] && !id2t[2])
#     stop("one of the two identifiers or types must be set to TRUE")
  
  setwd("~/workspace")
#   myBoF = read.csv(paste("benchmark", prname , "BoF", paste(prname, "BoF.csv", sep="-"), sep="/"),  sep = ",")
#   rownames(myBoF) <- myBoF[,1]
#   myBoF <- myBoF[,-1]
#   myBoF <- data.matrix(myBoF)
  
  myBoF = read.csv(paste("benchmark", prname , "BoF", paste(prname, "BoF.csv", sep="-"), sep="/"),  sep = ",", header=FALSE, stringsAsFactors = FALSE)
  

  identifiers <- unlist(myBoF[1,2:dim(myBoF)[2]])
  types <- unlist(myBoF[2,2:dim(myBoF)[2]])

  myBoF <- myBoF[-c(1,2),]


  rownames(myBoF) <- myBoF[,1]
  myBoF <- myBoF[,-1]
  myBoF <- data.matrix(myBoF)
  
  #combine the idetifiername and type name
  if (id2t[1] && id2t[2]){
    
    combined_names <- paste(identifiers, types, sep = ":")
    
    colnames(myBoF) <- combined_names
    
    return(myBoF)
  }

  #Only return the identifiers
  if (id2t[1]){
    unique_identifiers <- unique(identifiers)
    
    m <- matrix(0, nrow=dim(myBoF)[1], ncol = length(unique_identifiers))
    
    
    for (i in 1:length(unique_identifiers)){
      
      
      indices <- which(identifiers == unique_identifiers[i])
      
      
      s <- myBoF[,indices]
      
      if (is.vector(s))
        m[,i] <- s
      else
        m[,i] <- rowSums(s)
      
      
    }
    
    colnames(m) <- unique_identifiers
    rownames(m) <- rownames(myBoF)
    
    return(m)
    
  }
  #Only return the types
  if (id2t[2]){
    unique_types <- unique(types)
    
    m <- matrix(0, nrow=dim(myBoF)[1], ncol = length(unique_types))
    
    
    for (i in 1:length(unique_types)){
      
      
      indices <- which(types == unique_types[i])
      
      
      s <- myBoF[,indices]
      
      if (is.vector(s))
        m[,i] <- s
      else
        m[,i] <- rowSums(s)
      
      
    }
    
    colnames(m) <- unique_types
    rownames(m) <- rownames(myBoF)
    
    return(m)
    
  }


return(myBoF)
}

compute_semantic_similarity <- function(prname, Adj, eval.fun, weights){
  #remove Adj of types to identifier names
  #find the first occurence of a type name in the colnames(ADj), set everything else to 0
  names <- colnames(Adj)
  #   names <- rownames(Adj)
  
  #FIND THE STARTING INDEX OF TYPE NAMES  
  classTypeIndex <- which(unlist(gregexpr(pattern ="\\.",names)) > 0)[1]
  primitiveTypeIndices <- which(names %in% c("float", "int", "char", "byte", "void", "double", "boolean"))
  startIndex <- min(c(classTypeIndex, primitiveTypeIndices))
  
  
  print("starting Index")
  print(startIndex)
  
  #size of dictionary
  D <- startIndex -1
  #size of types
  V <- dim(Adj)[2] - startIndex +1
  
  #Also remove contaninmet dependencies from method names to local variables
  #   for (i in 1:dim(Adj)[1])
  #     for (j in 1:(startIndex-1))
  #       Adj[i,j] <- 0
  
  S <- Adj[startIndex:dim(Adj)[1], startIndex:dim(Adj)[2]] 
  dimnames(S) <- dimnames(Adj[startIndex:dim(Adj)[1], startIndex:dim(Adj)[2]])

  print("printing weights")
  print(weights)
  
  
  #DONE make this a higher function argument
  type_sim <- eval.fun(prname, S)
  
  
  if (sum(weights) > 1){
    r <- compute_combined_IPO_ISA_SN(Adj, startIndex, type_sim)
    
    return(list(sim=r$sim, normalized_sim=r$normalized_sim))
  }
  
  if (weights[1] > 0){
    ISA <- compute_ISA_SN(Adj, startIndex, type_sim)

  }
  else
    ISA <- list(sim=0, normalized_sim=0)
  
  
  
  if (weights[2] > 0)
    IPO <- compute_IPO_SN(Adj, startIndex, type_sim)
  else
    IPO <- list(sim=0, normalized_sim=0)
  
  #then combine the two Synsets and get one similarity matrix
 
  
  list(sim=(weights[1] * ISA$sim + weights[2] * IPO$sim), normalized_sim = (weights[1] * ISA$normalized_sim + weights[2] * IPO$normalized_sim) )
}

#combine combined ISA-IPO semantic similarity between terms for SN
#eval.fun is the evaluation function
compute_combined_IPO_ISA_SN <- function(Adj, startIndex, type_sim){
  names <- rownames(Adj)
  #   names <- rownames(Adj)
  
  outIndices <- list()
  outISA_IPO <- list()
  outdegree <- c()
  
  for (i in 1:(startIndex-1)){
    
    outIndices[[i]] <- c(which(Adj[startIndex:dim(Adj)[1], i]>0), which(Adj[i, startIndex:dim(Adj)[2]]>0))
    
    outISA_IPO[[i]] <- c(Adj[names(outIndices[[i]]),i], Adj[i,names(outIndices[[i]])])
    
    outdegree[i] <- sum(outISA_IPO[[i]])        
    
  }
  
  result <- compute_term_similarity(outIndices, outISA_IPO, outdegree, type_sim)
  
  dimnames(result$sim) <- dimnames(Adj[1:(startIndex-1), 1:(startIndex-1)])
  dimnames(result$normalized_sim) <- dimnames(Adj[1:(startIndex-1), 1:(startIndex-1)])
  
  result
  
}

#compute ISA semantic similarity between terms for SN
#eval.fun is the evaluation function
compute_IPO_SN <- function(Adj, startIndex, type_sim){
  
  names <- rownames(Adj)
  #   names <- rownames(Adj)
  
  outIndices <- list()
  outIPO <- list()
  outdegree <- c()
  
  for (i in 1:(startIndex-1)){
    
    outIndices[[i]] <- which(Adj[startIndex:dim(Adj)[1], i]>0)
    
    outIPO[[i]] <- Adj[names(outIndices[[i]]),i]
    
    outdegree[i] <- sum(outIPO[[i]])        
    
  }
  
  result <- compute_term_similarity(outIndices, outIPO, outdegree, type_sim)
  
  dimnames(result$sim) <- dimnames(Adj[1:(startIndex-1), 1:(startIndex-1)])
  dimnames(result$normalized_sim) <- dimnames(Adj[1:(startIndex-1), 1:(startIndex-1)])

  result
}
  
compute_term_similarity  <- function(outIndices, outlinks, outdegree, type_sim){
  
  #length of all parameters must mach
  len <- length(outIndices)
  
  sim <- matrix(0, nrow=len, ncol=len)
  normalized_sim <- matrix(0, nrow=len, ncol=len)
  
  for (i in 1:len){
    outIndices1 <- outIndices[[i]]
    if (length(outIndices1) < 1)
      next
    
    outlinks1 <- outlinks[[i]]
    outdegree1<- outdegree[i]
    
    for (j in 1:len){
      
      outIndices2 <- outIndices[[j]]
      if (length(outIndices2) < 1)
        next
      
      outlinks2 <- outlinks[[j]]
      outdegree2<- outdegree[j]
      
      if (length(intersect(outIndices1,outIndices2))>0){
        sim[i,j] <- 1
        normalized_sim[i,j] <- 1
        next
      }
      
      max_sim <- 0 
      max_normalized_sim <- 0
      
      for(k in 1:length(outIndices1)){
        for(l in 1:length(outIndices2)){
          
          sim_val <- type_sim[outIndices1[k], outIndices2[l]]
          
          normalized_sim_val <- 0.5 * sim_val * ((outlinks1[k]/outdegree1) + (outlinks2[l]/outdegree2))
          
          if (sim_val > max_sim)
            max_sim <- sim_val
          
          if (normalized_sim_val > max_normalized_sim)
            max_normalized_sim <- normalized_sim_val
          
        }
      }
      sim[i,j] <- max_sim
      normalized_sim[i,j] <- max_normalized_sim
    }
    
  }
  
  
#   sim <- fill_lower_diagonal(sim)
#   normalized_sim <- fill_lower_diagonal(normalized_sim)
  
  return(list(sim=sim, normalized_sim=normalized_sim))
  
  
}

#compute ISA conceptual density between terms for SN
#eval.fun is the evaluation function
compute_ISA_SN <- function(Adj, startIndex, type_sim){
  
  names <- rownames(Adj)
  #   names <- rownames(Adj)
  
  outIndices <- list()
  outISA <- list()
  outdegree <- c()
  
  for (i in 1:(startIndex-1)){
    
    outIndices[[i]] <- which(Adj[i, startIndex:dim(Adj)[2]]>0)
    
    outISA[[i]] <- Adj[i,names(outIndices[[i]])]
    
    outdegree[i] <- sum(outISA[[i]])        
      
  }
  
  result <- compute_term_similarity(outIndices, outISA, outdegree, type_sim)
  
  dimnames(result$sim) <- dimnames(Adj[1:(startIndex-1), 1:(startIndex-1)])
  dimnames(result$normalized_sim) <- dimnames(Adj[1:(startIndex-1), 1:(startIndex-1)])
  
  result

}


compute_Resnik_similarity <- function(prname, S){

  
  BoT <- load_BoF(prname, id2t = c(F,T))
  
  SN_type_names <- colnames(S)
  
  BoT_type_names <- colnames(BoT)
  
  common_type_names <- intersect(BoT_type_names, SN_type_names)

  #FIXME should I make the intersection
# BoT <- BoT[,typenames]
  
  sum_BoT <- colSums(BoT)
  
  information_content <- -log(sum_BoT /sum(sum_BoT))
  

sim_score <- matrix(0, nrow=dim(S)[1], ncol=dim(S)[2])
  
  for (i in 1:dim(S)[1])
    for(j in i:dim(S)[2])
      if (i != j)
      {
        
        #check if one is subtype of other
        LCHs <- compute_least_common_hypernym(S, i, j)
        
        if (length(LCHs) < 1)
          next
        
        maxSim <- 0
        
        for (k in 1 : length(LCHs)){
          
          LCH <- LCHs[k]
          
          if (!(SN_type_names[LCH] %in% BoT_type_names))
            next
          
          sim <- information_content[SN_type_names[LCH]]
          
          if (sim > maxSim)
            maxSim <- sim
          
        }
        
        
        sim_score[i,j] <- maxSim
        
      }
  
sim_score
  
}



compute_Lin_similarity <- function(prname, S){
  BoT <- load_BoF(prname, id2t = c(F,T))
  
  SN_type_names <- colnames(S)
  
  BoT_type_names <- colnames(BoT)
  
  common_type_names <- intersect(BoT_type_names, SN_type_names)
  
  #FIXME should I make the intersection
  # BoT <- BoT[,typenames]
  
  sum_BoT <- colSums(BoT)
  
  information_content <- -log(sum_BoT /sum(sum_BoT))
  
  
  sim_score <- matrix(0, nrow=dim(S)[1], ncol=dim(S)[2])
  
  for (i in 1:dim(S)[1])
    for(j in 1:dim(S)[2])
      if (i != j)
      {
        
        if (!(SN_type_names[i] %in% BoT_type_names && SN_type_names[j] %in% BoT_type_names))
          next
        
        
        #check if one is subtype of other
        LCHs <- compute_least_common_hypernym(S, i, j)
        
        if (length(LCHs) < 1)
          next
        
        maxSim <- 0
        
        for (k in 1 : length(LCHs)){
          
          LCH <- LCHs[k]
          
          if (!(SN_type_names[LCH] %in% BoT_type_names))
            next
          
          sim <- information_content[SN_type_names[LCH]]
          sim <- 2 * information_content[SN_type_names[LCH]] / (information_content[SN_type_names[i]] + information_content[SN_type_names[j]])
          
          if (sim > maxSim)
            maxSim <- sim
          
        }
        
        
        sim_score[i,j] <- maxSim
        
      }
  
  sim_score
  
}


compute_Wu_Palmer_similarity <- function(prname, S, rootNode ="java.lang.Object"){
  
  require(igraph)
  
  g <- graph.adjacency(S, weighted=TRUE, mode="directed", diag=F)
  
  #this could be average, max and min
  SP <- shortest.paths(g, mode="all")
  
  SP[is.infinite(SP)] <- 0 
  
  sim_matrix <- matrix(0, nrow=dim(S)[1], ncol=dim(S)[2])
  
  root <- which(rownames(S)==rootNode)
  
  for (i in 1:dim(S)[1])
    for(j in i:dim(S)[2])
      if (i != j)
      {
        
        #check if one is subtype of other
        LCHs <- compute_least_common_hypernym(S, i, j)
        
        if (length(LCHs) < 1)
          next
        
        
        maxSim <- 0
        
        for (k in 1 : length(LCHs)){
          
          LCH <- LCHs[k]
          depth <- compute_depth(S, root, LCH)
#           print(depth)
#           print(depth)
          distance1 <- SP[i, LCH]
#           print(distance1)
          distance2 <- SP[j, LCH]
# print(distance2)
          sim <- 2*depth/(distance1 + distance2 + 2*depth)
          
          if (sim > maxSim)
            maxSim <- sim
          
        }
        
        
        sim_matrix[i,j] <- maxSim
        
      }
  
  
  sim_matrix[lower.tri(sim_matrix)] <- t(sim_matrix)[lower.tri(sim_matrix)]
  sim_matrix
  
  
}


compute_Leacock_Chodorow_similarity <- function(prname, S, rootNode ="java.lang.Object"){
  
  require(igraph)
  
  g <- graph.adjacency(S, weighted=TRUE, mode="directed", diag=F)
  
  #this could be average, max and min
  SP <- shortest.paths(g, mode="all")
  
  SP[is.infinite(SP)] <- 0 
  
  root <- which(rownames(S)==rootNode)
  
  sim_matrix <- matrix(0, nrow=dim(S)[1], ncol=dim(S)[2])
  
  for (i in 1:dim(S)[1])
    for(j in 1:dim(S)[2])
      if (i != j)
      {
          
        max_depth <- max(compute_depth(S, root, i), compute_depth(S, root, j))
        distance <- SP[i, j]
          
        sim <- -log(distance/(2*max_depth))

        sim_matrix[i,j] <- sim
        
      }
  
  sim_matrix
  
  
}


compute_inverted_path_length <- function(prname, S, alpha=1){
 
  require(igraph)
  
  g <- graph.adjacency(S, weighted=TRUE, mode="directed", diag=F)
  
  #this could be average, max and min
  SP <- shortest.paths(g, mode="all")
  
  SP[is.infinite(SP)] <- 0 
  
  1 / (1 + SP^alpha)
}

#Input: Adj is a directed graph

#FIXME there has to be a better way to find the starting Index for the types
#TODO what to do with subtypes
compute_conceptual_density <- function(prname, S){
  
  branching_factor <- compute_branching_factor(S)
  
  CD <- matrix(0, nrow=dim(S)[1], ncol=dim(S)[2])

  for (i in 1:dim(S)[1])
    for(j in i:dim(S)[2])
#   for (i in 1:20)
#     for(j in 1:20)
      if (i != j)
        {
      
        #check if one is subtype of other
        LCHs <- compute_least_common_hypernym(S, i, j)
        
        if (length(LCHs) < 1)
          next
        
        maxCD <- 0
        
        for (k in 1 : length(LCHs)){
          
          LCH <- LCHs[k]
          
          cd <- calculate_conceptual_density(LCH, branching_factor)
          
          if (cd > maxCD)
            maxCD <- cd
          
        }
        
      
        CD[i,j] <- maxCD
        
      }
  
  CD <- fill_lower_diagonal(CD)

  CD
}


fill_lower_diagonal <- function(sim_matrix){
  sim_matrix[lower.tri(sim_matrix)] <- t(sim_matrix)[lower.tri(sim_matrix)]
  
  sim_matrix
  
}


#branching factor is a list ()

calculate_conceptual_density <- function(LCH, branching_factor){
    
  total_branching_factor <- branching_factor[[LCH]]$total_branching
#   print(total_branching_factor)
  total_children <- branching_factor[[LCH]]$total_children
#   print(total_children)
  avg_branching_factor <- total_branching_factor / total_children
  
  #Output: h for calculating the conceptual density
  estimate_depth <- function(avg_branching_factor){
    if (avg_branching_factor== 1)
      return(2)
    
    floor(log(2, base = avg_branching_factor))
    
  }
  
  
  h <- estimate_depth(avg_branching_factor)
  
  sum_squared_avg_branching_factor <- 0
  
  for (i in 0:h){
    
    sum_squared_avg_branching_factor <- sum_squared_avg_branching_factor + avg_branching_factor^i
  }
  
  sum_squared_avg_branching_factor/total_branching_factor
  
}


compute_transitive_closure <- function(graph){
  
  # reach[][] will be the output matrix that will finally have the 
  #shortest distances between every pair of vertices */
  #     int reach[V][V], i, j, k;
  
  if (dim(graph)[1] != dim(graph)[2])
    stop("Incompatible dimensions of the graph")
  
  V <- dim(graph)[1]
  
  reach <- matrix(0, nrow=V, ncol=V)
  
  
  # Initialize the solution matrix same as input graph matrix. Or
  #we can say the initial values of shortest distances are based
  #on shortest paths considering no intermediate vertex. */
  reach <- graph
  
  # Add all vertices one by one to the set of intermediate vertices.
  #---> Before start of a iteration, we have reachability values for
  #all pairs of vertices such that the reachability values 
  #consider only the vertices in set {0, 1, 2, .. k-1} as 
  #intermediate vertices.
  #----> After the end of a iteration, vertex no. k is added to the 
  #set of intermediate vertices and the set becomes {0, 1, .. k} */
  for (k in 1:V)
  {
    # Pick all vertices as source one by one
    for (i in 1:V)
    {
      # Pick all vertices as destination for the
      # above picked source
      for (j in 1:V)
      {
        # If vertex k is on a path from i to j,
        # then make sure that the value of reach[i][j] is 1
        reach[i][j] = reach[i][j] || (reach[i][k] && reach[k][j]);
      }
    }
  }
  
  # Print the shortest distance matrix
  return(reach);  
  
  
}



#Input: root, x, y should be index in adj
#Output: a list of nodes, one for each node, excluding 
compute_subhierarchy <- function(adj, root, x, y){
  require(igraph)

  g <- graph.adjacency(adj, mode='directed', diag=FALSE)  
  
  p <- get.shortest.paths(g, root, to = c(x,y), mode = c("in"))
  
  unique(unlist(p$vpath))
  
}

#FIXME: there could be another way of implementing this by finding all the reachable paths from a node, then finding the one
#minimal distance that reaches one possible root Node
#In case where a unique top node exists, the following algorithm: the shortest path from root node to the ther node suffices
#Input: Adjacency, 
#Output: 
compute_depth <- function(adj, root, x){
  
  require(igraph)
  
  g <- graph.adjacency(adj, mode='directed', diag=FALSE)
  
  shortest.paths(g, v=root, to=x, mode="in")
}


#Input: an adjacency graph Adj, LCA(x,y)
#Output single least common ancestors of x and y
compute_least_common_hypernym <- function(adj, x, y){
  
  #Let say you want to compute LCA(x,y) with x and y two nodes. Each node must have a value color and count, resp. initialized to white and 0.
  
  #Color all ancestors of x as blue (can be done using BFS)
  #Color all blue ancestors of y as red (BFS again)
  #For each red node in the graph, increment its parents' count by one
  #Each red node having a count value set to 0 is a solution.
  
  require(igraph)
  
  g <- graph.adjacency(adj, mode='directed', diag=FALSE,weighted=T)  
  
  V <- dim(adj)[1]
  
  blue <- subcomponent(g, x, "out")
  
  temp <- subcomponent(g, y, "out")
  
  
  red <- blue[blue %in% temp]
  
#   indices <- which(red)
  count <- rep(0, V)
  
  
  for (i in 1:length(red))
    count[which(adj[red[i],]>0)] <- count[which(adj[red[i],]>0)] + 1
  
#FIXME change this
if (length(red) < 1)
  return (c())

for (i in 1:length(red))
  if (count[red[i]] == 0 )
    return(which(rownames(adj) == red[i]$name))

}


#Input: adjacency matrix of a directed acyclic graph
#Output: computes branching factor for each node (type)
compute_branching_factor <- function(adj){
  
  require(igraph)
  
  g <- graph.adjacency(adj, mode='directed', diag=FALSE,weighted=T)  
  
  V <- dim(adj)[1]
  
  branching <- colSums(adj)
  
  r <- c()
  
  for (i in 1:V){
    reachable_nodes <- subcomponent(g, i, "in")
    
    
    total_children <- length(reachable_nodes)
    
    total_branching <- sum(branching[reachable_nodes])
    
    r[[i]] <- list(total_branching=total_branching, total_children=total_children)
    
  }
  
  r
    
}


#Input: adjacency matrix of a directed acyclic graph
#Precondition: Only works when it is a tree
#Output: computes branching factor for each node (type)
# compute_branching_factor <- function(adj){
#   
#   V <- dim(adj)[2]
#   
#   r=c()  
#   for (i in 1:V){
#     r[[i]] <- list(total_branching=0, total_children=0, inQueue=F)
#   }
#   
#   #Is it directed?!?  
#   #Find the leaves and add them to the queue
#   q=queue()
#   
#   leaves <- which(colSums(adj) < 1)
#   
#   push.queue(q, leaves)
#   
#   while (!is.empty.queue(q)) {
#     node <- pop.queue(q)
#     r[[node]]$inQueue <- FALSE
#     
# #     print(node)
#     
#     parents <- which(adj[node, ]>0)
#     
#     if (length(parents) < 1)
#       next
# 
#     for (i in 1:length(parents)){
#       
#       if (node == parents[[i]])
#         next
#         
# #       print(parents[[i]])
#       r[[parents[[i]]]]$total_branching <-   r[[node]]$total_branching + adj[node, parents[[i]]]
#       
#       r[[parents[[i]]]]$total_children <- r[[node]]$total_children + 1
# 
#     if (!r[[parents[[i]]]]$inQueue){
#         push.queue(q,parents[i])
#         
#         r[[parents[[i]]]]$inQueue <- TRUE
#       }
#       
#     }
#     
#   }
# 
#   r
#   
# }


normalized_LCU_kernel <- function(names){
  
  len <- length(names)
  
  m = matrix(0, nrow=len, ncol=len)
  
  for (i in 1:len)
    for(j in i:len)
      if (i != j)
        m[i,j] <- compute_normalized_LCU(names[i], names[j])
  
  m <- fill_lower_diagonal(m)
  
  m
}


normalized_LCS_kernel <- function(names){
  
  len <- length(names)
  
  m = matrix(0, nrow=len, ncol=len)
  
  for (i in 1:len)
    for(j in i:len)
      if (i != j)
        m[i,j] <- compute_normalized_LCS(names[i], names[j])
  
  m <- fill_lower_diagonal(m)
  
  m
}

#Input: two strings str1, str2
#Output: length_longest_common_substring(str1, str2)^2/(len(str1)*len(str2))

compute_normalized_LCU <-function(str1, str2){
  llcu <- length_longest_common_substring(str1, str2)
  
  llcu^2/(nchar(str1) * nchar(str2))
}

#Input: two strings str1, str2
#Output: LCS(str1, str2)^2/(len(str1)*len(str2))
compute_normalized_LCS <- function(str1, str2){
  llcs <- LCS(str1, str2)
  
  llcs^2/(nchar(str1) * nchar(str2))
  
}

#prnames <- c("")
run_in_batch_mode <- function(prnames){
  
  for(i in 1:length(prnames)){
    print("")
    print(prname)
    run_in_batch_mode_each_project(prname)
    
  }
  
}


run_in_batch_mode_each_project <- function(prname, size=0.25){
  
  run_each_setting <- function(weights, eval.fun=NULL, lex.fun =NULL,  Adj, priori.decomp, myBoF){
    require(igraph)
    require(gelato)
    
    sim_kernel <- NULL
    
    if (!is.null(eval.fun)){
    
      r <- compute_semantic_similarity(prname, Adj, eval.fun, weights)
      
      if (normalized)
        sim_kernel <- r$normalized_sim
      else
        sim_kernel <- r$sim
    }
    
    #Calculate the identifier names et
    names <- colnames(Adj)
      
    #FIND THE STARTING INDEX OF TYPE NAMES  
    classTypeIndex <- which(unlist(gregexpr(pattern ="\\.",names)) > 0)[1]
    primitiveTypeIndices <- which(names %in% c("float", "int", "char", "byte", "void", "double", "boolean"))
    startIndex <- min(c(classTypeIndex, primitiveTypeIndices))
      
      
    print("starting Index")
    print(startIndex)
      
    #size of dictionary
    D <- startIndex -1
      
    identifierNames <- rownames(Adj)[1:D]
    #Find common identifiernames between BoF and the Semantic Network
    identifierNames <- intersect(colnames(myBoF), identifierNames)
    
    #Filter out names shorter than 4
    identifierNames <- identifierNames[which(unlist(lapply(identifierNames, nchar))>4)]
    
    
    myBoF <- myBoF[,identifierNames]
    
    if (!is.null(sim_kernel))
      sim_kernel <- sim_kernel[identifierNames, identifierNames]
    else
      sim_kernel <- matrix(1, nrow=length(identifierNames), ncol=length(identifierNames))
    #remove classes with no identifiers, when combined with the semantic network
    
    myBoF <- myBoF[which(!apply(myBoF,1,FUN = function(x){all(x == 0)})),]
    myBoF <- myBoF[order(rownames(myBoF)), ]
    
    
    if (!is.null(lex.fun))
        string_kernel <- lex.fun(identifierNames)
    else
      string_kernel <- matrix(1, nrow=length(identifierNames), ncol=length(identifierNames))
    
    #Element-wise product  
    semantic <- string_kernel * sim_kernel
    
    # semantic <- SK
    
    #SVD to compute to USU^T
    USUt <- svd(semantic)
    S <- USUt$u %*% diag(sqrt(USUt$d))
    
    #LSA by reducing concepts
    # d <- rep(0, nrow(USUt$u))
    # D <- USUt$d[USUt$d > 0.7]
    # d[1:10] <- USUt$d[1:10]
    # semantic <- USUt$u %*% diag(d) %*% t(USUt$u)
    
    #diagonal matrix for term weighings
    #TODO CHECK if this is correct
    doc.freq <- colSums(myBoF>0)
    doc.freq[doc.freq == 0] <- 1
    w <- 1/log(nrow(myBoF)/doc.freq)
    R <- diag(w)
    
    
    #Compute cosine similarity
    Phi_d <- myBoF %*% R %*% S
    
    
    
    dimnames(Phi_d) <- dimnames(myBoF)
    
    # return(Phi_d)
    
    kernel <- compute_cosine_kernel(Phi_d)
    kernel <- kernel[order(rownames(kernel)), order(colnames(kernel))]
    
    
    
    
    #Fix priori decomposition 
    dummy_v <- rep(0, dim(kernel)[1])
    names(dummy_v) <- rownames(kernel)
    
    priori.decomp <- find.intersection(priori.decomp, dummy_v)
    priori.decomp <- normalizeVector(priori.decomp)
    priori.decomp <- priori.decomp[order(names(priori.decomp))]
    
    if(!all(rownames(kernel) == names(priori.decomp)))
      stop("names don't match!")
    
    #K number of clusters
    noc <- max(priori.decomp)
    
    print("printing the numer of groups:")
    print(max(priori.decomp))
    
    
    #find the intersection with the available classes
    
    #kmeans
    clusters <- kmeans(kernel, centers = noc, iter.max = 1500, nstart = 20000)$cluster
    # result <- spectral.clustering(addition.kernel, k, iter = 15)
    # result <- spectral.clustering(laplacian(addition.kernel, TRUE), k, iter = 40)
    
    clusters <- normalizeVector(clusters)
    
#     names(clusters) = rownames(kernel)
    
    # precision <- compute.precision(clusters, priori.decomp)
    # recall <- compute.recall(clusters, priori.decomp)
    f1.score <- compute.f1(clusters, priori.decomp)
    
    adjustedRI <- compute.AdjRI(clusters, priori.decomp)
    
    mojosim <- compute.MoJoSim(clusters, priori.decomp)
    
    return(list(mojosim=mojosim, f1.score=f1.score, adjustedRI=adjustedRI))
  }
  
  
  require(igraph)
  require(gelato)
  setwd("~/workspace")
  # Read the authoritative decomposition
  decomposition <- read.csv(paste("benchmark", prname ,"decomposition.csv", sep="/"), sep=",",  header = TRUE)
  priori.decomp <- decomposition$x
  names(priori.decomp) <- decomposition$X
  priori.decomp <- normalizeVector(priori.decomp)
  
  
  Adj <- load_SN(prname,make_symmetric = F)

  #Bag of Features
  myBoF <- load_BoF(prname, c(T,F))
  
  myBoF <- merge_names_by_lower_case(myBoF, 2)
  
  
  #Get the sample src code units
  src.code.units <- intersect(rownames(myBoF), names(priori.decomp))
  myBoF <- myBoF[src.code.units,]
  priori.decomp <- priori.decomp[src.code.units]
  
  if (size < 1)
    myBoF <- myBoF[get_sample_docs(prname, priori.decomp, size),]
  
  #Remove unused identifiernames
  
  myBoF <- myBoF[,which(!apply(myBoF,2,FUN = function(x){all(x == 0)}))]
  
  
  #populate different settings
  
  eval.funs <- c(compute_inverted_path_length, compute_Wu_Palmer_similarity, compute_Leacock_Chodorow_similarity, 
                 compute_conceptual_density, compute_Lin_similarity, compute_Resnik_similarity )
  
  weights <- list(c(0,1), c(1,0), c(0.5,0.5), c(1,1))
  
  results = list()
  
  for (i in 1:length(weights)){
    print(weights[[i]])
    results = c()
    for (j in 1:length(eval.funs)){
      
      r <- run_each_setting(weights[[i]], eval.funs[[j]], NULL, Adj, priori.decomp, myBoF)
      
      print.mojosim <- round(r$mojosim, 3)
      print.f1.score <- round(r$f1.score, 3)
      print.adjustedRI <-round(r$adjustedRI, 3)
      
      results[j] <- paste(print.mojosim, print.f1.score ,print.adjustedRI, sep="&")
      
      print(results[j])
      
    }
    
    
    if (i==1)
      write(results, file = paste("benchmark", prname ,"IPO.txt", sep="/"))
    else if (i==2)
      write(results, file = paste("benchmark", prname ,"ISA.txt", sep="/"))
    else if (i==3)
      write(results, file = paste("benchmark", prname ,"ISA-IPO.txt", sep="/"))
    else
      write(results, file = paste("benchmark", prname ,"Combined.txt", sep="/"))
    
    
  }

  lex.funs <- c(normalized_LCS_kernel, normalized_LCU_kernel, constant.string.kernel)

  for (i in 1:length(lex.funs)){
    
    r <- run_each_setting(c(1,1), NULL, lex.funs[[i]], Adj, priori.decomp, myBoF)
    
    print.mojosim <- round(r$mojosim, 3)
    print.f1.score <- round(r$f1.score, 3)
    print.adjustedRI <-round(r$adjustedRI, 3)
    
    result <- paste(print.mojosim, print.f1.score , print.adjustedRI, sep="&")
    
    print(result)
    
    
    if (i==1)
      write(result, file = paste("benchmark", prname ,"LCS.txt", sep="/"))
    else if (i==2)
      write(result, file = paste("benchmark", prname ,"LCU.txt", sep="/"))
    else
      write(result, file = paste("benchmark", prname ,"Constant.txt", sep="/"))
    
  }

}

#     projects <- list("apache-ant-1.9.3", "hadoop-0.20.2", "apache-log4j-1.2.17",
#"jdom-2.0.5", "jedit-5.1.0", "jfreechart-1.2.0", "jhotdraw-7.0.6", "junit-4.12" ,"weka-3.6.11")