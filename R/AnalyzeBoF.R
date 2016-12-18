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
                "local_", "exit", "prn", "while", "resolve",  "backup1", "backup2", "backups", "replaceWith", "replaceAll", "with", 
                "URI", "temp", "public", "private", "read", "write",
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
  
  
  
  
  
  
  Phi_d <- Phi_d[,which(!apply(Phi_d,2,FUN = function(x){length(which(x != 0)) <= 1}))]
  
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



compute_t_SNE <- function(mydata){
  library(tsne)
  tsne_mydata <- tsne(mydata, k=2, max_iter=500, epoch=100)
  
  rownames(tsne_mydata) <- rownames(mydata)
  
  tsne_mydata
}

# lapply(rank_euc, function(r) r[1:7])

#TODO find exclusive terms


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