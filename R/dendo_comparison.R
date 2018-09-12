
convertDependencies2Sim <- function(dependencies, beta=0.5, identifiers=c()){

  m1 <- dependencies[[1]]
  
  for (i in 2:length(dependencies)){
    m2 <- dependencies[[i]]
    tmp = rbind(as.data.frame(as.table(m1)), as.data.frame(as.table(m2)))
    m1 <- xtabs(Freq ~ Var1 + Var2, tmp)
  }
  
  require(igraph)
  #   Adj <- load_SN(prname)
  
  g <- build.graph(m1)
  r <- graph.diffusion(g, beta= beta, v=V(g))
  # dimnames(r$kernel) <- dimnames(Adj)
  dimnames(r$dist) <- dimnames(Adj)
  
  
  if (length(identifiers) > 0)
    r$dist <- r$dist[identifiers, identifiers]
  
  return(r$dist)
}



#Input: sim = list(sim_kernel, string_kernel, myBoF)
x <- function(sim){
  require(vegan)
  
  #unbox
  myBoF <- sim$myBoF
  string_kernel <- sim$string_kernel
  sim_kernel <- sim$sim_kernel
  
  
  non_duplicated <- !duplicated(t(myBoF))
  
  myBoF <- myBoF[, non_duplicated]
  
  sim_kernel <- sim_kernel[non_duplicated, non_duplicated]
  string_kernel <- string_kernel[non_duplicated, non_duplicated]
  
  semantic <- string_kernel * sim_kernel
  
  USUt <- svd(semantic)
  S <- USUt$u %*% diag(sqrt(USUt$d))
  
  #diagonal matrix for term weighings
  #TODO CHECK if this is correct
  doc.freq <- colSums(myBoF>0)
  doc.freq[doc.freq == 0] <- 1
  
  #   term.freq <- rowSums(myBoF)
  #   term.freq[term.freq == 0] <- 1
  
  #   w <- 1/log(nrow(myBoF)/doc.freq)
  w <- log(doc.freq/nrow(myBoF))
  R <- diag(w)
  
  #Compute enriched Phi_d
  enriched_Phi_d <- myBoF %*% R %*% S
  dimnames(enriched_Phi_d) <- dimnames(myBoF)
  enriched_Phi_d <- enriched_Phi_d[order(rownames(enriched_Phi_d)),]
  
#   names <- get.sample.names(colnames(enriched_Phi_d))
# names <- manual.names()
#   enriched_Phi_d <- enriched_Phi_d[,names]
  
  enriched_Phi_d <- enriched_Phi_d[,which(!apply(enriched_Phi_d,2,FUN = function(x){length(which(x != 0)) <= 1}))]
  enriched_Phi_d <- enriched_Phi_d[which(!apply(enriched_Phi_d,1,FUN = function(x){all(x == 0)})),]
  

  #Compute plain Phi_d


  #diagonal matrix for term weighings
  #TODO CHECK if this is correct
  doc.freq <- colSums(myBoF>0)
  doc.freq[doc.freq == 0] <- 1
  
  #   term.freq <- rowSums(myBoF)
  #   term.freq[term.freq == 0] <- 1
  
  #   w <- 1/log(nrow(myBoF)/doc.freq)
  w <- log(doc.freq/nrow(myBoF))
  R <- diag(w)

  Phi_d <- myBoF %*% R
  dimnames(Phi_d) <- dimnames(myBoF)
  Phi_d <- Phi_d[order(rownames(Phi_d)),]
#   Phi_d <- Phi_d[rownames(enriched_Phi_d),colnames(enriched_Phi_d)]
#   Phi_d <- Phi_d[!duplicated(Phi_d), ]

  Phi_d <- Phi_d[,which(!apply(Phi_d,2,FUN = function(x){length(which(x != 0)) <= 1}))]
  Phi_d <- Phi_d[which(!apply(Phi_d,1,FUN = function(x){all(x == 0)})),]
  
#   enriched_Phi_d <- enriched_Phi_d[, colnames(Phi_d)]


  names <- hand_picked_names()
  
  prefix_names <- function(names) {
    rootName = "dendroot"
    for(i in 1:length(names)){
      names[[i]] <- paste(rootName, paste("group", i, sep="-"), names[[i]],sep="/")
      
    }
    
    return(names)
  }

  allNames <- unlist(prefix_names(names))  
  names.dend <- build.dendrogam(allNames)
  
  names <- unlist(names)

  myDist <-  dist(decostand(t(Phi_d), method = "normalize"), method="euclidean")
  m <- as.matrix(myDist)
  m <- m[names, names]
  # print(rownames(m)[c(1,4,13,21)])
  dimnames(m) <- list(allNames, allNames)
  # print(rownames(m)[c(1,4,13,21)])

# r <- svd(Phi_d)
# dimnames_Phi_d <- dimnames(Phi_d)
# 
# Phi_d <- r$u %*% diag(length(c(r$d[1:8], rep(0, length(r$d) - 8)))) %*% t(r$v)
# dimnames(Phi_d) <- dimnames_Phi_d

  myDist <- as.dist(m)

  myEnrichedDist <- dist(decostand(t(enriched_Phi_d), method = "normalize"), method="euclidean")
  m <- as.matrix(myEnrichedDist)
  m <- m[names, names]
  dimnames(m) <- list(allNames, allNames)
  myEnrichedDist <- as.dist(m)


# hclusters <- agnes(x = myEnrichedDist, method = "complete")
# plot(hclusters, which.plots = 2, main = "", sub = "", xlab = "")
#ensure dims of myDist and myEnrichedDist match

#Quantitative analysis

# dendrogram from cluster 1 (single-linkage)
hc1 <- hclust(myDist, method="average")
# hc1 <- hclust(myDist, method="complete")
# plot(hc1)

# dendrogram from cluster 2 (complete-linkage)
hc2 <- hclust(myEnrichedDist, method="average")
# hc2 <- hclust(myEnrichedDist, method="complete")
# plot(hc2)

# correlation
cor(cophenetic(hc1),cophenetic(hc2))

# For a confidence level, use the "Mantel Test" from package vegan.
require(vegan)
mantel(cophenetic(hc1), cophenetic(hc2))

#Qualitative analysis
noc <- 7

require(dendextend)
dend1 <- as.dendrogram(hc1)
dend2 <- as.dendrogram(hc2)


dend1 %>% set("labels_color")  %>%  
  set("branches_k_lty", k=noc)  %>%  set("branches_k_color", k = noc) %>% plot

dend2 %>% set("labels_color")  %>%  
  set("branches_k_lty", k=noc)  %>%  set("branches_k_color", k = noc) %>% plot


tanglegram(dend1 , dend2)

tanglegram(dend1 , dend2, lab.cex = 1.5, lwd = 3, edge.lwd = 4, dLeaf = -0.01,columns_width = c(10,4,10), common_subtrees_color_branches = T, sort=T, 
           cex_main_left = 3, cex_main_right = 10,
#            highlight_distinct_edges = T,
           margin_inner= 12,center = TRUE,
#            dLeaf = -0.1, 
           # k_labels = noc, k_branches = noc, main_left = "Enriched BoF", main_right = "Plain BoF")
k_labels = noc, k_branches = noc, main_left = "Plain BoF", main_right = "Enriched BoF")



#MAKE COMPARISON WITH THE PACKAGE STRUCTURE

# compute tree distance 
treeDistance1 = compute_tree_edit_distance_for_hc(hc1, names.dend$graph)
treeDistance2 = compute_tree_edit_distance_for_hc(hc2, names.dend$graph)

clusters.tree1 <- ape::as.phylo(hc1)
clusters.tree2 <- ape::as.phylo(hc2)
priori.tree <- names.dend$tree

path.difference1 <- phangorn::path.dist(clusters.tree1, priori.tree, check.labels = T)
path.difference2 <- phangorn::path.dist(clusters.tree2, priori.tree, check.labels = T)


}


get.sample.names <- function(names, size=150){
  setwd("~/workspace")
  
  fileName <- paste("benchmark", prname , "BoF", paste(prname, "names.txt", sep="-"), sep="/")
  
  if (file.exists(fileName)){
    names = read.table(fileName, sep="\t", fill=FALSE, strip.white=TRUE)
    
    return(as.vector(names$V1))
  }
  
  names <- refine.names(names, size)
  
  
  write(names, fileName, sep="\t")
  
  
  return(names)
}


manual.names <- function(){
  names <- c("goToBufferEnd", "bufferSwitcher", "paintFoldShape", "applyDockingLayout", "foldPainter", "bufferListeners", "loadBufferSwitcher",   
    "getBufferSet", "getBufferSetManager",   "BSHINIT",  "getBshPrompt" ,  "getFoldPainter",
    "bufferTree","caretLine", "insertionSort", "insureNodesParsed" , "insureParsed" , "parseColor", "parseXML",
    "createRegexpEOLSpanRule", "createRegexpIndentRule" , "getPermissiveTextReader", "getTextReader",                    
    "clearRegister", "getRegister", "getRegisters", "insertRegister"  , "loadRegisters", "register",
    "registerList" , "registerName" , "Registers" , "registerTransferableService", "saveRegisters", "setRegister", "activatePlugin",
    "activatePluginIfNecessary", "addPlugin", "addPluginDockable", "addPluginJAR" ,  "addPluginJARsFromDirectory", "addPluginLocalizationProps",       
    "addPluginProps", "breakPlugin", "checkForObsoletePlugins" , "disablePlugin", "downloadingPluginList" , "getAllPluginEntries",
    "getDependentPlugins", "getNotLoadedPluginJARs" , "getPlugin", "getPluginHome",  "getPluginJAR", "getPluginJARs", "getPluginList",
    "getPlugins", "loadPluginSet", "pluginDetail", "createCustomMenu", "createMacrosMenu", "createMenuItem" , "createPopupMenu", 
    "encodingMenuItems", "getContextMenu", "beginCompoundEdit", "canChangeEditMode" , "getLayoutSize", "isHorizontalLayout", 
    "preferredLayoutSize", "addLayoutComponent", "contentTextArea", "gutterCurrentLineHighlightEnabled", "parseColor" , "highlight", 
    "acceptFile", "getFoldAtLine",  "jEditHome" , "removePluginJAR", "showBuffer", "addPluginLocalizationProps", "downloadXml",
    "getXMLEncoding", "getTextArea", "caretPositionChanged", "logLineCount", 
    "getDependentPlugins", "isStartupDone", "getEndColumn", "getSearchDialog", "addCurrentToHistory",
    "camelCasedWords", "getColumnCount",  "timerIncrementalSearch", "setTokenMarker", "userKeymapFolder", "getLineSegment",
    "getToolTipText", "searchField", "gutterClickActions", "gutterComponents", "gutterCurrentLineHighlightEnabled", "gutterEnabled",  
    "invocationHandler", "gzipURL", "highlightDigits" , "findCompletion", "systemKeymapFolder",  "searchBar",
    "addFoldStyleChooser", "indentLine", "toggleWordWrap", "getResourcePath", "macroHandlers",  "getMainRuleSet", "loadFavorites", "clipText",
    "mappingFeedbackListener", "currentEncoding",   "serviceName",  "getNextVisibleLine",  "scrollpane", "getIgnoreWhitespace",  "deepIndent",
    "boldFont", "escapeRule", "exitOnEOF", "abbreviate" , "addRuleSet" ,"getLogListModel",
    "invokeAction",  "invokeMethod", "invokeObjectMethod", "invokeReadNextChar", "invokeStaticMethod",  "declaringInterpreter",
    "Interpreter", "escapeRule", "getEscapeRule", "getMainRuleSet" , "getRuleSetAtOffset", "resolveExpectedJavaField", "resolveExpectedJavaMethod",
    "resolveImports", "resolveJavaMethod", "httpEnabled", "httpHost"  , "addBrowserListener",   "awtQueue",  "dispatch",
    "runInDispatchThread",  "preview", "previewStatusBar",  "userDir", "userInput", "userKeymapFolder", "keyEventInterceptor","processKeyEvent",
    "translateKeyEvent", "toolbarBox", "captionBox", "openFile" , "newFile" ,  "zipFile", "jEditFileList", "formatFileSize", "registerName",
    "escapeRegexp", "regionMatches", "registersXML", "createRegexpEOLSpanRule", "getCompletions",  "focusedComponent",  
    "cacheDockableWindows", "needFullRepaint", "painter", "paintFoldShape", "replaceAll" ,  "getReplaceFromRemoveInsert",  "replace",
    "performOperationsInAWTThread",  "performOperationsInWorkThread",  "gzipURL", "console", "getNextToken", "matchToken", "TOKEN_TYPES",
    "perspectiveXML","registersXML" )
  
  names <- unique(names)
  
  return(names)
}
"action, box, component, event, button, layout, GUI"

refine.names <- function(names, size){
  
  names <- cleanse.names(names)
  
  patterns <- c("BSH", "buffer", "tar", "textArea", "menu", "parse",            
                 "regular", "regex", "XML", "dispatch", "microstar", "reader", "register", "receive",                
                "plugin", "paint", "caret", "gutter", "layout", "edit")
  
  indices <- c()
  
  for(i in 1:length(patterns)){
    
    temp <- which(unlist(lapply(names, function(x) grepl(pattern= patterns[i] ,x, ignore.case = T ) )))
    
    print(length(temp))
    
    indices <- c(indices, temp )
  }
  
  indices <- unique(indices)
  
  all_names <- names[indices]
  
  half_all_names <- sample(all_names, ceiling(2 * (size/3)))
  
  names <- names[-indices]
  
  other_half_all_names <- sample(names, ceiling(size/3))
  
  return(c(half_all_names, other_half_all_names))
}

hand_picked_names <- function() {
  names <- list( 
    #BSH
    c("BSHINIT", 
    "isJavaBaseAssignable",  #"inNativeCode", 
    "Interpreter", "invoke", "resolveJavaMethod", "isWrapperType"), 
          
 #Regular Expression
 c("escapeRule" , #"getMainRuleSet", 
 "startRegexp", "endRegexp", #"regexp", 
 #"expression", 
 "pattern", "terminateChar", "matchType") , #"MATCH_TYPE_CONTEXT", "MATCH_TYPE_RULE"  , 
 #"getLeadingWhiteSpace", #"match","escapeRegexp", 


#Text Area
 c("caretLine", "autoIndent" , "findMatchingBracket" ,  #"selectToMatchingBracket", #"lineHighlight" ,"highlight", "indentLine" ,
"getLineCount" ,  "caret", #"getCaretPosition", "getCaretLine",
"getSelectedText") , #"moveCaretPosition", 


#XML support
#"parseXML",  
#"Registers" , "registerTransferableService", "saveRegisters", "setRegister",
            
#User Interface
#"preferredLayoutSize" , 
c("processKeyEvent",  "toolbar", "menubar" ,
#"showPopupMenu", 
"needFullRepaint", "focusedComponent", "addDockableWindow"), #, "actionBar", "searchField"
#"editPane", 
#,"OPEN_DIALOG", "SAVE_DIALOG",

#Core
c("settingsDirectory", "queueAWTRunner", "createVFSSession", #"zipFile", "fileVFS", #"createVFSSession", , #"sourceFile" , 
"invokeAction", "buffer", "handleMessage"), 
#"bufferCount", #"bufferLoaded", "bufferOpened",  #"queueAWTRunner", 
# "updateBufferStatus", 
# "getKeyEventInterceptor",  "addBufferListener", #"MESSAGE", "performOperationsInAWTThread",

#plugins
c("activatePlugin"  , #"addPluginJAR", #"PluginOptions", 
# "download", 
"author", "description", "version", "pluginSet", "getPluginJAR"),

#Macros
# "loadMacros",  "runScript", 

#Search & Replace
c("doBackwardSearch", "doForwardSearch", 
"searchField",  "hyperSearch" , "replace" , "replaceSelection")
)

# pd <- c(rep(1, 7), rep(2,3), rep(3, 7), rep(4, 5), rep(5, 9), rep(6, 7), rep(7, 8), rep(8, 2), rep(9, 6))
# pd <- c(rep(1, 7), rep(2,14), rep(3, 13), rep(4, 5), rep(5, 13), rep(6, 16), rep(7, 9), rep(8, 2), rep(9, 6))
# pd <- c(rep(1, 7), rep(2,14), rep(3, 13), rep(4, 5), rep(5, 11), rep(6, 16), rep(7, 9), rep(8, 2), rep(9, 6))
# pd <- c(rep(1, 6), rep(2,13), rep(3, 13), rep(4, 4), rep(5, 6), rep(6, 11), rep(7, 8), rep(8, 6))

# pd <- c(rep(1, 5), rep(2,11), rep(3, 10), rep(4, 6), rep(5, 6), rep(6, 8), rep(7, 4)

# f1.score <- compute.f1(c1, pd) 
# adjustedRI <- compute.AdjRI(c1, pd)
# mojosim <- compute.MoJoSim(c1, pd)
# 
# r1 <- list(mojosim = mojosim, f1.score=f1.score, adjustedRI=adjustedRI)
  
 return(names) 
}


plot.arc.diagram <- function(myBoF){
  
  require(arcdiagram)
  require(igraph)
  
  myBoF <- myBoF[which(!apply(myBoF,1,FUN = function(x){all(x == 0)})),]
  
  
  nIds <- dim(myBoF)[2]
  identifierMatrix <- matrix(0, nrow = nIds, ncol = nIds, dimnames = list(colnames(myBoF), colnames(myBoF)) )
  
  s <- apply(myBoF, 1, function(x) which(x!=0))
  
  for (k in 1:length(s)) {
    names <- s[[k]]
    for (i in 1:length(names))
      for (j in i:length(names))
        if (i != j){
          currentVal <- identifierMatrix[names[i], names[j]]
          print("sdsd")
          identifierMatrix[names[i], names[j]] <- currentVal + 1
        }
  }
  
  identifierMatrix <- fill_lower_diagonal(identifierMatrix)
  g <- igraph::graph.adjacency(identifierMatrix, mode = "undirected", add.colnames = T, add.rownames = T)
  
  edgeList = get.edgelist(g)
  
  set.seed(120)
  arcplot(edgeList, ordering=sample(1:dim(identifierMatrix)[1]), labels=colnames(identifierMatrix),
          lwd.arcs=4*runif(10,.5,2), col.arcs=hsv(runif(9,0.6,0.8),alpha=0.4),
          show.nodes=TRUE, pch.nodes=21, cex.nodes=runif(10,1,3), 
          col.nodes="gray80", bg.nodes="gray90", lwd.nodes=2)
}


cleanse.names <- function(names){
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
  
  indices <- c()
  
  for(i in 1:length(patterns)){
    
    temp <- which(unlist(lapply(names, function(x) grepl(pattern= patterns[i] ,x, ignore.case = T ) )))
    
    indices <- c(indices, temp )
  }
  
  indices <- unique(indices)
  # "BSH_PACKAGE", "textArea", "regex", "buffer" , "getXml", "submenu", "parse"
  
  names <- names[-indices]
  
  return(names)
}

