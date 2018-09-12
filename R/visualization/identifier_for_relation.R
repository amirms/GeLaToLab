getters <-  c(
    "getBooleanProperty",
    
    "getIntegerProperty",
    
    "getSelectedIndex",
    "getMaximum" ,
    
    "getPaper",
    "getValue",
    "getShowHiddenFiles" ,
    

    "getDefault",
    "getStructureHighlightColor",
    "getIgnoreCase",
    "getSearchString",
    "getListener",
    "getFontHandler",
    "getScope"
    
  );

setters <-  c(
  "setBooleanProperty",
  
  "setIntegerProperty",
  
  "setSelectedIndex",
  "setMaximum" ,
  
  "setPaper",
  "setValue",
  "setShowHiddenFiles" ,
  "setDefault",
  "setStructureHighlightColor",
  "setIgnoreCase",
  "setSearchString",
  "setListener",
  "setFontHandler",
  "setScope"
  
  
)

identifiers_for_relations <- c(
  getters,setters
)


identifiersPerClass <- list(
  list("getMaximum", "setMaximum"),
  list("getBorder", "setBorder"),
  list("getVisibleRowCount", "setVisibleRowCount"),
  #list("getDefaultRenderer", "setDefaultRenderer"),
  list("getKeyEventInterceptor", "setKeyEventInterceptor"),
  list("getReshowDelay", "setReshowDelay"),
  list("getFirstLine", "setFirstLine"),
  #list("getMethod", "setMethod"),
  list("getSelectedItem", "setSelectedItem"),
  #list("getStyles", "setStyles"),
  list("getMultipleSelectionColor", "setMultipleSelectionColor"),
  list("getSelectionColor", "setSelectionColor"),
  list("getSelectionPath", "setSelectionPath"),
  #list("getTerminateChar", "setTerminateChar"),
  list("getAutoWrapAround", "setAutoWrapAround"),
  list("getDefault", "setDefault"),
  list("getDefaultMax", "setDefaultMax"),
  list("getDTDHandler", "setDTDHandler"),
  list("getAutoReload", "setAutoReload"),
  #list("getCurrentLineForeground", "setCurrentLineForeground"),
  #list("getColumnWidth", "setColumnWidth"),
  list("getRepeatCount", "setRepeatCount"),
  #list("getFeature", "setFeature"),
  #list("getIntegerProperty", "setIntegerProperty"),
  #list("getMaximumSize", "setMaximumSize"),
  list("getIgnoreCase", "setIgnoreCase")
  
)




#abbreviated, domain-specific vocabulary and encrypted forms.

finGettersAndSetters <- function(names){
getters = names[grepl("^get",names)]
setters = names[grepl("^set",names)]


getterRoot = unique(unlist(strsplit(getters, 'get')))
getterRoot = getterRoot[-1]

results = c()
for(i in 1:length(getterRoot)){
  if (paste("set", getterRoot[i],sep="") %in% setters){
    results <- c(results, paste("get", getterRoot[i],sep=""), paste("set", getterRoot[i],sep=""))
  }
}

return(results)

}


plot = function(){
  library(Rtsne)
  
  identifiers = results[101:120]
  # identifiers = unlist(identifiersPerClass)
  
  ## calling the installed package
  # train<- read.csv(file.choose(),sep=";") ## Choose the train.csv file downloaded from the link above  
  train<- read.csv("C:\\Users\\AmirM\\workspace\\org.graph.cnn\\res\\jedit-5.1.0\\identifier2vec.bern", sep=";");
  
  train <- train[,-101]
  # names <- read.csv(file.choose())
  names <- read.csv("C:\\Users\\AmirM\\workspace\\org.graph.cnn\\data\\jedit-5.1.0\\identifier2id.txt")
  
  
  names <- apply(names,1, function(n) strsplit(n, ";")[[1]][1])
  
  
  rownames(train) <- names

  classes = getClasses(identifiersPerClass, names)
  
  train$label = names
  
  ## Curating the database for analysis with both t-SNE and PCA
  # train$label = unlist(lapply(seq(1:dim(train)[1]), function(x) paste("id", x, sep="_")))
  Labels<-train$label
  train$label<-as.factor(train$label)
  ## for plotting
  colors = rainbow(length(unique(train$label)))
  names(colors) = unique(train$label)
  
  # colors = rainbow(max(classes))
  # names(colors) = unique(classes)
  

    x <- prcomp(train[,-dim(train)[2]], center = T, scale. = T)
    tsne <- list(Y=x$x[identifiers,1:2])

  ## Plotting
  plot(tsne$Y, t='n', main="tsne")
  text(tsne$Y, labels=rownames(tsne$Y))
  # text(tsne$Y, labels=train$label, col=colors[classes])
  
  
  ## Plotting using ggplot2
  library(ggplot2)
  
  data = as.data.frame(tsne$Y)
  # color = Topics)
  pc1 <- ggplot(data, aes(x=PC1, y=PC2)) + theme_bw() + theme(#axis.title.x=element_text("PC1"), axis.title.y=element_text("PC2"),
                                                                      # axis.text.x=element_blank(),
                                                                      # axis.ticks.x=element_blank(),
                                                                      panel.grid.major = element_blank(),
                                                                      panel.grid.minor = element_blank(), 
                                                                      axis.line = element_line(colour = "black"))
  
  # pc1 + geom_point()
  # pc1 + geom_point(shape = 1, size = 4)
  
  pc2 <- pc1 + geom_point(size = 3, color="blue")#, stroke = 2.5)
  library("ggrepel")
  
  # x <- list(c(3,4), c(1,5), c(2,15), c(6,17), c(7,13), c(8,20), c(9,16), c(10,26), c(11,19), c(12,21), c(14, 24), c(18,23 ), c(22,25));
  
  
  g <- c("G1", "G1", "G2", "G2", "G3", "G3", "G4", "G4", "G5", "G5",
         "G6", "G6", "G7", "G7", "G8","G8", "G9", "G9", "G10", "G10",
         "G11", "G11", "G12", "G12", "G13", "G13", "G14", "G14", "G15", "G15",
         "G16", "G16", "G17", "G17") 
  
  
  pc3 <- pc2 + geom_line(aes(group=g), linetype = 2, color="red", size = 1)
  
  (pc4 <- pc3 +     geom_text_repel(aes(label = identifiers),
                                    color = "gray20", force=1))
}

