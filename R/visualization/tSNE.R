visualize_tsne=function(identifiersPerClass=list(), applyTSNE=TRUE){
  
  library(Rtsne)
  
  identifiers = unlist(identifiersPerClass);

## calling the installed package
# train<- read.csv(file.choose(),sep=";") ## Choose the train.csv file downloaded from the link above  
train<- read.csv("C:\\Users\\AmirM\\workspace\\org.graph.cnn\\res\\jedit-5.1.0\\identifier2vec.bern", sep=";");

train <- train[,-101]
# names <- read.csv(file.choose())
names <- read.csv("C:\\Users\\AmirM\\workspace\\org.graph.cnn\\data\\jedit-5.1.0\\identifier2id.txt")


names <- apply(names,1, function(n) strsplit(n, ";")[[1]][1])


rownames(train) <- names

# unknown_words <- read.csv("C:\\Users\\AmirM\\workspace\\org.graph.cnn\\res\\jedit-5.1.0\\unmatchedWords.bern")
# unknown_words <- apply(unknown_words,1, function(n) strsplit(n, ";")[[1]][1])

# word2vec <- read.csv("C:\\Users\\AmirM\\workspace\\org.graph.cnn\\res\\jedit-5.1.0\\word2vec.bern", sep=";");
# word2vec <- word2vec[,-301]
# 
# words <- read.csv("C:\\Users\\AmirM\\workspace\\org.graph.cnn\\data\\jedit-5.1.0\\word2id.txt")
# words <- apply(words,1, function(n) strsplit(n, ";")[[1]][1])
# 
# rownames(word2vec) <- words

if (length(identifiers) > 0){
  indices <- which(names%in%identifiers)

  # rnd = sample(length(1:length(names)),500 , TRUE)
  # indices <- unique(c(rnd, indices))
  train <- train[indices,]
  names <- names[indices]
}

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

if (applyTSNE){
## Executing the algorithm on curated data
tsne <- Rtsne(train[,-dim(train)[2]], dims = 2, perplexity=3, verbose=TRUE, max_iter = 15000)
# exeTimeTsne<- system.time(Rtsne(train[,-1], dims = 2, perplexity=10, verbose=TRUE, max_iter = 500))
} else{ #apply PCS
  x <- prcomp(train[,-dim(train)[2]], center = T, scale. = T)
  tsne <- list(Y=x$x[,1:2])
}
## Plotting
plot(tsne$Y, t='n', main="tsne")
text(tsne$Y, labels=train$label, col=colors[train$label])
# text(tsne$Y, labels=train$label, col=colors[classes])


## Plotting using ggplot2
library(ggplot2)

Topics <- unlist(lapply(classes, function(c) {
 if (c == 1){
   return("BeanShell")
 }  else if (c==2){
   return("Regular Expression")
 }  else if (c==3){
   return("Text Area")
 } else if (c==4){
   return("User Interface")
 } else if (c==5){
   return("Core")
 } else if (c==6){
   return("Plugins")
 } else { # c==7
   return("Search & Replace")
   
 }
  
}))

data = as.data.frame(tsne$Y)
# color = Topics)
pc1 <- ggplot(data, aes(x=data$V1, y=data$V2)) + theme_bw() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
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


g <- c("G1", "G3", "G2", "G2", "G1", "G4", "G5", "G6", "G7", "G8",
                  "G9", "G10", "G5", "G11", "G3","G7", "G4", "G12", "G9", "G6",
                  "G10", "G13", "G12", "G11", "G13", "G8") 


pc3 <- pc2 + geom_line(aes(group=g), linetype = 2, color="red", size = 1)

(pc4 <- pc3 +     geom_text_repel(aes(label = names),
              color = "gray20", force=1))
#,              data = subset(dat, Country %in% pointsToLabel)))
}


getClasses <- function(identifiers=list(), names=c()){
  if (length(identifiers) == 0){
    return(classes=seq(1:length(names)))
  }
  
  classes = rep(0, length(names))
  
  for( i in 1:length(identifiers)){
    x = identifiers[[i]];
    
    indices = which(names %in% x);
    
    classes[indices] = i;
  }
  
  return(classes)
}


# find anatonym identifier names: getter and setter identifiers


# 
findClosestWords = function(word2vec, words){
  K = cos.sim(t(word2vec))
  
  dimnames(K) = list(rownames(word2vec), rownames(word2vec));
  
  result = list();
  for(i in 1:length(words)){
    word = words[i];
    
    if(word %in% rownames(K)) {
    result[[i]] = sort(K[word,], decreasing = T)[1:4];
    }
  }
  
  
  return(result);
}
