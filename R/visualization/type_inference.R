infer_types <- function() {
  library(igraph)
  library(readtext)
  
  g <- igraph::make_empty_graph();
  
  #read files
  directory = "C:\\workspace\\org.graph.cnn\\data\\jedit-5.1.0"
  files = c("test.ee.txt", "train.ee.txt", "valid.ee.txt");
  
  d = lapply(files, function(f) readtext(paste(directory, f, sep="\\")))
  z <- unlist(lapply(d, function(x) strsplit(x$text, '\n') [[1]]))
  z <- gsub("#", "$", z)
  
  all <- read.table(text = z,  sep=";",
                    col.names=c("head", "tail", "relation"))
  
  
  isAssignedTo <- all[which(all$relation== "isAssignedTo"),]
  
  freshIdentifiers = isAssignedTo$tail[startsWith(as.vector(isAssignedTo$tail), "_freshI")]  
  indices <- which(all$head %in% freshIdentifiers | all$tail %in% freshIdentifiers)
  
  freshIsAssignedToData <- all[indices,]
  unqiueFreshIsAssignedToIdentifiers <- unique(freshIdentifiers)
  
  
  
  isInstanceOf <- all[which(all$relation== "isInstanceOf"),]
  
  freshIdentifiers = isInstanceOf$head[startsWith(as.vector(isInstanceOf$head), "_freshI")]  
  indices <- which(all$head %in% freshIdentifiers  |  all$tail  %in%  freshIdentifiers )
  
  freshIsInstanceOfData <- all[indices,]
  
  unqiueFreshIdentifiers <- unique(freshIdentifiers)
  
  
  identifier2vec <- read.csv("C:\\workspace\\org.graph.cnn\\res\\jedit-5.1.0\\identifier2vec.transr.ee.bern", sep=";", header = F);
  identifier2vec <- data.matrix(identifier2vec[,-101])
  
  identifier2id <- read.csv("C:\\workspace\\org.graph.cnn\\data\\jedit-5.1.0\\identifier2id.txt", sep=";", header = F);

  getIdentifier2Id <- function(identifierName) {
    identfierId <- as.vector(identifier2id[which(identifier2id$V1 == identifierName),]$V2) + 1  
    return(identfierId)
  }
  
  getIdentifierName <- function(identifierId) {
    identifierName <- as.vector(identifier2id[which(identifier2id$V2 == (identifierId -1) ),]$V1)
    
    return(identifierName)
  }
  
  
  head_relation_identifier_matrices <- read.csv("C:\\workspace\\org.graph.cnn\\res\\jedit-5.1.0\\headRelationIdentifierMatrix.transr.ee.bern", sep=";", header = F);
  head_relation_identifier_matrices <- head_relation_identifier_matrices[,-51]
  
  tail_relation_identifier_matrices <- read.csv("C:\\workspace\\org.graph.cnn\\res\\jedit-5.1.0\\tailRelationIdentifierMatrix.transr.ee.bern", sep=";", header = F);
  tail_relation_identifier_matrices <- tail_relation_identifier_matrices[,-51]
  
  id2relation <- read.csv("C:\\workspace\\org.graph.cnn\\data\\jedit-5.1.0\\relation2id.txt", sep=";", header = F)
  isOfTypeId <- id2relation[which(id2relation$V1 == "isOfType"),]$V2 +1
  
  
  isOfTypeRelationIndices <- ((isOfTypeId -1) * 100 + 1):(isOfTypeId*100);
  
  head_isOfType_relation_matrix <- data.matrix(head_relation_identifier_matrices[isOfTypeRelationIndices,]);
  tail_isOfType_relation_matrix <- data.matrix(tail_relation_identifier_matrices[isOfTypeRelationIndices,]);
  
  relation2vec <- read.csv("C:\\workspace\\org.graph.cnn\\res\\jedit-5.1.0\\relation2vec.transr.ee.bern", sep=";", header = F);
  relation2vec <- relation2vec[,-51]
  isOfTypeRelation2vec <- as.numeric(relation2vec[isOfTypeId,])
  
  ##################INSTACE_OF####################
  # freshes <- c("_freshIdentifier711", "_freshIdentifier815", "_freshIdentifier730", "_freshIdentifier555", "_freshIdentifier328")
  
  freshNames <- c("_freshIdentifier555", "_freshIdentifier730","_freshIdentifier815",
                              "_freshIdentifier587", "_freshIdentifier827", "_freshIdentifier829", "_freshIdentifier497", "_freshIdentifier788" )
  
  freshIndices <- unlist(lapply(freshNames, function(freshIdentifier) getIdentifier2Id(freshIdentifier)))
 
  fresh_head_isOfType <- identifier2vec[freshIndices, ] %*% head_isOfType_relation_matrix;
  rownames(fresh_head_isOfType) <- freshNames
  # 
  # types <- c("org.gjt.sp.jedit.textarea.RectParams", "java.lang.String", "org.gjt.sp.jedit.bsh.ParseException",
  #            "org.gjt.sp.jedit.textarea.LineCharacterBreaker", "org.gjt.sp.jedit.textarea.Rect", 
  #            "java.util.Vector", "java.util.Hashtable", "java.util.List", "javax.swing.JTextField", "javax.swing.JLabel", "javax.swing.JMenu",
  #            "javax.swing.JPanel", "org.gjt.sp.jedit.OptionGroup", "java.awt.Insets",  "org.gjt.sp.util.IntegerArray") 
  
  typeNames <- c("java.lang.Integer", "int", "java.lang.String", "org.gjt.sp.jedit.textarea.Rect",
             "org.gjt.sp.jedit.textarea.LineCharacterBreaker", "double" )
  
  typesIndices <- unlist(lapply(typeNames, function(typeIdentifier) getIdentifier2Id(typeIdentifier)))
  
  types_tail_isOfType <- identifier2vec[typesIndices, ] %*% tail_isOfType_relation_matrix;
  
  rownames(types_tail_isOfType) <- typeNames
  
  
  actual_fresh_types = fresh_head_isOfType + isOfTypeRelation2vec
  rownames(actual_fresh_types) = paste("type of ", rownames(actual_fresh_types), sep="")
  
  
  # combined <- rbind(fresh_head_isOfType, actual_fresh_types, types_tail_isOfType)
  # combined <- rbind(fresh_head_isOfType, ( identifier2vec %*% head_isOfType_relation_matrix), types_tail_isOfType)
  
  # a <- prcomp(combined, center=F, scale.=F)$x
  # names = c(rownames(a))

  library(ggplot2)
  library(ggrepel)
  
  pc1 <- ggplot(as.data.frame(combined), aes(x=PC1, y=PC2)) + theme_bw() + theme(#axis.title.x=element_text("PC1"), axis.title.y=element_text("PC2"),
        # axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(), axis.title = element_text(size = 14),
       axis.line = element_line(colour = "black"),
       text = element_text(size = 14))
  
  Identifiers <- factor(colors)
  pc2 <- pc1 + geom_point(aes(color = Identifiers), size = 3)#, color=c("red", "blue", "green"))#, stroke = 2.5)
  
  
  
  # g <- c("G1", "G2", "G3", "G4", "G5", "G1", "G2", "G3", "G4", "G5") 
  
  
  pc3 <- pc2 + geom_line(aes(group=g), linetype = 2, color="red", size = 1)
  
  pc3 <- pc3 + geom_line(aes(group=g2), linetype = 2, color="black", size = 1)
  
  
  pc4 <- pc3 + geom_text_repel(aes(label = rownames(combined)), size=5,
                                    color = "gray20", force=1)
  
}

get_top_types <- function(unqiueFreshIdentifiers, freshIsInstanceOfData) {
  
  for(i in 500:700) {
    freshIdentifier <- as.vector(unqiueFreshIdentifiers[i])
    freshIndices <- which(freshIsInstanceOfData$head == freshIdentifier  |  freshIsInstanceOfData$tail ==  freshIdentifier )
    
    if (length(freshIndices) <= 5) {
      # print(freshIdentifier)
      # print(freshIsInstanceOfData[freshIndices,])
      
      print(i)
      
      freshId <- getIdentifier2Id(freshIdentifier)  # it is +1 id
      
      hs <- as.numeric(identifier2vec[freshId,])
      
      scores <- apply(identifier2vec, 1, function(ts) compute_struct_energy_function(hs, head_isOfType_relation_matrix, isOfTypeRelation2vec, as.numeric(ts), tail_isOfType_relation_matrix))
      
      a <- sort(scores, decreasing = T, index=T)
      
      topIdentifiers <- unlist(lapply(a$ix[1:10], function(index) getIdentifierName(index)))
      
      print(freshIsInstanceOfData[freshIndices,])
      
      print("Top 3 Predicted Types")
      print(topIdentifiers)
      
    }
  }
}


compute_struct_energy_function <- function(hs, Mrh, r, ts, Mrt) {
  
  # -Math.abs(e2_vec[ii] - e1_vec[ii] - relation_vec[rel][ii]);
  return(sum(-1 * abs(ts %*% Mrt- r - hs %*% Mrh )))
  
  # return(hs %*% Mrh + r - ts %*% Mrt )
}

compute_struct_energy_function2 <- function(hs, Mrh, r, ts, Mrt) {
  
  # -Math.abs(e2_vec[ii] - e1_vec[ii] - relation_vec[rel][ii]);
  return(sum(abs(hs %*% Mrh + r - ts %*% Mrt  )))
  
  # return(hs %*% Mrh + r - ts %*% Mrt )
}