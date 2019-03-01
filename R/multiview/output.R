PROJECTS = list(list("apache-ant-1.9.3", "Apache Ant"),list("hadoop-0.20.2", "Apache Hadoop"),list("apache-log4j-1.2.17", "Apache Log4j"),
                list("eclipse-jdt-core-3.8","Eclipse JDT Core"), list("jdom-2.0.5","JDOM"), list("jedit-5.1.0","JEdit"),
                list("jfreechart-1.2.0","JFreeChart"), list("jhotdraw-7.0.6","JHotDraw"), list("junit-4.12","JUnit"), list("weka-3.6.11","Weka"))

#Input: model type: one of the "ContextModel", "DG", "BoF", "String", "SemanticRelatedness"

generate.multiview.clustering.table <- function() {
  resultsDirectory = "MULTIVIEW/Results"
  filename = "All_EVAL.txt"
  NoOfFields = 9
  generate.multiview.table(resultsDirectory, filename, NoOfFields)
}

generate.multiview.cfg.recommendation.table <- function() {
  resultsDirectory = "MULTIVIEW/Recommender/Results"
  filename = "CFG_PRED.txt"
  NoOfFields = 7
  generate.multiview.table(resultsDirectory, filename, NoOfFields)
}

generate.multiview.freq.recommendation.table <- function() {
  resultsDirectory = "MULTIVIEW/Recommender/Results"
  filename = "FREQ_PRED.txt"
  NoOfFields = 7
  generate.multiview.table(resultsDirectory, filename, NoOfFields)
}

generate.multiview.lex.recommendation.table <- function() {
  resultsDirectory = "MULTIVIEW/Recommender/Results"
  filename = "LEX_PRED.txt"
  NoOfFields = 7
  generate.multiview.table(resultsDirectory, filename, NoOfFields)
}

generate.multiview.table <- function(resultsDirectory, filename, noOfFields){

  line=""
  txt =""
  for(i in 1:length(PROJECTS)) {
    setwd("~/workspace")
    project = PROJECTS[[i]]
    
    projectDirectory = project[[1]]
    projectName = project[[2]]
    
    path = paste("benchmark", projectDirectory, resultsDirectory, sep="/")
    destfile = paste(path, filename, sep="/")

    print(destfile)
    
    line=projectName
    if (file.exists(destfile)) {

      scores <- read.table(destfile)
      line <- paste(line, paste(scores$V1, scores$V2, scores$V3, scores$V4, scores$V5, scores$V6, scores$V7, scores$V8) , sep="&")
    } else {
      line <- paste(line, rep("&", noOfFields-1),sep="&")[1]
    }
    
    line <- paste(line, "\\\\ \\hline \n", sep="")
    txt <- paste(txt, line)
  }
  
  setwd("~/workspace")
  dir.create(file.path(getwd(), paste("benchmark", "Multiview/Results", sep="/")), showWarnings = FALSE)
  write(txt, file = paste("benchmark/Multiview/Results", filename, sep="/"))

  return(txt)
}