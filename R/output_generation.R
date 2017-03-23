
models = list(list(type="ContextModel", metrics=list()),
              list(type="String", metrics=list("LCS", "LCU", "Constant")), 
              list(type="SemanticRelatedness", metrics=list("BoF", "IPL", "WP", "LC", "CD")),
              list(type="DG", metrics=list()), 
              list(type="BoF", metrics=list())
              )
PROJECTS = list(list("apache-ant-1.9.3", "Apache Ant"),list("hadoop-0.20.2", "Apache Hadoop"),list("apache-log4j-1.2.17", "Apache Log4j"),
                list("eclipse-jdt-core-3.8","Eclipse JDT Core"), list("jdom-2.0.5","JDOM"), list("jedit-5.1.0","JEdit"),
                list("jfreechart-1.2.0","JFreeChart"), list("jhotdraw-7.0.6","JHotDraw"), list("junit-4.12","JUnit"), list("weka-3.6.11","Weka"))

EVAL_SCORES = c("Bk", "Gamma", "Coph") #The order must be the same as exported text

#Input: model type: one of the "ContextModel", "DG", "BoF", "String", "SemanticRelatedness"

generate.latex.table <- function(model, metrics =c("MoJoSim", "F1", "PMI"), make_comparison=NULL, highlight_best=T){
  all_Types = unlist(lapply(models, function(model) model$type))
  stopifnot(model %in% all_Types)
  
  resultsDirectory = ""
  
  if (model == "ContextModel")
    resultsDirectory = model
  else if ((model == "String") || (model == "BoF") || (model=="SemanticRelatedness"))
    resultsDirectory = "EnrichedBoF"
    
  resultsDirectory = paste("Results", resultsDirectory, sep="/")
  
  model_metrics = unlist(lapply(models, function(m) if (m$type == model) return(m$metrics) ))
  
  allResults = list()
  
  for(i in 1:length(PROJECTS)) {
    setwd("~/workspace")
    project = PROJECTS[[i]]
    
    projectDirectory = project[[1]]
    projectName = project[[2]]
    
    dir.create(file.path(getwd(), paste("benchmark", projectDirectory, resultsDirectory, sep="/")), showWarnings = FALSE)

    setwd(paste("benchmark", projectDirectory, resultsDirectory, sep="/"))
    
    files = list.files(".")
    
    results = list()
    if (length(files) > 0) {
      for (j in 1:length(files)){
        
        scores = read.table(files[j])
        
        indices <- unlist(lapply(model_metrics, function(m) grepl(m, files[j])))

        if (any(indices)) {
          metric <- model_metrics[which(indices)]
           results[[length(results)+1]] <- list(metric=metric, scores=scores)
          
        }
        
      }
    }
    
    allResults[[length(allResults) + 1]] <- list(project=projectName, result = results)
    
  }
  
  tex <- convertResultsToTable(allResults, model_metrics, make_comparison, highlight_best)
  return(tex)
  
}

convertResultsToTable <- function(results, model_metrics, make_comparison, highlight_best=1){
  txts=""
  for(i in 1:length(results)){
    prname = results[[i]]$project
    result = results[[i]]$result
    
    noOfScores <- length(EVAL_SCORES)
    maxScores <- rep(0, length(EVAL_SCORES))
    
    line <- prname
    
    allScores = list()
    
    baselineScore = NULL
    
    for(j in 1:length(model_metrics)) {
      scores <- findScoresForMetric(result, model_metrics[j])
      currentScores <- rep(NULL, len)
      
      if (!is.null(scores)){
        
        print("Scores")
        print(scores)
        scores <- toString(scores$V1)
        
        splitted_scores <- as.numeric(unlist(strsplit(scores, "&")))
        
        for(i in 1:noOfScores) {
        
          currentScore <- splitted_scores[i]
        
          if (maxScores[i] > currentScore) {
            maxScores[i] <- currentScore
          }
          
          currentScores[i] <- currentScore
        }
        
        #Is this the index for baseline comparison
        if (!is.null(highlight_best) && j == highlight_best) {
          baselineScore = currentScores
        }
      
      }
      
      allScores[[length(allScores) + 1]] <- list(currentScores)
    }
    
    # print("printing length of MJ_F1")
    # print(MJ_F1)
    
    for (s in 1:length(allScores)) {
      score = allScores[[s]]
      

      for(i in m:length(scores)){
        cell = ""
        
        if(!is.null(score[i])){
          
          cell <- printScoreCell(score[i], maxScores[i], baselineScore[i], s==highlight_best)
          
        }
        line <- paste(line, cell, sep="&")
      }
      
      line <- paste(line, cell, sep="&")
      
    }
    
    line <- paste(line, "\\\\ \\hline \n", sep="")
    txts <- paste(txts, line)
  }
  

  return(txts)
  
}

printScoreCell <- function(score, maxScore=NULL, baselineScore=NULL, isSameAsBaseline=F) {
  
  cell=""
  
  if (!is.null(maxScore) && score == maxScore) {
    cell <- paste(cell, "\\cellcolor[gray]{0.8}{", sep="")
  }
  
  cell <- paste(cell, score, sep="")
  if ((!is.null(baselineScore)) && (!isSameAsBaseline)){
    x = (((score - baselineScore) / baselineScore) * 100)
    percentage <- format(round(x, 2), nsmall = 2)
    if (x> 0) {
      percentage <- paste("+", percentage, "%", sep="")
    }
    cell <- paste(cell, "(", percentage, ")", sep="")
    
  }
  
  
  if (!is.null(maxScore) && score == maxScore) {
    cell <- paste(cell, "}", sep="")
  }
  
  return(cell)
  
}

findScoresForMetric <- function(results, metric) {

  if(length(results) > 0){
    for (i in 1:length(results)){
      r <- results[[i]]
      if (r$metric == metric)
        return(r$scores)
    }
  }
      
  return(NULL)
}

output_semantic_relatedness_results <- function() {
  model = "SemanticRelatedness"
  tex <- generate.latex.table(model = model, highlight_best = 1)
  print_table_results(tex, paste(model,".txt", sep=""))
}

print_table_results <- function(result, fileName){
  setwd("~/workspace")
  
  write(result, file = paste("benchmark", "Results", fileName, sep="/"))
}