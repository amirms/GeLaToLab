
ALL_MODELS = list(list(type="ContextModel", metrics=list("DIFF_STRING", "Enriched_Identifier_Type_BoF", "Enriched_Type_BoF")),
              list(type="String", metrics=list("LCS", "LCU", "CONSTANT")), 
              list(type="SemanticRelatedness", metrics=list("TYPE_IPL", "TYPE_WP", "TYPE_LC", "TYPE_CD")),
              list(type="DG", metrics=list("TYPE_BoF", "Enriched_Identifier_Type_BoF", "DG")), 
              list(type="BoF", metrics=list("BoF", "TYPE_STRING_CD", "DIFF_STRING"))
              )
PROJECTS = list(list("apache-ant-1.9.3", "Apache Ant"),list("hadoop-0.20.2", "Apache Hadoop"),list("apache-log4j-1.2.17", "Apache Log4j"),
                list("eclipse-jdt-core-3.8","Eclipse JDT Core"), list("jdom-2.0.5","JDOM"), list("jedit-5.1.0","JEdit"),
                list("jfreechart-1.2.0","JFreeChart"), list("jhotdraw-7.0.6","JHotDraw"), list("junit-4.12","JUnit"), list("weka-3.6.11","Weka"))

EVAL_SCORES_INDICES = list(list(4, F, F), list(6, F, F))  # c("Bk", "PATH") The order must be the same as exported text

#Input: model type: one of the "ContextModel", "DG", "BoF", "String", "SemanticRelatedness"

generate.latex.table <- function(model, eval_scores_indices =EVAL_SCORES_INDICES, make_comparison=NULL, highlight_best=T){
  all_Types = unlist(lapply(ALL_MODELS, function(model) model$type))
  stopifnot(model %in% all_Types)
  
  resultsDirectory = ""
  
  if ((model == "ContextModel") || (model == "DG"))
    resultsDirectory = "ContextModel"
  else if ((model == "String") || (model == "BoF") || (model=="SemanticRelatedness"))
    resultsDirectory = "EnrichedBoF"
    
  resultsDirectory = paste("Results", resultsDirectory, sep="/")
  
  model_metrics = unlist(lapply(ALL_MODELS, function(m) if (m$type == model) return(m$metrics) ))
  
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
  
  tex <- convertResultsToTable(allResults, model_metrics, eval_scores_indices, make_comparison, highlight_best)
  return(tex)
  
}

# input: all results per project, and result per metric 
# all model-metrics to include: list("BoF", "IPL", "WP", "LC", "CD")
# eval_scores_indices: the indices corresponding to the scores of interest
# make_comparison: the index of the metric used for baseline for X, Y (v%)
# highlight_best: highlight the metric with the best result
#
convertResultsToTable <- function(results, model_metrics, eval_scores_indices, make_comparison=1, highlight_best=T){
  txts=""
  for(i in 1:length(results)){
    prname = results[[i]]$project
    result = results[[i]]$result
    
    print(paste("Project: ", prname))
    
    noOfScores <- length(eval_scores_indices)
    maxScores <- rep(0, noOfScores)
    
    eval_indices <- unlist(lapply(eval_scores_indices, function(idx) idx[1]))
    normalize_indices <- unlist(lapply(eval_scores_indices, function(idx) idx[2]))
    is_max <- unlist(lapply(eval_scores_indices, function(idx) idx[3]))
    
    # minimum value for normalized indices is -1
    maxScores[which(normalize_indices)] <- -1
    maxScores[which(!is_max)] <- Inf
    
    line <- prname
    
    allScores = list()
    
    baselineScore = NULL
    
    for(j in 1:length(model_metrics)) {
      scores <- findScoresForMetric(result, model_metrics[j])
      currentScores <- rep(NULL, noOfScores)
      
      if (!is.null(scores)){
        
        print("Scores")
        print(scores)
        scores <- toString(scores$V1)
        
        splitted_scores <- as.numeric(unlist(strsplit(scores, "&")))
        
        for(k in 1:noOfScores) {
          
          index <- eval_indices[k]
          currentScore <- splitted_scores[index]
          
          if (!is.na(currentScore)){
            if (is_max[k]) {
              maxScores[k] <- max(maxScores[k], currentScore)
            } else {
              maxScores[k] <- min(maxScores[k], currentScore)
            }
            currentScores[k] <- currentScore
          }
        }
        
        #Is this the index for baseline comparison
        if (!is.null(make_comparison) && j == make_comparison) {
          baselineScore = currentScores
        }
      }
      
      allScores[[length(allScores) + 1]] <- list(currentScores)
    }
    
    # print("printing length of MJ_F1")
    # print(MJ_F1)
    
    for (s in 1:length(allScores)) {
      # line <- paste(line, "&", sep="")
      
      scores = unlist(allScores[[s]])
      print(scores)
      

      for(j in 1:length(scores)){
        cell = ""
        
        if(!is.null(scores[j])){

          cell <- printScoreCell(scores[j], maxScores[j], baselineScore[j], s==make_comparison, normalize_indices[j], highlight_best)
          
        }
        line <- paste(line, cell, sep="&")
      }
      
    }
    
    line <- paste(line, "\\\\ \\hline \n", sep="")
    txts <- paste(txts, line)
  }
  

  return(txts)
  
}

printScoreCell <- function(score, maxScore=NULL, baselineScore=NULL, isSameAsBaseline=F, normalize=T, highlight_best=T) {
  
  cell=""
  
  print(score)
  print(maxScore)
  
  isMaxScore <- !is.null(maxScore) && score == maxScore
  
  if (isMaxScore && highlight_best) {
    cell <- paste(cell, "\\cellcolor[gray]{0.8}{", sep="")
  }

  if (score > 100)
    roundedScore <- format(round(score, 0), nsmall = 0)
  else
    roundedScore <- score
  
  cell <- paste(cell, "$", roundedScore, "$", sep="")
  if ((!is.null(baselineScore)) && (!isSameAsBaseline)){
    
    if (normalize){
      score <- score + 1
      baselineScore <- baselineScore + 1
    }
    
    x = (((score - baselineScore) / baselineScore) * 100)
    percentage <- format(round(x, 2), nsmall = 2)
    
    if (x> 0) {
      percentage <- paste("+", percentage, "\\%", sep="")
    } else {
      percentage <- paste(percentage, "\\%", sep="")
    }
    cell <- paste(cell, " $(", percentage, ")$", sep="")
  }
  
  if (isMaxScore && highlight_best) {
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
  tex <- generate.latex.table(model = model, make_comparison = NULL, highlight_best = T)
  print_table_results(tex, paste(model,".txt", sep=""))
}

output_string_results <- function() {
  model = "String"
  tex <- generate.latex.table(model = model, make_comparison = NULL, highlight_best = T)
  print_table_results(tex, paste(model,".txt", sep=""))
}

output_context_vector_results <- function() {
  model = "ContextModel"
  tex <- generate.latex.table(model = model, make_comparison = NULL, highlight_best = T)
  print_table_results(tex, paste(model,".txt", sep=""))
}

output_relatedness_diffusion_comparison_results <- function() {
  model = "BoF"
  tex <- generate.latex.table(model = model, make_comparison = 1, highlight_best = T)
  print_table_results(tex, paste(model,".txt", sep=""))
}

output_boit_dg_comparison_results <- function() {
  model = "DG"
  tex <- generate.latex.table(model = model, make_comparison = 1, highlight_best = T)
  print_table_results(tex, paste(model,".txt", sep=""))
}

print_table_results <- function(result, fileName){
  setwd("~/workspace")
  
  write(result, file = paste("benchmark", "Results", fileName, sep="/"))
}