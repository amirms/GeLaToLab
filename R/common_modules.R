compute_common_modules <- function(prname){

  modules_names1 <- compute_common_modules_identifiers(prname)
  modules_names2 <- compute_common_modules_types(prname)
  modules_names3 <- compute_common_modules_identifier_types(prname)
  modules_names4 <- compute_common_modules_dependency_graphs(prname)

  common_names <- intersect_all(modules_names1, modules_names2, modules_names3, modules_names4)
  old_common_names = c()
  
  while(!all(common_names %in% old_common_names)) {
    old_common_names <- common_names
    modules_names1 <- compute_common_modules_identifiers(prname, old_common_names)
    modules_names2 <- compute_common_modules_types(prname, old_common_names)
    modules_names3 <- compute_common_modules_identifier_types(prname, old_common_names)
    modules_names4 <- compute_common_modules_dependency_graphs(prname, old_common_names)
    
    common_names <- intersect_all(modules_names1, modules_names2, modules_names3, modules_names4)
  }
  
  #WRITE
  setwd("~/workspace")
  write(common_names, file = paste("benchmark", prname ,"MODULE_NAMES.txt", sep="/"))
  
  return(common_names)
}

intersect_all <- function(a,b,...){
  Reduce(intersect, list(a,b,...))
}

compute_common_modules_dependency_graphs <- function(prname, module_names =c(), dirname="weka"){
  setwd("~/workspace")
  setwd(paste("benchmark", prname , "DG", sep="/"))
  
  dependencies = list()
  
  pattern <- "*.csv"
  filenames = list.files(path = dirname, pattern = pattern, all.files = FALSE,
                         full.names = TRUE, recursive = TRUE,
                         ignore.case = TRUE, include.dirs = TRUE, no.. = TRUE)
  
  filenames <- unlist(lapply(filenames, function(fn) substr(fn, 1, regexpr(".csv", fn) -1)))
  
  
  print("The length of filenames before eliminating small packages")
  print(length(filenames))
  if (length(module_names) == 0) {
    module_names <- filenames
  } 
  
  filenames <- eliminate_small_packages(module_names)
  print("The length of filenames after eliminating small packages")
  print(length(filenames))
  
  return(filenames)
}

compute_common_modules_identifier_types <- function(prname, module_names =c()){
  setwd("~/workspace")
  
  # Read the authoritative decomposition
  # decomposition <- read.csv(paste("benchmark", prname ,"decomposition.csv", sep="/"), sep=",",  header = TRUE)
  # priori.decomp <- decomposition$x
  # names(priori.decomp) <- decomposition$X
  # priori.decomp <- normalizeVector(priori.decomp)
  
  #Bag of Features
  myBoF_data <- load_BoF(prname, c(T,T)) 
  myBoF <- myBoF_data$myBoF
  myBoF <- merge_names_by_lower_case(myBoF, 2)
  
  print("The dimension of myBoF before eliminating small packages")
  print(dim(myBoF))
  if (length(module_names) == 0) {
    module_names <- rownames(myBoF)
  } 
  
  myBoF <- myBoF[eliminate_small_packages(module_names),]
  print("The dimension of myBoF after eliminating small packages")
  print(dim(myBoF))
  
  #Remove unused identifier_type names 
  myBoF <- myBoF[,which(!apply(myBoF,2,FUN = function(x){all(x == 0)}))]
  #Remove empty classes/interfaces
  myBoF <- myBoF[which(!apply(myBoF,1,FUN = function(x){all(x == 0)})),]
  
  return(rownames(myBoF))
}

compute_common_modules_types <- function(prname, module_names =c()){
  setwd("~/workspace")

  #Bag of Features
  myBoF <- load_BoF(prname, c(F,T)) 
  myBoF <- merge_names_by_lower_case(myBoF, 2)
  
  print("The dimension of myBoF before eliminating small packages")
  print(dim(myBoF))
  if (length(module_names) == 0) {
    module_names <- rownames(myBoF)
  } 
  
  myBoF <- myBoF[eliminate_small_packages(module_names),]
  print("The dimension of myBoF after eliminating small packages")
  print(dim(myBoF))
  
  #Remove unknown type
  unknownIdx <- which(colnames(myBoF) == "Unknown")
  myBoF <- myBoF[,-unknownIdx]
  
  #Remove unused type names 
  myBoF <- myBoF[,which(!apply(myBoF,2,FUN = function(x){all(x == 0)}))]
  myBoF <- myBoF[which(!apply(myBoF,1,FUN = function(x){all(x == 0)})),]
  
  return(rownames(myBoF))
}

compute_common_modules_identifiers <- function(prname, module_names =c()){
  
  setwd("~/workspace")
  
  #Bag of Features
  myBoF <- load_BoF(prname, c(T,F)) 
  myBoF <- merge_names_by_lower_case(myBoF, 2)
  
  print("The dimension of myBoF before eliminating small packages")
  print(dim(myBoF))
  if (length(module_names) == 0) {
    module_names <- rownames(myBoF)
  } 
  
  myBoF <- myBoF[eliminate_small_packages(module_names),]
  print("The dimension of myBoF after eliminating small packages")
  print(dim(myBoF))
  
  #Remove unused identifiernames 
  myBoF <- myBoF[,which(!apply(myBoF,2,FUN = function(x){all(x == 0)}))]
  
  #Filter out names shorter than 5
  identifierNames <- colnames(myBoF)
  identifierNames <- identifierNames[which(unlist(lapply(identifierNames, nchar))>=5)]
  myBoF <- myBoF[,identifierNames]
  
  #Remove empty classes/interfaces
  myBoF <- myBoF[which(!apply(myBoF,1,FUN = function(x){all(x == 0)})),]
  
  return(rownames(myBoF))
}