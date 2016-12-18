# (LCU) Longest common substring
length_longest_common_substring <- function(String1, String2) {
  
  s1 <- unlist(strsplit(String1,split="")) # rozbicie stringu na tablice znakow
  s2 <- unlist(strsplit(String2,split=""))
  
  num <- matrix(0,nchar(String1),nchar(String2) )  	
  maxlen <- 0
  
  for (i in 1:nchar(String1)) {
    
    for (j in 1:nchar(String2)) {
      
      if (s1[i] == s2[j]) {
        if ((i==1) || (j==1)) { 
          num[i,j] <- 1
        } 
        else {
          num[i,j] <- 1+num[i-1,j-1]
        }
        if (num[i,j] > maxlen) {
          maxlen <- num[i,j]
        }
      }
    }
  }
  
  maxlen			
}

#Input: two strings str1, str2
#Output: length_longest_common_substring(str1, str2)^2/(len(str1)*len(str2))

compute_normalized_LCU <-function(str1, str2){
  llcu <- length_longest_common_substring(str1, str2)
  
  llcu^2/(nchar(str1) * nchar(str2))
}


#Input: two strings str1, str2
#Output: length of longest common subsequence LCS(str1, str2)^2/(len(str1)*len(str2))
compute_normalized_LCS <- function(str1, str2){
  require(qualV)
  
  llcs <- LCS(strsplit(str1, "")[[1]], strsplit(str2, "")[[1]])$LLCS
  
  llcs^2/(nchar(str1) * nchar(str2))
  
}


normalized_LCU_kernel <- function(names){
  names <- tolower(names)
  len <- length(names)
  
  m = matrix(0, nrow=len, ncol=len)
  
  for (i in 1:len)
    for(j in i:len)
      if (i != j)
        m[i,j] <- compute_normalized_LCU(names[i], names[j])
  
  m <- fill_lower_diagonal(m)
  
  dimnames(m) <- list(names, names)
  
  m
}


normalized_LCS_kernel <- function(names){
  
  names <- tolower(names)
  len <- length(names)
  
  m = matrix(0, nrow=len, ncol=len)
  
  for (i in 1:len)
    for(j in i:len)
      if (i != j)
        m[i,j] <- compute_normalized_LCS(names[i], names[j])
  
  m <- fill_lower_diagonal(m)
  
  dimnames(m) <- list(names, names)
  
  m
}

