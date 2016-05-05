LCS = function(s1, s2)
{
  
  s1 <- tolower(s1)
  s2 <- tolower(s2)

  cm = matrix(0, nrow=nchar(s1)+1, ncol=nchar(s2)+1)
  
  for (i in 2:(nchar(s1)+1))
    cm[i,1] = 0
  for (j in 2:(nchar(s2)+1))
    cm[1, j] = 0
  
  for (i in 2:(nchar(s1)+1))
    for (j in 2:(nchar(s2)+1))
    {
      if (substr(s1, i-1,i-1) == substr(s2, j-1,j-1))
        cm[i, j] = cm[i - 1, j - 1] + 1
      else
      {
        cm[i, j] = max(cm[i - 1, j], cm[i, j - 1])
        
      }
      
    }
  
  cm[nchar(s1)+1, nchar(s2)+1]; 
  
#   cm
}
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