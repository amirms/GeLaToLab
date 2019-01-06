nmf <- function(prname = "jedit-5.1.0", noOfTopics = 8){
  require(NMF)
  
  BoW <- load_BoW(prname)
  # apply tf-idf mechanism and eliminate features that are lower than some threshold,
  # then remove those features from BoW
  x <- apply_tf_idf(BoW)
  dimnames(x) <- dimnames(BoW)
  BoW <- x
  
  result <- nmf(BoW, noOfTopics, "frobenius")
  
  H <- coef(result)
  
  x <- apply(H, 1, function(x) names(sort(x, decreasing=T)[1:10]))
  apply(x, 1, function(z) paste(z, collapse="&"))
}