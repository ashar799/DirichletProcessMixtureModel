calR2error = function(c, Time, time.predicted ) {
  numclust <- table(factor(c, levels = 1:K))
  activeclust <- which(numclust!=0)
  
  cis <- c(0)
  
  for (j in 1:length(activeclust)) {
    
  }
  
  
  list('cindex' = cis)
}
