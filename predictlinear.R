predictlinear = function(c, Y, That, Time, beta0, betahat, sigma2 ) {
  
  Ytemp <- matrix(NA, nrow = N, ncol = D)
  numclust <- table(factor(c, levels = 1:K))
  activeclass<- which(numclust!=0)
  
  prediction <- c(0)
  
  
  for ( i in 1:length(activeclass)) {
    
    clust <- which(c == activeclass[i])
    
    Ytemp[clust,1:D] <- scale(Y[clust,1:D], center = TRUE, scale = TRUE)
  }
  
  
  for ( h in 1:N){  
    prediction[h]<-  beta0[c[h]] + betahat[c[h],1:D ] %*% Ytemp[h,1:D] 
  }


  list('predictlinear' = prediction) 
  
  
}
