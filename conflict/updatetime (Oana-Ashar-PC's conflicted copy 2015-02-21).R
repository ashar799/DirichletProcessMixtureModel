updatetime = function(c, Y, Time,That, beta0, betahat, sigma2 ) {
  
  
  Ytemp <- matrix(NA, nrow = N, ncol = D)
  numclust <- table(factor(c, levels = 1:K))
  activeclass<- which(numclust!=0)
  
  
  ## Updating the Estimated Survival Time
  for ( h in 1:N) {
    
    if (Time[h,2]==0) {
      
      ## The columns of Y have to be centred
      clust <- which(c == c[h])
      
      if (length(clust)==1){
        for ( k in 1:D){
          Ytemp[h,k] <- 0
        }
        
      } else { tempmatrix <- Y[clust,1:D]
      for ( k in 1:D){
        Ytemp[h,k] <- Y[h,k] - mean(tempmatrix[,k])
      }
      }
      
      
      
      
      That[h]<- rtruncnorm(1, a = Time[h,1], b = Inf, mean = beta0[c[h]] + betahat[c[h],1:D ] %*% Ytemp[h,1:D] , sd = sqrt(sigma2[c[h]]) )
      
    } else {
      That[h] <- Time[h,1] }
  }
  
  list('time' = That) 
  
  
}
