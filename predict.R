predict = function(c, Y, That, Time, beta0, betahat, sigma2, Ystar ) {
  
  
  Ystarscaled <- scale(Ystar, center = TRUE, scale = TRUE )
  numclust <- table(factor(c, levels = 1:K))
  activeclass<- which(numclust!=0)
  
  
  
  weight  <- matrix(0, nrow = nrow(Ystar), ncol = length(activeclass))
 
  
for (l in 1:nrow(Ystar))
  for ( i in 1:length(activeclass)) {
    weight[l,i] <-  dMVN(as.vector(t(Ystar[l,1:D])), mean = mu[activeclass[i],1:D], Q = S[activeclass[i],1:D,1:D]) 
  }
  
  
  tstar <- rep(0, nrow(Ystar))
  ct <-    rep(0, nrow(Ystar))
  
  for ( h in 1:nrow(Ystar)){  
      ct[h] <- sample(x = activeclass,1, prob = weight[h,], replace = TRUE)
       tstar[h] <- rnorm(1, mean = beta0[activeclass[ct]] + betahat[activeclass[ct],1:D ] %*% Ystarscaled[h,1:D], sd = sqrt(sigma2[activeclass[ct]]))
    }
    
    
  
  
  list('tstar' = tstar,'ct'= ct) 
  
  
}
