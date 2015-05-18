loglikelihood = function(c,Y,mu,S,alpha,That, beta0, betahat, sigma2, lambda2, tau2, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.dat) {
  
  disclass <- table(factor(c, levels = 1:K))
  activeclass <- which(disclass!=0)
  
  loglikelihood <- rep(0, length(activeclass))
  
  
  for (j in 1:length(activeclass)) {
    clust <- which(c==activeclass[j])
    
    Ytemp <-  matrix(NA, nrow = length(clust), ncol = D)
    
    if (length(clust)==1){
      Ytemp <- matrix(0, nrow =1, ncol =D)
      
    } else {
      Ytemp <- scale(Y[clust,1:D], center = TRUE, scale = TRUE)
    }
   
   loglikelihood[j] <- 0
    for ( l in 1:length(clust)) {
    if (Time[clust[l],2]==1){
      loglikelihood[j] <- loglikelihood[j] +  log(dMVN(x = as.vector(t(Y[clust[l],1:D])), mean = mu[activeclass[j],1:D],  Q = S[activeclass[j],1:D,1:D])) + log(dnorm(x = That[clust[l]], mean = beta0[activeclass[j]] + betahat[activeclass[j],1:D] %*% as.vector(t(Ytemp[l,1:D])), sd = sqrt(sigma2[activeclass[j]]) ))
    } else{
      loglikelihood[j] <- loglikelihood[j] + log(dMVN(x = as.vector(t(Y[clust[l],1:D])), mean = mu[activeclass[j],1:D], Q = S[activeclass[j],1:D,1:D])) + log(dtruncnorm(x = That[clust[l]], a = Time[clust[l],1], b = max(Time[,1]), mean = beta0[activeclass[j]] + betahat[activeclass[j],1:D] %*% as.vector(t(Ytemp[l,1:D])), sd = sqrt(sigma2[activeclass[j]]) ))
    }
    }
    
    
  }
  
  return(sum(loglikelihood))
  
  
}