kmeansinit = function(Y,time, F,K, beta, W, epsilon, ro, r, si, N, D, sig2.dat) {
  
  
  c <- 0
  mu = matrix(data = NA, nrow = K, ncol = D)
  S = array(data = NA, dim =c(K,D,D))
  lambda2 <- numeric(K)
  tau2 = matrix(data = NA, nrow = K, ncol = D)
  betahat = matrix(data = NA, nrow = K, ncol = D)
  sigma2 <- rep(NA, K)
  beta0 <- rep(NA, K)
  That <-  numeric(N)
  
  source('priordraw.R')
  G <- F
  k.data <- kmeans(Y,G)
  c <- k.data$cluster
  
  
  prior.numclust <- table(factor(c, levels = 1:K))
  prior.activeclass<- which(prior.numclust!=0)
  
  ### The means are set using the k-means
  for ( i in 1:length(prior.activeclass)){
    mu[prior.activeclass[i],1:D] <-  k.data$centers[i,1:D] 
 
    lclust <- which(c == prior.activeclass[i])
    if (length(lclust) > D){
      Ysc <- scale(Y[lclust,1:D], center = TRUE, scale =TRUE)
      lm.data <- lm(time[lclust] ~ Ysc)
      sigma2[prior.activeclass[i]] <-  var(lm.data$residuals)
      ind <- D +1
      betahat[prior.activeclass[i],] <-  lm.data$coefficients[2:ind]
      beta0[prior.activeclass[i]] <-  mean(time[lclust])
    }
    else{
      prior <- priordraw(beta, W, epsilon, ro, r, si,N,D, sig2.dat)
      sigma2[prior.activeclass[i]] <- prior$sigma2
      betahat[prior.activeclass[i],1:D] <- prior$betahat 
      beta0[prior.activeclass[i]] <-  mean(time[lclust])
      }
    
    lambda2[prior.activeclass[i]] <- priordraw(beta, W, epsilon, ro, r, si,N,D, sig2.dat)$lambda2
    tau2[prior.activeclass[i], 1:D] <-  priordraw(beta, W, epsilon, ro, r, si,N,D, sig2.dat)$tau2
    S[prior.activeclass[i],1:D,1:D] <-  priordraw(beta, W, epsilon, ro, r, si,N,D, sig2.dat)$Sigma
  }
  
 ## Deleting those values which are no longer relevant
 g <- table(factor(c, levels = 1:K))
 inactive <- which(g==0)
 
 for ( i in 1:length(inactive)){
   mu[inactive[i],1:D]  <- NA 
   S[inactive[i],1:D,1:D]  <- NA  
   beta0[inactive[i]] <- NA 
   sigma2[inactive[i]] <- NA
   betahat[inactive[i],1:D] <- NA 
   lambda2[inactive[i]] <- NA
   tau2[inactive[i], 1:D] <- NA
 }
 
 
  
list('c' = c, 'mu'=mu, 'beta0'=beta0, 'betahat'= betahat, 'sigma2' =sigma2, 'lambda2' = lambda2, 'tau2'= tau2, 'S' =S)  
  
}