kmeansBlasso = function(Y,That, F,K, beta, W, epsilon, ro, r, si, N, D, sig2.dat, c, mu, S, beta0, betahat, sigma2, lambda2, tau2 ) {
  
  source('priordraw.R')
  G <- F
  k.data <- kmeans(Y,F,nstart =10)
  c <- k.data$cluster
  
  #### Under special cases
  ###### c <- c.true
  
  
### We either initialize the cluster assignments from the time data or from molecular data based on AIC criterion
#   source('kmeansAIC.R')
#   aic.data <- flexmix(Y ~ 1, k = 2, model = FLXMCmvnorm())
#   
#   ### For the time data points we use first a LASSO type model to fit and th
#   object <- cv.glmnet(Y,That)
#   ind.good <- which(coef(object)!=0)
#   Y.good <- Y[,ind.good]
#   ### Fit a FlexMix model on this subsetted Y
#   aic.time <- flexmix(That ~ Y.good, k =2)
#   
#   if( AIC(aic.data) < AIC(aic.time)){
#     c <- clusters(aic.data)
#   } else {
#     c <- clusters(aic.time)
#   }
#   
  
  
  prior.numclust <- table(factor(c, levels = 1:K))
  prior.activeclass<- which(prior.numclust!=0)
  
  ### The means are set using the k-means
  for ( i in 1:length(prior.activeclass)){
    mu[prior.activeclass[i],1:D] <-  k.data$centers[i,1:D] 
    
    S[prior.activeclass[i],1:D,1:D] <-  priordraw(beta, W, epsilon, ro, r, si,N,D, sig2.dat)$Sigma
    
    lclust <- which(c == prior.activeclass[i])
    
    reg.blas <- 0
    
    sum <- c(0)
    
    coeff <- 0
    
    Ytemp <-  matrix(NA, nrow = length(lclust), ncol = D)
    
    Ytemp <- scale(Y[lclust,1:D], center = TRUE, scale = TRUE)
    
    
    ### Part where I use the MONOMVN PACKAGE
    
    Ttemp <- as.vector(That[lclust])
    
    ntemp <- length(lclust)
      
    reg.blas <- blasso(Ytemp, Ttemp, T = 300,thin =  50, RJ = TRUE, mprior = 0.0 ,normalize = TRUE, verb = 0)
      
    sum <- summary(reg.blas, burnin= 100)
      
    ## Selecting those features which are relevant
    
    coeff <- unlist(lapply(strsplit(sum$coef[3,], split = ":"), function(x) as.numeric(unlist(x)[2])))
      
    beta0[prior.activeclass[1]] <- coeff[1]
      
    indexplusone <- D+1
      
    ind <- 2:indexplusone
      
    betahat[prior.activeclass[i], ] <- coeff[ind]
      
    ta <- unlist(lapply(strsplit(sum$tau2i[3,], split = ":"), function(x) as.numeric(unlist(x)[2])))
      
    tau2[prior.activeclass[i],] <- ta
      
    sigma2[prior.activeclass[i]] <- sum$s2[3]
      
    lambda2[prior.activeclass[i]] <- sum$lambda2[3]
    
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