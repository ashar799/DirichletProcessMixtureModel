posttime = function(c, That, lambda2,tau2,sigma2,beta0, betahat, Y, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.data ) {
  
  numclust <- table(factor(c, levels = 1:K))
  activeclass<- which(numclust!=0)
  
  for (j in 1:length(activeclass)) {
    
    ## A Temporary matrix that needs to store the standardized regressors
    
    clust <- which(c==activeclass[j])
    
    Ytemp <-  matrix(NA, nrow = length(clust), ncol = D)
    
    if (length(clust)==1){
     Ytemp <- matrix(0, nrow =1, ncol =D)
      
    } else {
      Ytemp <- scale(Y[clust,1:D], center = TRUE, scale =TRUE)
             }
    
    tempvector <- as.vector(That[clust])
    tempmean <- mean(tempvector)
    tmpscl <- scale(tempvector, center = TRUE, scale =FALSE)
    tempmatrix <- Ytemp
    tempnumber <- length(tempvector)
    
    
    tempD <- matrix( 0, nrow = D, ncol =D)
    for ( i in 1:D ) {
      tempD[i,i] <- tau2[activeclass[j],i]
    }
    
    
    
    
    
    
    ## For updating the sparsity prior
    lambda2[activeclass[j]] <- rgamma(1, shape = r+D, rate = si + tr(tempD) )
    
    #For updating tau2
    
    for ( h in 1:D)  {
      tau2[activeclass[j], h] <- (rinv.gaussian(1,mu= sqrt(lambda2[activeclass[j]] * sigma2[activeclass[j]]/ (betahat[activeclass[j],h])^2), lambda = lambda2[activeclass[j]]))^-1
    } 
    
    #For updating sigma2
    ## For updating the sigma2 parameter we need temporary matrices
    
    tempprod <- NA
    
    tempscalesigma1 <- as.vector(tmpscl - Ytemp %*% betahat[activeclass[j], ])
    
    tempprod <- tempscalesigma1 %*% tempscalesigma1
    
    tempscalesigma2 <- NA
    
    tempscalesigma2 <- t(betahat[activeclass[j], ] %*% solve(tempD) %*% betahat[activeclass[j], ] )
    
    
    sigma2[activeclass[j]] <- rinvgamma(1, shape = 1+ 0.5 * (tempnumber +D -1), scale = 1 + (0.5* (tempprod + tempscalesigma2 )) )
    ## This is because the error of the model may make it computationally infeasible
    
    
    ## For updating Betahat we need some matrices
    tempD <- matrix( 0, nrow = D, ncol =D)
    for ( i in 1:D ) {
      tempD[i,i] <- tau2[activeclass[j],i]
    }
    
    tempA <-   matrix(NA, nrow = D, ncol = D)

    tempA <- t(Ytemp) %*% Ytemp + solve(tempD)
    
    
    betahat[activeclass[j],] <- mvrnorm(1, mu = solve(tempA) %*% t(tempmatrix) %*% tmpscl, Sigma=  sigma2[activeclass[j]] * solve(tempA))
    
    
    beta0[activeclass[j]] <- rnorm(1, mean = tempmean, sd= sqrt(sigma2[activeclass[j]]/tempnumber))
}







list('beta0' = beta0,'sigma2' = sigma2, 'betahat' = betahat, 'lambda2' = lambda2, 'tau2' =  tau2 )
}
