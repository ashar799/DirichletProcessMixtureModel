### This function takes the posterior parameters AND predicts the time for the new points
#### The fundamental assumption is that EACH NEW TEST POINT IS CONDITIONALLY INDEPENDENT on the OTHER POINTS
#### We predict value of one point at a time
### The final output is Time for the new samples, ONE AT A TIME
### Just need to collect samples for S

predictchineseAFTtime = function(Y.new){
  
  
  c.new.list <- list(0)
  ## The number of posterior samples
  That.new <- time.new 
  post.time  = matrix(NA,nrow = nrow(Y.new), ncol = Nps)
  print("GOING THROUGH MCMC Samples")
  pb <- txtProgressBar(min = 1, max = Nps , style = 3)
  
  for (count in 1:Nps){
    
    ## Assign the parameters to the posterior sample
    ctemp <- c.list[[count]]
    mu <- mu.list[[count]]
    S <- S
    beta0 <- beta0.list[[count]]
    betahat  <- betahat.list[[count]] 
    sigma2  <- sigma2.list[[count]] 
    lambda2 <- lambda2.list[[count]] 
    tau2 <- tau2.list[[count]] 
    Ytemp <- Y
    Ytemp.scaled <- matrix(NA, nrow = N, ncol = D)
    g <- table(factor(ctemp, levels = 1:K))
    activeclass <- which(g!=0)
    ## The table function helps converting the data point specific indicator variables to class specific indicator variables
    kminus <- length(activeclass)
    ## Two Auxilary Variables
    ## The name of the auxilary variables are taken to be one and two more than the maximum value in the already active cluster set
    activeclass <- append(activeclass, max(activeclass)+1)
    activeclass <- append(activeclass, max(activeclass)+1)
    active <- activeclass 
    ### Assigning values to parameters 
    
    source('priordraw.R')
    priortwo <- NA
    priortwo <- priordraw(beta, W, epsilon, ro, r, si,N,D, sig2.dat)
    mu[active[kminus+1],1:D]  <- priortwo$mu  
    S[active[kminus+1],1:D,1:D]  <- priortwo$Sigma[1:D,1:D]  
    beta0[active[kminus+1]] <- priortwo$beta0 
    sigma2[active[kminus+1]] <- priortwo$sigma2
    betahat[active[kminus+1],1:D] <- priortwo$betahat 
    lambda2[active[kminus+1]] <- priortwo$lambda2 
    tau2[active[kminus+1], 1:D] <- priortwo$tau2
    
    
    source('priordraw.R')
    priorthree <- NA
    priorthree <- priordraw(beta, W, epsilon, ro, r, si,N,D, sig2.dat)
    mu[active[kminus+2],1:D]  <- priorthree$mu  
    S[active[kminus+2],1:D,1:D]  <- priorthree$Sigma[1:D,1:D]  
    beta0[active[kminus+2]] <- priorthree$beta0 
    sigma2[active[kminus+2]] <- priorthree$sigma2
    betahat[active[kminus+2],1:D] <- priorthree$betahat 
    lambda2[active[kminus+2]] <- priorthree$lambda2 
    tau2[active[kminus+2], 1:D] <- priorthree$tau2
    
    #######################################################
    ctemp.new = c(0)
    
    
    ## This can't be parallelized !!!!!
    for(l in 1:N.new)  {
      
      
      posteriortime <- matrix(NA, nrow = length(active), ncol = 1)
      posteriornormtime <- matrix(NA, nrow = length(active), ncol = 1)
      
      Y.new.sc <- matrix(0, nrow = N.new, ncol =D)
      
      ## Calculating the Expectations and also the normalization constant for the Expectation
      for (j in 1:kminus) {
        
        clust <- which(ctemp == active[j])
        
        obj.t <- scale(Y[clust,1:D], center = TRUE, scale = TRUE)
        
        for ( h in 1:D){
          Y.new.sc[l,h] <- (Y.new[l,h] - attr(obj.t,"scaled:center")[h])/(attr(obj.t,"scaled:scale")[h])
        }
        
        posteriortime[j] <- g[active[j]]  * dMVN(as.vector(t(Y.new[l,1:D])), mean = mu[active[j],1:D], Q = S[active[j],1:D,1:D]) *  ( beta0[active[j]] + betahat[active[j],1:D] %*% as.vector(t(Y.new.sc[l,1:D])) )
        
        posteriornormtime[j] <- g[active[j]]  * dMVN(as.vector(t(Y.new[l,1:D])), mean = mu[active[j],1:D], Q = S[active[j],1:D,1:D]) 
        }
      
      res <- try(dMVN(as.vector(t(Y[l,1:D])), mean = mu[active[kminus+1],1:D], Q= S[active[kminus+1],1:D,1:D]), silent=TRUE)
      if (class(res) == "try-error"){
        posteriortime[kminus+1] <- 0
        posteriornormtime[kminus+1] <- 0
      } else{
        posteriortime[kminus+1] <- alpha * dMVN(as.vector(t(Y[l,1:D])), mean = mu[active[kminus+1],1:D], Q= S[active[kminus+1],1:D,1:D]) *  (beta0[active[kminus+1]] + betahat[active[kminus+1],1:D] %*% as.vector(t(Ytemp[l,1:D])) )
        posteriornormtime[kminus+1] <-   alpha * dMVN(as.vector(t(Y[l,1:D])), mean = mu[active[kminus+1],1:D], Q= S[active[kminus+1],1:D,1:D]) 
      }
      
      res2 <- try(dMVN(as.vector(t(Y[l,1:D])), mean = mu[active[kminus+2],1:D], Q= S[active[kminus+2],1:D,1:D]), silent=TRUE)
      if (class(res) == "try-error"){
        posteriortime[kminus+2] <- 0
        posteriornormtime[kminus+2] <- 0
        
      } else{
        posteriortime[kminus+2] <- alpha * dMVN(as.vector(t(Y[l,1:D])), mean = mu[active[kminus+2],1:D], Q= S[active[kminus+2],1:D,1:D]) *  (beta0[active[kminus+2]] + betahat[active[kminus+2],1:D] %*% as.vector(t(Ytemp[l,1:D])) )
        posteriornormtime[kminus+2] <-   alpha * dMVN(as.vector(t(Y[l,1:D])), mean = mu[active[kminus+2],1:D], Q= S[active[kminus+2],1:D,1:D])      
      }
      
      
      ## Calculating the actual time for one Set of Parameter for one data point
       post.time[l,count] <-  sum(posteriortime)/ sum(posteriornormtime)
     
      
    }
    
    Sys.sleep(0.1)
    setTxtProgressBar(pb, count)
    
    
  }
  
  #### To calculate average values over MCMC samples
  post.time.avg <<- apply(post.time[,1:20],1,mean)
  
 
  
}
