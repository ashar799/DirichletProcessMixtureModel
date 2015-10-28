#### THIS function predicts the Class of the new Data Points
#### It is Based on the PredictChineseAFT function

predictCLASS = function(Y.new, time.new){
  
  N.new <<- nrow(Y.new)
  c.new.list <- list(0)
  ## The number of posterior samples
  Nps = as.integer(iter/ iter.thin)
  That.new <- time.new 
  
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
      
      
      posterior <- matrix(NA, nrow = length(active), ncol = 1)
      Y.new.sc <- matrix(0, nrow = N.new, ncol =D)
      
      ## Calculating the probabalities for drawing the value of c_i from the active classes
      for (j in 1:kminus) {
        
        clust <- which(ctemp == active[j])
        
        obj.t <- scale(Y[clust,1:D], center = TRUE, scale = TRUE)
        
        for ( h in 1:D){
          Y.new.sc[l,h] <- (Y.new[l,h] - attr(obj.t,"scaled:center")[h])/(attr(obj.t,"scaled:scale")[h])
        }
        
        posterior[j] <- g[active[j]] /(N-1+alpha) * dMVN(as.vector(t(Y.new[l,1:D])), mean = mu[active[j],1:D], Q = S[active[j],1:D,1:D]) *  dnorm(x = That.new[l], mean = beta0[active[j]] + betahat[active[j],1:D] %*% as.vector(t(Y.new.sc[l,1:D])), sd = sqrt(sigma2[active[j]]) )
      }
      
      res <- try(dMVN(as.vector(t(Y[l,1:D])), mean = mu[active[kminus+1],1:D], Q= S[active[kminus+1],1:D,1:D]), silent=TRUE)
      if (class(res) == "try-error"){
        posterior[kminus+1] <- 0
      } else{
        posterior[kminus+1] <- (0.5 * alpha) /(N-1+alpha) * dMVN(as.vector(t(Y[l,1:D])), mean = mu[active[kminus+1],1:D], Q= S[active[kminus+1],1:D,1:D]) *  dnorm(x = That.new[l], mean = beta0[active[kminus+1]] + betahat[active[kminus+1],1:D] %*% as.vector(t(Ytemp[l,1:D])), sd = sqrt(sigma2[active[kminus+1]]) )
      }
      
      res2 <- try(dMVN(as.vector(t(Y[l,1:D])), mean = mu[active[kminus+2],1:D], Q= S[active[kminus+2],1:D,1:D]), silent=TRUE)
      if (class(res) == "try-error"){
        posterior[kminus+2] <- 0
      } else{
        posterior[kminus+2] <- (0.5 * alpha) /(N-1+alpha) * dMVN(as.vector(t(Y[l,1:D])), mean = mu[active[kminus+2],1:D], Q = S[active[kminus+2],1:D,1:D]) *  dnorm(x = That.new[l], mean = beta0[active[kminus+2]] + betahat[active[kminus+2],1:D] %*% as.vector(t(Ytemp[l,1:D])), sd = sqrt(sigma2[active[kminus+2]]) )
      }
      
      
      #     posterior[kminus+1] <- (0.5 * alpha) /(N-1+alpha) * dMVN(as.vector(t(Y[l,1:D])), mean = mu[active[kminus+1],1:D], Q= S[active[kminus+1],1:D,1:D]) *  dnorm(x = That[l], mean = beta0[active[kminus+1]] + betahat[active[kminus+1],1:D] %*% as.vector(t(Ytemp[l,1:D])), sd = sqrt(sigma2[active[kminus+1]]) )
      #     
      #     posterior[kminus+2] <- (0.5 * alpha) /(N-1+alpha) * dMVN(as.vector(t(Y[l,1:D])), mean = mu[active[kminus+2],1:D], Q = S[active[kminus+2],1:D,1:D]) *  dnorm(x = That[l], mean = beta0[active[kminus+2]] + betahat[active[kminus+2],1:D] %*% as.vector(t(Ytemp[l,1:D])), sd = sqrt(sigma2[active[kminus+2]]) )
      #     
      
      ## Calculating the normalization constant for probabilities
      normalization <- sum(posterior) 
      
      if (normalization < 1e-200 || normalization ==Inf){
        ctemp.new[l] <- sample(active, 1, prob = rep(1,length(active)), replace = TRUE)
      } else {  
        ctemp.new[l] <- sample(active, 1, prob= posterior, replace = TRUE)
      }
      
    }
    
    c.new.list[[count]] <- ctemp.new
    Sys.sleep(0.1)
    setTxtProgressBar(pb, count)
    
  }
  
  #### To calculate the posterior probabilities
  posteriorprob <- matrix(0, nrow = N.new, ncol = kminus+ 1)
  rownames(posteriorprob) <- rownames(Y.new)
  for ( i in 1:N.new){
    temp.c <- c(0)
    for ( j in 1:Nps){
      temp.c[j] <- c.new.list[[j]][i] 
    }
    for ( v in 1:kminus){
      posteriorprob[i,v] <- length(which(temp.c ==v))
    }
    posteriorprob[i,kminus+1] <-  length(which(temp.c ==kminus+1)) + length(which(temp.c ==kminus+2))
  }
  
  posteriorprob <- posteriorprob/Nps
  
  posteriorprob <<- posteriorprob
  
}
