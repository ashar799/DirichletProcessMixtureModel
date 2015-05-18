posteriorAFTonlyGMM = function(c,Y,mu,S,alpha, K, epsilon, W, beta, ro,r, si,D ,N, sig2.dat) {
  
  
  Ytemp <- matrix(NA, nrow = N, ncol = D)
  ctemp <- c
  
  ## This can't be parallelized !!!!!
  for(l in 1:N)  {
    
    temp <- ctemp[l]
    cminus <- ctemp
    cminus[l] <- NA
    ## The table function helps converting the data point specific indicator variables to class specific indicator variables
    g <- table(factor(cminus, levels = 1:K))
    active <- which(g!=0)
    
    
    
    kminus <- length(active)
    ## Two Auxilary Variables
    ## The name of the auxilary variables are taken to be one and two more than the maximum value in the already active cluster set
    active <- append(active, max(active)+1)
    active <- append(active, max(active)+1)
    
    
    
    ## If the observation was singelton (i.e no other point was associated with it then we assign to kminus +1 parameter)
    if(length(which(cminus==temp))==0)  
    {
      ## The kminus+1 parameter gets the value of the temporary variable
      ctemp[l] <- active[kminus+1]
      mu[active[kminus+1],1:D] <- mu[temp,1:D]
      S[active[kminus]+1,1:D,1:D] <- S[temp,1:D,1:D]
     
      ## Also the second auxilary variable should be drawn from the prior distribution  
      
      source('priordraw.R')
      priorone <- NA
      priorone <- priordraw(beta, W, epsilon, ro, r, si,N,D, sig2.dat)
      mu[active[kminus+2],1:D]  <- priorone$mu  
      S[active[kminus+2],1:D,1:D]  <- priorone$Sigma  
      
    } else {
      
      ## We have to deal with centred matrices
      clust <- which(ctemp == temp)
      
      source('priordraw.R')
      priortwo <- NA
      priortwo <- priordraw(beta, W, epsilon, ro, r, si,N,D, sig2.dat)
      mu[active[kminus+1],1:D]  <- priortwo$mu  
      S[active[kminus+1],1:D,1:D]  <- priortwo$Sigma[1:D,1:D]  
      
      
      source('priordraw.R')
      priorthree <- NA
      priorthree <- priordraw(beta, W, epsilon, ro, r, si,N,D, sig2.dat)
      mu[active[kminus+2],1:D]  <- priorthree$mu  
      S[active[kminus+2],1:D,1:D]  <- priorthree$Sigma[1:D,1:D]  
      
    }
    
    
    
    #######################################################
    
    
    posterior <- matrix(NA, nrow = length(active), ncol = 1)
    
    
    
    
    ## Calculating the probabalities for drawing the value of c_i from the active classes
    for (j in 1:kminus) {
      
      posterior[j] <- g[active[j]] /(N-1+alpha) * dMVN(as.vector(t(Y[l,1:D])), mean = mu[active[j],1:D], Q = S[active[j],1:D,1:D]) 
    }
    
    
    posterior[kminus+1] <- (0.5 * alpha) /(N-1+alpha) * dMVN(as.vector(t(Y[l,1:D])), mean = mu[active[kminus+1],1:D], Q= S[active[kminus+1],1:D,1:D]) 
    
    posterior[kminus+2] <- (0.5 * alpha) /(N-1+alpha) * dMVN(as.vector(t(Y[l,1:D])), mean = mu[active[kminus+2],1:D], Q = S[active[kminus+2],1:D,1:D]) 
    
    
    ## Calculating the normalization constant for probabilities
    normalization <- sum(posterior) 
    
    if (normalization < 1e-200){
      ctemp[l] <- sample(active, 1, prob = rep(1,length(active)), replace = TRUE)
    } else {  
      ctemp[l] <- sample(active, 1, prob= posterior, replace = TRUE)
    }
    
  }
  
  c <- ctemp
  ## Delete those observations that are not associcated with no data point
  g <- table(factor(c, levels = 1:K))
  inactive <- which(g==0)
  
  for ( i in 1:length(inactive)){
    mu[inactive[i],1:D]  <- NA 
    S[inactive[i],1:D,1:D]  <- NA  
   
  }
  
  
  
  list('indicator' = c,'mean' = mu,'precision' = S)
  
}
