posteriorhyper = function(c, Y, mu, S, epsilon, W, beta, ro ) {
  
  numclust <- table(factor(c, levels = 1:K))
  activeclust <- which(numclust!=0)
  nactive <- length(activeclust)
  InvCov <- solve(cov(Y))
  meandata <- apply(Y, 2, mean )
  meandata <-  as.matrix(meandata)
#   This was to be used incase the lo likelihood was needed
#   logS <- c(rep(0, N))
#   for ( i in 1:N ) {
#   logS[i] <- ( (D/2) * ((ro + nj)/(ro + nj +1)) - ((D/2) * log(3.14)) + lgamma((beta+nj+1)/2) - lgamma((beta+nj+1-D)/2) + ((beta+ nj)/2 * log(abs(det(Wst))) - ((beta + nj +1)/2)* log( abs(det(Wst +  ((ro + nj)/(ro + nj+ 1)) * (Y[i,1:D] - epsilonstar)%o% Y[i,1:D] - epsilonstar)                           )))
# 
#   }
  
# Update the Epsilon paramter
   sum.precision <- matrix(0, nrow = D, ncol =D)
   sum.mean.precision <- matrix(0, nrow = D, ncol =1)
   for ( z in 1:nactive) {
   sum.precision <- sum.precision + ro * S[activeclust[z],1:D, 1:D]
   sum.mean.precision <-  sum.mean.precision + ro* S[activeclust[z],1:D, 1:D] %*% as.matrix(mu[activeclust[z],1:D]) 
  }
   precision.epsilon <- InvCov + sum.precision
   mean.epsilon <- solve(precision.epsilon) %*% ( InvCov %*% meandata + sum.mean.precision) 

  epsilon <- mvrnorm(n=1, mu = as.vector(mean.epsilon), Sigma = solve(precision.epsilon)) 
  
# Update the ro paramter
  
   sum.ro <- 0
   for ( z in 1:nactive) {
   sum.ro <- sum.ro + t(as.matrix(mu[activeclust[z],1:D]- epsilon)) %*% S[activeclust[z],1:D, 1:D] %*% as.matrix(mu[activeclust[z],1:D]- epsilon)
   }
  ro <- rgamma(1, shape = (nactive/2 + 0.25 ), scale = (as.numeric(sum.ro) +0.5)^-1)

# Update W the Wishart parameter
 
  sum.w <- sum.precision <- matrix(0, nrow = D, ncol =D)
  for ( z in 1:nactive) {
  sum.w <- sum.w + beta * S[activeclust[z],1:D, 1:D]
  }

  
res <- try(rWishart(n = 1, df = beta * nactive + D, Sigma =  solve(D * InvCov + beta * sum.w )), silent=TRUE)
if (class(res) == "try-error"){
  W = W
} else{
  W = rWishart(n = 1, df = beta * nactive + D, Sigma =  solve(D * InvCov + beta * sum.w ))
}

  list('epsilon' = epsilon,'W' = W , 'ro' = ro) 
}
