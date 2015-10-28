rm(list =ls())
 D =2
 K = 50
 ro = 1.6
 beta =3.3

epsilon.real <- c(6, 60)
W.real <- matrix(c(4,-5,-5, 81), nrow = 2, ncol =2)


## Initialization of the parameters for Gaussian Mixture
mu = matrix(data = NA, nrow = K, ncol = D)
S = array(data = NA, dim =c(K,D,D))



for ( j in 1:K){
  S[j,1:D,1:D] <- rWishart(1, beta, solve((beta*W.real)))
  mu[j,1:D] <- mvrnorm(n=1, mu = epsilon.real, Sigma = solve(ro*S[j,1:D,1:D]))
}



## Initialize the desired variables
W = array(data = NA, dim =c(D,D))
epsilon <- matrix(data = NA, nrow = 1, ncol = D)

## BURN-IN period

for (o in 1:50) {
  
  ######################### THE HYPERPARAMETERS OF THE GMM #################################  
  
  # Update the Epsilon paramter
  sum.precision <- matrix(0, nrow = D, ncol =D)
  sum.mean.precision <- matrix(0, nrow = D, ncol =1)
  
  for ( z in 1:K) {
    sum.precision <- sum.precision + ro * S[z,1:D, 1:D]
    sum.mean.precision <-  sum.mean.precision + ro* S[z,1:D, 1:D] %*% as.matrix(mu[z,1:D]) 
  }
  precision.epsilon <-  sum.precision
  mean.epsilon <- solve(precision.epsilon) %*% (sum.mean.precision) 
  
  epsilon <- mvrnorm(n=1, mu = as.vector(mean.epsilon), Sigma = solve(precision.epsilon)) 
  
  
  
# Update W the Wishart parameter
  
  sum.w <- sum.precision <- matrix(0, nrow = D, ncol =D)
  for ( z in 1:K) {
    sum.w <- sum.w + beta * S[z,1:D, 1:D]
  }
  
  W <- rWishart(n = 1, df = beta * K , Sigma =  solve(sum.w ))
  
} 





f =1
for (o in 1:iter) {
  
  
  
  for (j in 1:length(activeclust)) {
    clust <- which(c==activeclust[j])
     tempmatrix <- as.matrix(Y[clust,1:D])
      tempmean <- apply(tempmatrix,2,mean)
      ## Some temporary matrices for calcuating the posterior
      ## The number of elements in the class
      temp.number <- as.numeric(numclust[activeclust[j]])
      tempmat1 <-  matrix(apply(as.matrix(apply(tempmatrix,1,function(x) (x -tempmean)%o% (x-tempmean))),1,sum), nrow=D,ncol=D)
      tempmat2 <- ((ro* temp.number)/(ro + temp.number)) * ((tempmean - epsilon ) %o% (tempmean - epsilon))
      
      Zi <- matrix(beta*W, nrow =2, ncol=2) + tempmat1 + tempmat2
      ## Posterior estimate for precision matrix
      S[activeclust[j],1:D,1:D] <- rWishart(1, (beta + temp.number), solve(Zi))
      ## Posterior estimate for the mean
      mu[activeclust[j],1:D] <- mvrnorm(n=1, mu = (ro*epsilon + (temp.number)* tempmean )/(ro + temp.number ), Sigma = solve((ro+ temp.number)*S[activeclust[j],1:D,1:D]))
    }
   
  
  
  
  ######################### THE HYPERPARAMETERS OF THE GMM #################################  
  
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
#   ro <- rgamma(1, shape = (nactive/2 + 0.25 ), scale = (as.numeric(sum.ro) +0.5)^-1)
#   
  # Update W the Wishart parameter
  
  sum.w <- sum.precision <- matrix(0, nrow = D, ncol =D)
  for ( z in 1:nactive) {
    sum.w <- sum.w + beta * S[activeclust[z],1:D, 1:D]
  }
  
  W <- rWishart(n = 1, df = beta * nactive + D, Sigma =  solve(D * InvCov + beta * sum.w ))
  
  
  if( o %% 10 == 0) {
  
  mu1est[f,1:D] <- mu[1,1:D]
  mu2est[f,1:D] <- mu[2,1:D]
  mu3est[f,1:D] <- mu[3,1:D]
  
  S1est[f,1:D,1:D] <- S[1,1:D,1:D]
  S2est[f,1:D,1:D] <- S[2,1:D,1:D]
  S3est[f,1:D,1:D] <- S[3,1:D,1:D]
  
  roest[f] <- ro
  
  epsilonest[f,1:D] <- epsilon
  
  West[f,1:D,1:D] <- W[1:D,1:D,1]
  f = f+1
  
  } 
  
 print(o/iter)
  
} 



## Analyzing the converged outputs


