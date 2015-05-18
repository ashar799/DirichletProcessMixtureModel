rm(list =ls())
epsilon.real <- c(6, 60)
W.real <- matrix(c(4,-5,-5, 81), nrow = 2, ncol =2)


S1 <- matrix(NA, nrow =2, ncol =2)
S2 <- matrix(NA, nrow =2, ncol =2)
S3 <- matrix(NA, nrow = 2, ncol =2)


ro.real = 1.6
beta.real = 4

S1[1:2,1:2] <-  rWishart(1, beta.real , solve(beta.real*W.real))
S2[1:2,1:2] <-  rWishart(1, beta.real , solve(beta.real*W.real))
S3[1:2,1:2] <-  rWishart(1, beta.real , solve(beta.real*W.real))

m1 <- mvrnorm(1, mu = epsilon.real, Sigma = solve(ro.real*S1))
m2 <- mvrnorm(1, mu = epsilon.real, Sigma = solve(ro.real*S2))
m3 <- mvrnorm(1, mu = epsilon.real, Sigma = solve(ro.real*S3))

y1 <- mvrnorm(n =120, mu = m1, Sigma = solve(S1))
y2 <- mvrnorm(n =150, mu = m2, Sigma = solve(S2))
y3 <- mvrnorm(n =130, mu = m3, Sigma = solve(S3))


Y <-  rbind(y1, y2, y3)
K= 3
D = NCOL(Y)
N = NROW(Y)
iter = 100

#### Initialization


alpha <- 1.8
c <- c(rep(1,120),rep(2,150),rep(3,130))

beta  = 4
ro = 1.6
epsilon = as.vector(apply(Y,2,mean))
W = cov(Y)

## Initialization of the parameters for Gaussian Mixture
mu = matrix(data = NA, nrow = K, ncol = D)
S = array(data = NA, dim =c(K,D,D))



disclass <- table(factor(c, levels = 1:3))
activeclass <- which(disclass!=0)

for ( j in 1:length(activeclass)){
  S[activeclass[j],1:D,1:D] <- rWishart(1, beta, solve((beta*W)))
  mu[activeclass[j],1:D] <- mvrnorm(n=1, mu = epsilon, Sigma = solve(ro*S[activeclass[j],1:D,1:D]))
}




cognate <- NA

numclust <- table(factor(c, levels = 1:K))
activeclust <- which(numclust!=0)
nactive <- length(activeclust)
InvCov <- solve(cov(Y))
meandata <- apply(Y, 2, mean )
meandata <-  as.matrix(meandata)

## Temporary matrices for storing the results

mu1est <- matrix(NA, nrow = 1000, ncol =D)
mu2est <- matrix(NA, nrow = 1000, ncol =D)
mu3est <- matrix(NA, nrow = 1000, ncol =D)

S1est <- array(NA, dim =c(1000,D,D))
S2est <- array(NA, dim =c(1000,D,D))
S3est <- array(NA, dim =c(1000,D,D))


epsilonest <- matrix(NA, nrow = 1000, ncol =D)

West <- array(NA, dim =c(1000,D,D))

roest <- c(NA)

## BURN-IN period

for (o in 1:50) {
  
  
  
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


