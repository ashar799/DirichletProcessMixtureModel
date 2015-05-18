####### To check if the GMM is working or not #####################33


rm(list = ls())
setwd('/home/bit/ashar/Dropbox/Code/DPmixturemodel/DPplusAFT')

library(MASS)
library(mixtools)
library(matrixcalc)
library(stats)
library(Runuran)
library(truncnorm)
library(Matrix)
library(MCMCpack)
library(psych)
library(VGAM)
library(MixSim)
library(statmod)
library(flexclust)
library(survcomp)
library(mixAK)
library(mclust)
#################################### SIMULATED DATA PROPERTIES ####################################################
## Number of points
N = 120

## Number of Clusters
F =3

## Distribution of the points within three clusters

p.dist = c(0.3,0.4,0.3)

## Total Number of features

D = 50

rel.D =5

irrel.D =45


########################################Simulated Data##############################################
A <- MixSim(MaxOmega = 0.10 ,K = F, p = rel.D, int =c(0,100))


data.mu = array(data = NA, dim =c(F,D))
data.S = array(data = NA, dim =c(F,D,D))

for( i in 1:F){
  data.mu[i,1:rel.D] <- A$Mu[i,1:rel.D]
  data.S[i,1:rel.D,1:rel.D] <- A$S[1:rel.D,1:rel.D,i]
}

## Initializing the points
Y.list <- list(0)
for ( i in 1:F){
  Y.list[[i]] <- matrix(0, nrow =  as.integer(N * p.dist[i]), ncol = D)
}

## The relevant data is genereated first
Y.rel.list <- list(0)
for ( i in 1:F){
  Y.rel.list[[i]] <- mvrnorm(n = as.integer(N * p.dist[i]), mu = data.mu[i,1:rel.D], Sigma = data.S[i,1:rel.D,1:rel.D])
  
}

Y.irrel.list <- list(0)
for ( i in 1:F){
  mean <- runif(irrel.D,0,100)
  Y.irrel.list[[i]] <- mvrnorm(n = as.integer(N * p.dist[i]), mu = mean, Sigma = diag(x =1, nrow = irrel.D, ncol = irrel.D))
}

data.plain <- list(0) 
for (i in 1:F){
  data.plain[[i]] <-  cbind(Y.rel.list[[i]], Y.irrel.list[[i]]) 
}





############################################### MAKING Y from the clusters data #####################3
Y <- c(0)
for (i in 1:F){
  Y <- rbind(Y, data.plain[[i]])
} 
Y <- Y[-1,]

####################
c.true <- as.factor(c(rep(1, as.integer(N * p.dist[1])),rep(2,as.integer(N * p.dist[2])), rep(3,as.integer(N * p.dist[3]))))

c <- c.true
##########################
## The true labels

for (i in 1:3){
  data.mu[i ,6:50] <- apply(Y.irrel.list[[i]],2, mean)
  
 }

K = F

shape.alpha <- 2
rate.alpha <- 1
## Hyperparameters for the GMM
beta  = D+2
ro = 0.5

## Empirical Bayes Estimate of the Hyperparameters
epsilon = as.vector(apply(Y,2,mean))
W = cov(Y)


## Initialization of the parameters for Gaussian Mixture
mu = matrix(data = NA, nrow = K, ncol = D)
S = array(data = NA, dim =c(K,D,D))

## Randomly initialize mu and S from the Prior Distribution

source('priordraw.R')
disclass <- table(factor(c, levels = 1:K))
activeclass <- which(disclass!=0)
for ( j in 1:length(activeclass)){
  
  S[activeclass[j],1:D,1:D] <- rWISHART(1, beta, solve((beta*W)))
  mu[activeclass[j],1:D] <- as.matrix(rMVN(n=1, mean = epsilon, Q = ro*S[activeclass[j],1:D,1:D])$x)
  
}

## Now Let's Calculate the Gibb's Update
source('posteriorGMMparametrs.R')
 mean <- list(0)

for (i in 1:100){
  param <- posteriorGMMparametrs(c,Y,mu,S, alpha,K, epsilon, W, beta, ro,N,D )
  mean[[i]] <- param$mean
  }

mu.test <- matrix(data = 0, nrow = K, ncol = D)
for ( i in 1:100 ){
  mu.test <- mu.test +mean[[i]]
} 
mu.test <-  0.01 * mu.test
########################################### GMM is WORKING #############################################################
##################################################################################################
#################################################################################################
#################################################################################################

## TIME PARAMETER UPDATE IS WORKING ##################################################

rm(list = ls())
setwd('/home/bit/ashar/Dropbox/Code/DPmixturemodel/DPplusAFT')

library(MASS)
library(mixtools)
library(matrixcalc)
library(stats)
library(Runuran)
library(truncnorm)
library(Matrix)
library(MCMCpack)
library(psych)
library(VGAM)
library(MixSim)
library(statmod)
library(flexclust)
library(survcomp)
library(mixAK)
library(mclust)
library(monomvn)


#################################### SIMULATED DATA PROPERTIES ####################################################
## Number of points
N = 120

## Number of Clusters
F =3

## Distribution of the points within three clusters

p.dist = c(0.3,0.4,0.3)

## Total Number of features

D =30

## Total Percentage of irrelevant feature

prob.noise.feature = 0.80

## Total Percentage of censoring

prob.censoring = 0.05


## Overlap between Cluster

prob.overlap = 0.10

## Percentage of Noise/Overlap in Time Data

prob.noise = 0.10

## Actual Number of Components and dimension and dimension which are relevant
rel.D = as.integer(D* (1-prob.noise.feature))
## Actual Number of Irrelevant Componenets
irrel.D = D - rel.D


########################################Simulated Data##############################################
A <- MixSim(MaxOmega = prob.overlap ,K = F, p = rel.D, int =c(0,100))


data.mu = array(data = NA, dim =c(F,D))
data.S = array(data = NA, dim =c(F,D,D))

for( i in 1:F){
  data.mu[i,1:rel.D] <- A$Mu[i,1:rel.D]
  data.S[i,1:rel.D,1:rel.D] <- A$S[1:rel.D,1:rel.D,i]
}

## Initializing the points
Y.list <- list(0)
for ( i in 1:F){
  Y.list[[i]] <- matrix(0, nrow =  as.integer(N * p.dist[i]), ncol = D)
}

## The relevant data is genereated first
Y.rel.list <- list(0)
for ( i in 1:F){
  Y.rel.list[[i]] <- mvrnorm(n = as.integer(N * p.dist[i]), mu = data.mu[i,1:rel.D], Sigma = diag(x =1, nrow = rel.D, ncol = rel.D) )
  
}



Y.rel.scale.list <- list(0)
for ( i in 1:F){
  Y.rel.scale.list[[i]] <- scale(Y.rel.list[[i]], center = TRUE, scale = FALSE)
}



Y.irrel.list <- list(0)
for ( i in 1:F){
  mean <- runif(irrel.D,0,100)
  Y.irrel.list[[i]] <- mvrnorm(n = as.integer(N * p.dist[i]), mu = mean, Sigma = diag(x =1, nrow = irrel.D, ncol = irrel.D))
}

## The true labels
for (i in 1:3){
  data.mu[i ,6:30] <- apply(Y.irrel.list[[i]],2, mean)
}

for (i in 1:3){
data.S[i,1:D,1:D] <- diag(x =1, nrow = D, ncol = D)
}

### Combining the data with relevant and irrelevant columns

data.plain <- list(0) 
for (i in 1:F){
  data.plain[[i]] <-  cbind(Y.rel.list[[i]], Y.irrel.list[[i]]) 
}


############################################### MAKING Y from the clusters data #####################3
Y <- c(0)
for (i in 1:F){
  Y <- rbind(Y, data.plain[[i]])
} 
Y <- Y[-1,]

#######################################

## Empirical Bayes Estimate of the Hyperparameters
epsilon = as.vector(apply(Y,2,mean))
W = cov(Y)
K =F

## Initialization of the parameters for Gaussian Mixture
mu = matrix(data = NA, nrow = K, ncol = D)
S = array(data = NA, dim =c(K,D,D))

for (i in 1:3){
  mu[i,1:D] <- data.mu[i ,1:D]
}

for (i in 1:3){
  S[i,1:D,1:D] <- data.S[i,1:D,1:D]  
}

##############################################################
###### WITH THE PARMATERS OF THE GMM KNOWN #####################
####### WE SUMULATE THE TIME DATA ################################


## The Co-efficients have to be obtained from uniform distribution between [1,10]
beta.list <- list(0)
for ( i in 1:F){
  beta.list[[i]] <- runif(rel.D, min =1, max = 10)
}


## The pure time is generated
time.pur.list <- list(0)
for ( i in 1:F){
  time.pur.list[[i]] <- t(beta.list[[i]]) %*% t(Y.rel.scale.list[[i]])
}

## Simulating Time Data which is ONE dimensional
time.cluster <- MixSim(MaxOmega = prob.noise, K = F, p = 1, int =c(0,100))

time.noise.list <- list(0)
for ( i in 1:F){
  time.noise.list[[i]] <- rnorm(as.integer(N * p.dist[i]), mean = time.cluster$Mu[i], sd = 1)
}

time.list <- list(0)
for ( i in 1:F){
  time.list[[i]] <- time.pur.list[[i]] + time.noise.list[[i]]
}

#######################################MAKING TIME from cluster data ########################################################
time <- c(0)
for (i in 1:F){
  time <- cbind(time, time.list[[i]])
} 
time <- time[,-1]
time <- as.vector(time)

c.true <- as.factor(c(rep(1, as.integer(N * p.dist[1])),rep(2,as.integer(N * p.dist[2])), rep(3,as.integer(N * p.dist[3]))))

#################################################
##### Initializing the Time Parameters ##########
#################################################
epsilon = as.vector(apply(Y,2,mean))
W = cov(Y)
beta = D+2
ro = 0.5
r =1
si = 1.78

lambda2 <- numeric(K)
tau2 = matrix(data = NA, nrow = K, ncol = D)
betahat = matrix(data = NA, nrow = K, ncol = D)
sigma2 <- rep(NA, K)
beta0 <- rep(NA, K)
That <-  numeric(N)

Ysc <- scale(Y[1:N,1:D], center = TRUE, scale =TRUE)
lm.data <- lm(time ~ Ysc)
sig2.dat <-  var(lm.data$residuals)


c <- c.true

source('priordraw.R')
disclass <- table(factor(c, levels = 1:K))
activeclass <- which(disclass!=0)
for ( j in 1:length(activeclass)){
  
  priorone <- priordraw(beta, W, epsilon, ro, r, si, N, D, sig2.dat)  
  beta0[activeclass[j]] <- priorone$beta0 
  sigma2[activeclass[j]] <- priorone$sigma2
  betahat[activeclass[j],1:D] <- priorone$betahat 
  lambda2[activeclass[j]] <- priorone$lambda2 
  tau2[activeclass[j], 1:D] <- priorone$tau2
}

############ Running 100 Gibb's Sampling to check ###########################33

source('posteriortimeparameters.R')

lambda2.list <- list(0)
tau2.list <- list(0)
betahat.list <- list(0)
sigma2.list <- list(0)
beta0.list <- list(0)

That <- time



for (i in 1:1000){
  paramtime <- posteriortimeparameters(c, That, lambda2,tau2,sigma2,beta0, betahat, Y, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.data)
  beta0 <- paramtime$beta0
  betahat <- paramtime$betahat
  sigma2 <- paramtime$sigma2
  lambda2 <- paramtime$lambda2
  tau2 <- paramtime$tau2
  
  beta0.list[[i]] <-  beta0
  betahat.list[[i]] <-  betahat
  sigma2.list[[i]] <- sigma2
  lambda2.list[[i]] <- lambda2
  tau2.list[[i]] <- tau2
  
  print(i/1000) 
  
}

beta0.test <- matrix(data = 0, nrow = K, ncol = 1)
betahat.test <- matrix(data = 0, nrow = K, ncol =D)
sigma2.test <- matrix(data = 0, nrow = K, ncol =1)
lambda2.test <- matrix(data = 0, nrow = K, ncol =1)
tau2.test <-  matrix(data = 0, nrow = K, ncol =D)



for ( i in 1:1000 ){
  beta0.test <- beta0.test + beta0.list[[i]]
  betahat.test <- betahat.test + betahat.list[[i]]
  sigma2.test <-  sigma2.test + sigma2.list[[i]]
  lambda2.test <- lambda2.test + lambda2.list[[i]]
  tau2.test <- tau2.test + tau2.list[[i]]
} 

beta0.test <-  0.001 * beta0.test
betahat.test <- 0.001 * betahat.test
sigma2.test <- 0.001*sigma2.test
lambda2.test <- 0.001* lambda2.test
tau2.test <- 0.001 * tau2.test

############################################################################
### Changed the Posterior Update for the Time Model Parameters ############


### Assume that TIME DATA IS GENERATED ONLY USING CENTRED REGRESSORS AND NOT SCALED REGRESSORS###



#########################################################
#### Using MONOMVN package ##############################


reg.blas <- blasso(data.plain[[1]], time.list[[1]], normalize = FALSE)
s <- summary(reg.blas, burnin=200)

