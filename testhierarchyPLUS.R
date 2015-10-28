rm(list =ls())

library(rms)
library(pec)
library(ipred)
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
library(survival)


D = 100
K = 50
ro = 1.6
beta = 105

W.real <- diag(runif(100,5,100))


## Initialization of the parameters for Gaussian Mixture
S = array(data = NA, dim =c(K,D,D))



for ( j in 1:K){
  S[j,1:D,1:D] <- rWishart(1, beta, solve((beta*W.real)))
 
}



## Initialize the desired variables
W = array(data = NA, dim =c(D,D))


## BURN-IN period

for (o in 1:50) {
  
  ######################### THE HYPERPARAMETERS OF THE GMM #################################  
  
  
#### THE NEW STRUCTURE OF W INVOLVES THAT IT BE A DAIGONAL MATRIX WITH alpha_a variables
  alpha_a <- c(rep(0,D))
  
  for ( i in 1:D){
    alpha_a[i] <- rgamma(n =1, shape = 0.5* beta*K, rate = 0.5*beta* sum(S[1:K,i,i]))    
  }
  W <- diag(alpha_a)
  
} 


