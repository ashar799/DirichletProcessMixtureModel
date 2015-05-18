rm(list =ls())
setwd('/home/bit/ashar/Dropbox/Code/DPmixturemodel/DPplusAFT')
library(MASS)
library(mvtnorm)
library(ars)
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
#################################### SIMULATED DATA####################################################


# Y =  faithful
## Actual Number of Components and dimension and dimension which are relevant
F = 2
## Total number of dimensions
D = 5
## Total number of relevant dimensions
rel.D = D-2

A <- MixSim(BarOmega = 0.01 ,K = F, p = rel.D, int =c(0,100))


data.mu = array(data = NA, dim =c(F,D))
data.S = array(data = NA, dim =c(F,D,D))

for( i in 1:F){
  data.mu[i,1:rel.D] <- A$Mu[i,1:rel.D]
  data.S[i,1:rel.D,1:rel.D] <- A$S[1:rel.D,1:rel.D,i]
}

## Generating the points
y1 <- mvrnorm(n =120, mu = data.mu[1,1:rel.D], Sigma = data.S[1,1:rel.D,1:rel.D])
y2 <- mvrnorm(n =150, mu = data.mu[2,1:rel.D], Sigma = data.S[2,1:rel.D,1:rel.D])

y <- rbind(y1,y2)

y1 <- scale(y1)
y2 <- scale(y2)

## factor of the points
factor <- c(0)
index1 <- 1:nrow(y1)
t1 <- nrow(y1)+1
t2 <- nrow(y2) + nrow(y2)
index2 <- t1 : t2
factor[index1] <- 1
factor[index2] <- 2



beta1 <- c(-1,2,4)
beta2 <- c(10,2,1)


time1.pur <- (t(beta1) %*% t(y1))
time2.pur <- (t(beta2) %*% t(y2))



## To add Noise Variance and Beta0

time1.noise <- rnorm(n = 120, mean = 3, sd = 1 )
time2.noise <- rnorm(n = 150, mean = 8, sd = 2 )



time1 <- time1.pur + time1.noise
time2 <- time2.pur + time2.noise




## Making the data with non relevant features

y1.extra <- mvrnorm(n =120, mu = c(-2,20), Sigma = matrix(c(5,-5,-5, 30), nrow = 2, ncol =2))
y2.extra  <- mvrnorm(n =150, mu = c(0,40), Sigma = matrix(c(4, 0,0 , 100), nrow = 2, ncol =2))


y.extra <- rbind(y1.extra, y2.extra)

#######################################################################################################

Y <- cbind(y, y.extra)

#######################################################################################################

time <- c(time1, time2)

#######################################################################################################


censoring <- rbinom(n = NROW(Y), size =1, prob =0.8)


iter = 50

################################# CALLING MAIN FUNCTION ###################################################

result <- NA

source('dirichletmixture.R')
result <- dirichletmixture(Y, time, censoring, iter,F) 
