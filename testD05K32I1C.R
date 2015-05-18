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
F = 3
D = 5
rel.D = D-2

A <- MixSim(BarOmega = 0.01 ,K = F, p = rel.D, int =c(0,100))


data.mu = array(data = NA, dim =c(F,D))
data.S = array(data = NA, dim =c(F,D,D))

for( i in 1:F){
  data.mu[i,1:rel.D] <- A$Mu[i,1:rel.D]
  data.S[i,1:rel.D,1:rel.D] <- A$S[1:rel.D,1:rel.D,i]
}

## Generating the points
y1 <- mvrnorm(n =40, mu = data.mu[1,1:rel.D], Sigma = data.S[1,1:rel.D,1:rel.D])
y2 <- mvrnorm(n =60, mu = data.mu[2,1:rel.D], Sigma = data.S[2,1:rel.D,1:rel.D])
y3 <- mvrnorm(n =50, mu = data.mu[3,1:rel.D], Sigma = data.S[3,1:rel.D,1:rel.D])
y <- rbind(y1,y2,y3)



y1 <- scale(y1)
y2 <- scale(y2)
y3 <- scale(y3)

beta1 <- c(-1,2,4)
beta2 <- c(10,2,1)
beta3 <- c(5,2,5)

time1.pur <- (t(beta1) %*% t(y1))
time2.pur <- (t(beta2) %*% t(y2))
time3.pur <- (t(beta3) %*% t(y3))


## To add Noise Variance and Mean we again use the package
B <- MixSim(BarOmega = 0.01 ,K = F, p = 1, int =c(0,100))

time.mean <- c(0,0,0)
time.sd <- c(0,0,0)

for( i in 1:F){
  time.mean[i] <- B$Mu[i,1]
  time.sd[i] <-   sqrt((B$S[1,1,i])^-1)
}


time1.noise <- rnorm(n = 40, mean = time.mean[1], sd = time.sd[1] )
time2.noise <- rnorm(n = 60, mean = time.mean[2], sd = time.sd[2] )
time3.noise <- rnorm(n = 50, mean = time.mean[3], sd =  time.sd[3])


time1 <- time1.pur + time1.noise
time2 <- time2.pur + time2.noise
time3 <- time3.pur + time3.noise

########################################################################################

time.real <-  as.vector(cbind(time1, time2, time3))

#######################################################################################






## Making the data with non relevant features

y.extra <- mvrnorm(n =150, mu = c(12,50), Sigma = matrix(c(5,-5,-5, 30), nrow = 2, ncol =2))


#######################################################################################################

Y <- cbind(y, y.extra)

#######################################################################################################

time <- c(time1, time2, time3)

#######################################################################################################


censoring <- rbinom(n = NROW(Y), size =1, prob =0.98)

right.censoring.time <- min(time)  


index.time <- which(censoring==0)
for ( q in 1:length(index.time)){
 time[index.time[q]] <- right.censoring.time
  
}


iter = 1000

################################# CALLING MAIN FUNCTION ###################################################

result <- NA

source('dirichletmixture.R')
result <- dirichletmixture(Y, time, censoring, iter,F) 
