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

#################################### SIMULATED DATA####################################################


# Y =  faithful
## Actual Number of Components and dimension and dimension which are relevant
F = 3
D = 10
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
y3 <- mvrnorm(n =130, mu = data.mu[3,1:rel.D], Sigma = data.S[3,1:rel.D,1:rel.D])


y <- rbind(y1,y2,y3)



y1 <- scale(y1, center = TRUE, scale = TRUE)
y2 <- scale(y2, center= TRUE, scale = TRUE)
y3 <- scale(y3, center= TRUE, scale = TRUE)

beta1 <- c(-2,5,3,1,9,-5,1,6)
beta2 <- c(-9,2,3,-1,1,-1,-7,7)
beta3 <- c(-1,1,-1,-1,1,1,7,7)

time1.pur <- (t(beta1) %*% t(y1))
time2.pur <- (t(beta2) %*% t(y2))
time3.pur <- (t(beta3) %*% t(y3))



time1.mean <- 8
time2.mean <- 18
time3.mean <- 12

time1.noise <- rnorm(n = 120, mean = time1.mean, sd = 1 )
time2.noise <- rnorm(n = 150, mean = time2.mean, sd = 2 )
time3.noise <- rnorm(n = 130, mean = time3.mean, sd = 1.5 )



time1 <- time1.pur + time1.noise
time2 <- time2.pur + time2.noise
time3 <- time3.pur + time3.noise

## Making the data with non relevant features

y1.extra <- mvrnorm(n =120, mu = c(-2,20), Sigma = matrix(c(5,-5,-5, 30), nrow = 2, ncol =2))
y2.extra  <- mvrnorm(n =150, mu = c(0,40), Sigma = matrix(c(4, 0,0 , 100), nrow = 2, ncol =2))
y3.extra <-  mvrnorm(n =130, mu = c(10,10), Sigma = matrix(c(3, 12, 12, 80), nrow = 2, ncol =2))
y.extra <- rbind(y1.extra, y2.extra, y3.extra)

#######################################################################################################

Y <- cbind(y, y.extra)

#######################################################################################################

time <- c(time1, time2, time3)

#######################################################################################################


censoring <- rbinom(n = NROW(Y), size =1, prob =0.8)


iter = 300
################################# CALLING MAIN FUNCTION ###################################################
# Q = detectCores()
# cl<-makeCluster(Q-1)
# registerDoParallel(cl)

result <- NA

source('dirichletmixture.R')
result <- dirichletmixture(Y, time, censoring, iter,F) 

