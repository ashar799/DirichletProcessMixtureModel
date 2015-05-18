rm(list =ls())
setwd('/home/bit/ashar/Dropbox/Code/DPmixturemodel/DPplusAFT')
library(MASS)
library(mvtnorm)
library(matrixcalc)
library(stats)
library(Runuran)
library(truncnorm)
library(Matrix)
library(MCMCpack)
library(psych)
library(VGAM)


#################################### SIMULATED DATA####################################################
# Y =  faithful
y1 <- mvrnorm(n =120, mu = c(0,15), Sigma = matrix(c(0.9,-5,-5, 30), nrow = 2, ncol =2))
y2 <- mvrnorm(n =150, mu = c(3,80), Sigma = matrix(c(2,-2,-2, 100), nrow = 2, ncol =2))
y3 <- mvrnorm(n =130, mu = c(10,45), Sigma = matrix(c(3, 12, 12, 80), nrow = 2, ncol =2))


Y <-  rbind(y1, y2, y3)

y1 <- scale(y1, center = TRUE, scale = TRUE)
y2 <- scale(y2, center= TRUE, scale = TRUE)
y3 <- scale(y3, center= TRUE, scale = TRUE)

beta1 <- c(8,5)
beta2 <- c(1, 10)
beta3 <- c(3, 3)

time1.pur <- (t(beta1) %*% t(y1))
time2.pur <- (t(beta2) %*% t(y2))
time3.pur <- (t(beta3) %*% t(y3))


## Now I generate Noise which also has a betahat term

time1.mean <- 12
time2.mean <- 59
time3.mean <- 3



time1.noise <- rnorm(n = 120, mean = time1.mean, sd = 3 )
time2.noise <- rnorm(n = 150, mean = time2.mean, sd = 8 )
time3.noise <- rnorm(n = 130, mean = time3.mean, sd = 1 )


time1 <- time1.pur + time1.noise
time2 <- time2.pur + time2.noise
time3 <- time3.pur +  time3.noise

N =nrow(Y)

censoring <- rep(1,N)
time <- c(time1, time2, time3)


iter = 300
F =3


################################# CALLING MAIN FUNCTION ###################################################

result <- NA

source('dirichletmixture.R')
result <- dirichletmixture(Y, time, censoring, iter,F) 


