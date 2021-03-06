rm(list = ls())
setwd('/home/bit/ashar/Dropbox/Code/DPmixturemodel/DPplusAFT')
#setwd("C:/Users/Oana-Ashar/Desktop/Dropbox/Code/DPmixturemodel/DPplusAFT")

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


## Run my model on some toy data set
N = 200
D = 5
dminus = 4

## Generating the covaiance structure 
Sigm <- matrix(0, nrow = D, ncol =D)

L=  as.integer(D/2)
Q = D- L
for ( i in 1:L){
  for (j in 1:L){
    Sigm[i,j]=  3.5
  }
}
for ( i in Q:D){
  for (j in Q:D){
    Sigm[i,j]=  3.5
  }
}

diag(Sigm) <- 4

xto <- mvrnorm(n = N, mu = c(rep(4,D)), Sigma = Sigm)

## Generating the mixing distribution
x <- xto[,1]
y <- rep(0,N)
c.true <- rep(0,N)
  
for ( i in 1:N){
  pi1 <-  dnorm(x[i],mean = 4, sd = 0.5)
  pi2 <-  dnorm(x[i],mean = 6, sd = 0.5)
  c.true[i] <- sample(2, size =1, prob= c(pi1,pi2))
  if ( c.true[i] ==1){
    y[i] <- rnorm(1,mean = x[i], sd = 0.25)
  } else{
    y[i] <-  rnorm(1,mean = 4.5 + 0.1*x[i], sd = sqrt(1/8))
  }
  
}

time <- y
censoring <- rep(1,N)

Time <-  cbind(time,censoring)


Y <-xto

K = as.integer(N/2)

source('rchinese.R')
## Initialization of all the hyperparameters and 
shape.alpha <- 2
rate.alpha <- 1
alpha  = rgamma(1, shape = shape.alpha, rate = rate.alpha )
beta  = D+1
ro = 0.5

## Empirical Bayes Estimate of the Hyperparameters

epsilon = as.vector(apply(Y,2,mean))
W = cov(Y)
c <-  rchinese(N,alpha)
f <- table(factor(c, levels = 1:max(c)))

## Initialization of the parameters for Gaussian Mixture
mu = matrix(data = NA, nrow = K, ncol = D)
S = array(data = NA, dim =c(K,D,D))


#Sparsity controlling parameter
r =1
si = 1.78

lambda2 <- numeric(K)
tau2 = matrix(data = NA, nrow = K, ncol = D)
betahat = matrix(data = NA, nrow = K, ncol = D)
sigma2 <- rep(NA, K)
beta0 <- rep(NA, K)
That <-  numeric(N)

## Fitting a linear model to the whole model
Ysc <- scale(Y[1:N,1:D], center = TRUE, scale =TRUE)
lm.data <- lm(time ~ Ysc)
sig2.dat <-  var(lm.data$residuals)


## Set Some Initial Values for the Cluster Parameters

source('priordraw.R')
disclass <- table(factor(c, levels = 1:K))
activeclass <- which(disclass!=0)
for ( j in 1:length(activeclass)){
  
  priorone <- priordraw(beta, W, epsilon, ro, r, si, N, D, sig2.dat)  
  mu[activeclass[j],] <- priorone$mu
  S[activeclass[j],1:D,1:D]  <- priorone$Sigma  
  beta0[activeclass[j]] <- priorone$beta0 
  sigma2[activeclass[j]] <- priorone$sigma2
  betahat[activeclass[j],1:D] <- priorone$betahat 
  lambda2[activeclass[j]] <- priorone$lambda2 
  tau2[activeclass[j], 1:D] <- priorone$tau2
}

# The Time has to be initialized
source('updatetime.R')
ti <- updatetime(c, Y, Time,That, beta0, betahat, sigma2)
That <- ti$time

F =2

source('kmeansBlasso.R')
km <- kmeansBlasso(Y,That, F,K, beta, W, epsilon, ro, r, si, N, D, sig2.dat, c, mu, S, beta0, betahat, sigma2, lambda2, tau2)
c.init <- km$c
c <- c.init
mu <- km$mu
S <- km$S
sigma2 <- km$sigma2
betahat <- km$betahat
beta0 <- km$beta0
lambda2 <- km$lambda2
tau2 <- km$tau2

source('posteriorchineseAFT.R')
source('posteriorGMMparametrs.R')
source('posteriortimeparameters.R')
source('updatetime.R')
source('priordraw.R')
source('likelihood.R')
source('posterioralpha.R') 

cognate <- NA
param <- NA
paramtime <- NA
loglike<- rep(0, iter)  
timeparam <- NA
time.predicted <- c(0)
cindex <- c(0)

print(loglikelihood(c,Y,mu,S,alpha,That, beta0, betahat, sigma2, lambda2, tau2, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.dat) )

likli <- c(0)
o.initi = 1
iter.burnin = 200


print("BURNIN...PHASE")
for (o in o.initi:iter.burnin) {
  ################## PARAMETERS OF THE DP Mixture Model ######################################################
  ## Updating the parameters based on the observations 
  param <- posteriorGMMparametrs(c,Y,mu,S, alpha,K, epsilon, W, beta, ro,N,D )
  mu <- param$mean
  S <- param$precision
  paramtime <- posteriortimeparameters(c, That, lambda2,tau2,sigma2,beta0, betahat, Y, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.data)
  beta0 <- paramtime$beta0
  betahat <- paramtime$betahat
  sigma2 <- paramtime$sigma2
  lambda2 <- paramtime$lambda2
  tau2 <- paramtime$tau2
  
  ########################## THE HYPERPARAMETERS OF THE GMM #################################  
  source('posteriorhyper.R')  
  #  Updating the hyper paramters
  hypercognate <- posteriorhyper (c, Y, mu, S, epsilon, W, beta, ro )
  epsilon <- hypercognate$epsilon
  tmpW <- hypercognate$W
  W <- matrix(as.matrix(tmpW),nrow = D, ncol =D)
  ro <- hypercognate$ro
  
  ################# INDICATOR VARIABLE ##################################################################
  ## Updating the indicator variables and the parameters
  source('posteriorchineseAFT.R')
  cognate <- posteriorchineseAFT(c,Y,mu,S,alpha,That, beta0, betahat, sigma2, lambda2, tau2, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.dat)
  c <- cognate$indicator
  mu <- cognate$mean
  S <- cognate$precision
  beta0 <- cognate$beta0
  betahat <- cognate$betahat
  sigma2 <- cognate$sigma2
  lambda2 <- cognate$lambda2
  tau2 <- cognate$tau2
  
  ########################### The Concentration Parameter #################################################################
  
  
  source('posterioralpha.R') 
  # Updating the concentration parameter
  alpha <- posterioralpha(c, N, alpha, shape.alpha, rate.alpha)
  
  
  ##################### Print SOME Statistics #####################################################
  
  likli[o] <- loglikelihood(c,Y,mu,S,alpha,That, beta0, betahat, sigma2, lambda2, tau2, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.dat)
  print(likli[o])
  
  print(o/iter.burnin)
} 


## considering that we are now in the posterior region we can try to sample from the posterior and make predictions
## We make predictions on a uniform Ystar in this case between 0 and 1

N2 =  35
x1star <- seq(from = 0, to = 10, by = 10/(N2-1) )
# x.irrelstar <- mvrnorm(n = N2, mu = c(rep(4,dminus)), Sigma = diag(x =1, nrow = dminus, ncol = dminus))
# 
x.irrelstar = xto[1:N2,2:D]
xstar = cbind(x1star,x.irrelstar)
print("POSTERIOR...PREDICTION")


iter  = 100
thin =  2
o =1
Tstar <- matrix(0, nrow =N2, ncol = as.integer(iter/thin))

for (o in 1:iter) {
  ################## PARAMETERS OF THE DP Mixture Model ######################################################
  ## Updating the parameters based on the observations 
  param <- posteriorGMMparametrs(c,Y,mu,S, alpha,K, epsilon, W, beta, ro,N,D )
  mu <- param$mean
  S <- param$precision
  paramtime <- posteriortimeparameters(c, That, lambda2,tau2,sigma2,beta0, betahat, Y, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.data)
  beta0 <- paramtime$beta0
  betahat <- paramtime$betahat
  sigma2 <- paramtime$sigma2
  lambda2 <- paramtime$lambda2
  tau2 <- paramtime$tau2
  
  ########################## THE HYPERPARAMETERS OF THE GMM #################################  
  source('posteriorhyper.R')  
  #  Updating the hyper paramters
  hypercognate <- posteriorhyper (c, Y, mu, S, epsilon, W, beta, ro )
  epsilon <- hypercognate$epsilon
  tmpW <- hypercognate$W
  W <- matrix(as.matrix(tmpW),nrow = D, ncol =D)
  ro <- hypercognate$ro
  
  ################# INDICATOR VARIABLE ##################################################################
  ## Updating the indicator variables and the parameters
  source('posteriorchineseAFT.R')
  cognate <- posteriorchineseAFT(c,Y,mu,S,alpha,That, beta0, betahat, sigma2, lambda2, tau2, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.dat)
  c <- cognate$indicator
  mu <- cognate$mean
  S <- cognate$precision
  beta0 <- cognate$beta0
  betahat <- cognate$betahat
  sigma2 <- cognate$sigma2
  lambda2 <- cognate$lambda2
  tau2 <- cognate$tau2
  
  ########################### The Concentration Parameter #################################################################
  
  
  source('posterioralpha.R') 
  # Updating the concentration parameter
  alpha <- posterioralpha(c, N, alpha, shape.alpha, rate.alpha)
  

  ##################### GEt the predictions after every thin iteration #####################################################
  
  source('predict.R')
  if((o %% thin) == 0){
    to <-  o/thin
    resultat <- predict(c, Y, That, Time, beta0, betahat, sigma2, xstar) 
    Tstar[,to]  <- resultat$tstar
    
  }
  likli[o] <- loglikelihood(c,Y,mu,S,alpha,That, beta0, betahat, sigma2, lambda2, tau2, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.dat)
  
  
  print(o/iter)
} 

## Predictions
Tstar <- Tstar[,1:20]

postmat <- apply(Tstar,1,mean)
postsd <- apply(Tstar,1,sd)

pdf(file = "Wade.pdf")
p2D <- ggplot() + ggtitle("Predicted values") + geom_point(aes(x = xstar[,1], y = postmat)) + geom_errorbar(aes(x=xstar[,1], ymin=postmat-postsd, ymax=postmat+postsd), width=0.25) + geom_point(aes(x = Y[,1], y = time, color = factor(c)))+ scale_color_manual(values = c("red", "blue", "green")) + labs(x = "x", y= "y")
dev.off()
