#################################################################################################
#################################################################################################
################## THE MAIN FUNCTION ###########################################################
#################################################################################################

dirichletmixture = function(Y, time, censoring, iter, F ) {

setwd('/home/bit/ashar/Dropbox/Code/DPmixturemodel/DPplusAFT')
  
Time <- cbind(time, censoring) 
  
D = NCOL(Y)
N = NROW(Y)
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


## Setting the parameters from the prior
source('priorparameter.R')
prior <- priorparameter(c,S, mu, lambda2,tau2,sigma2,beta0, betahat,K, epsilon, W, beta, ro, D, r, si, Time, N)
mu <- prior$mu 
S <- prior$Sigma  
sigma2 <- prior$sigma2
betahat <- prior$betahat 
lambda2 <- prior$lambda2
tau2 <- prior$tau2
beta0 <- prior$beta0


## Initialization part for the parmaters of AFT Model with k-means
source('kmeansinit.R')
km <- kmeansinit(Y,time, N, F,D, K, mu, sigma2, betahat, beta0)
c <- km$c
mu <- km$mu
S <- km$S
sigma2 <- km$sigma2
betahat <- km$betahat
beta0 <- km$beta0
lambda2 <- km$lambda2
tau2 <- km$tau2

## Fitting a linear model to the whole model
Ysc <- scale(Y[1:N,1:D], center = TRUE, scale =TRUE)
lm.data <- lm(time ~ Ysc)
sig2.dat <-  var(lm.data$residuals)

# The Time has to be initialized
source('updatetime.R')
ti <- updatetime(c, Y, Time,That, beta0, betahat, sigma2)
That <- ti$time


# That <- Time[,1]

## MCMC sampling 
source('posteriorchineseAFT.R')
source('posteriorGMMparametrs.R')
source('posteriortimeparameters.R')
source('updatetime.R')
source('priordraw.R')
source('posterioralpha.R')
source('likelihood.R')
source('posttime.R')

cognate <- NA
param <- NA
paramtime <- NA
loglike<- rep(0, iter)  
timeparam <- NA

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
#   source('posteriorhyper.R')  
#   #  Updating the hyper paramters
#   hypercognate <- posteriorhyper (c, Y, mu, S, epsilon, W, beta, ro )
#   epsilon <- hypercognate$epsilon
#   W <- hypercognate$W
#   W <- matrix(as.matrix(W),nrow = D, ncol =D)
#   ro <- hypercognate$ro
#   
  
  ############################# PARMATERS OF THE TIME MODEL ######################################

 source('posttime.R') 
 timeparam <- posttime(c, That, lambda2,tau2,sigma2,beta0, betahat, Y, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.data ) 
 beta0 = timeparam$beta0
 sigma2 = timeparam$sigma2
 betahat <- timeparam$betahat 
 lambda2 <- timeparam$lambda2 
 tau2 <- timeparam$tau2

  
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


######################## The Censored Times ###########################################################
source('updatetime.R')
# Updating the Time Variable
ti <- NA
ti <- updatetime(c, Y, Time,That, beta0, betahat, sigma2)
That <- ti$time

 

## Value of the Log-likelihood
source('likelihood.R')
loglike[o] <-loglikelihood(c,Y,mu,S,alpha,That, beta0, betahat, sigma2, lambda2, tau2, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.dat) 

print(o/iter) 
print(loglike[o])

 } 
  
return(list('c' = c, 'time' = That, 'beta0' = beta0,'sigma2' = sigma2, 'betahat' = betahat, 'lambda2' = lambda2, 'tau2' =  tau2 ))

}
