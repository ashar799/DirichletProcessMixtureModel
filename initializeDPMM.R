initializeDPMM = function(){
  

################################# GIBBS SAMPLING  ###################################################

Time <- cbind(time, censoring) 

K <<-  as.integer(N/5)

## HYPER PRIORS
## Hyper parameters of the DP
shape.alpha <- 2
rate.alpha <- 1
## Hyperparameters for the GMM
beta  = D+1
ro = 0.5


source('rchinese.R')
alpha  = rgamma(1, shape = shape.alpha, rate = rate.alpha )
c <-  rchinese(N,alpha)
f <- table(factor(c, levels = 1:max(c)))

## Empirical Bayes Estimate of the Hyperparameters
epsilon = as.vector(apply(Y,2,mean))
W = diag(diag(cov(Y)))


## Initialization of the parameters for Gaussian Mixture
mu = matrix(data = NA, nrow = K, ncol = D)
S = array(data = NA, dim =c(K,D,D))


#Sparsity controlling hyperparameter of the BAYESIAN LASSO MODEL
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
  mu[activeclass[j],] <- (priorone$mu) 
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

## Check What is the RMSE
# source('calcrmse.R')
# er_random <<- calcrmse(time.real,That)$rmse


## Initialization part for the parmaters of AFT Model with k-means and Bayesian Lasso
source('kmeansBlasso.R')
source('flexmixBLASSO.R')
F = k
fm <- kmeansBlasso(Y,That, F,K, beta, W, epsilon, ro, r, si, N, D, sig2.dat, c, mu, S, beta0, betahat, sigma2, lambda2, tau2)
#F =2
# fm2 <- flexmixBlasso(Y,That, F,K, beta, W, epsilon, ro, r, si, N, D, sig2.dat, c, mu, S, beta0, betahat, sigma2, lambda2, tau2)

c <- fm$c
mu <- fm$mu
S <- fm$S
sigma2 <- fm$sigma2
betahat <- fm$betahat
beta0 <- fm$beta0
lambda2 <- fm$lambda2
tau2 <- fm$tau2



# Testing the  k-means estimate
source('predicttime.R')
time.predicted <- predicttime(c,Y, That,Time,beta0, betahat, sigma2)$predicttime

## Check RMSE
#source('calcrmse.R')
#er_kmeansblasso <<- calcrmse(time.real,time.predicted)$rmse


##See How close is the predicted time with the real time
#wilcox.test(as.vector(time.predicted), as.vector(time.real), paired = TRUE)


## Adjusted Initial Rand INDEX measure
randindexi <<- adjustedRandIndex(c.true,as.factor(c))



assign("Time", Time, envir = .GlobalEnv)
assign("r", r, envir = .GlobalEnv)
assign("si", si, envir = .GlobalEnv)
assign("shape.alpha", shape.alpha, envir = .GlobalEnv)
assign("rate.alpha", shape.alpha, envir = .GlobalEnv)
assign("alpha", alpha, envir = .GlobalEnv)
assign("beta", beta, envir = .GlobalEnv)
assign("ro", ro, envir = .GlobalEnv)
assign("That", That, envir = .GlobalEnv)
assign("c", c, envir = .GlobalEnv)
assign("epsilon", epsilon, envir = .GlobalEnv)
assign("W", W, envir = .GlobalEnv)
assign("mu", mu, envir = .GlobalEnv)
assign("S", S, envir = .GlobalEnv)
assign("beta0", beta0, envir = .GlobalEnv)
assign("betahat", betahat, envir = .GlobalEnv)
assign("sigma2", sigma2, envir = .GlobalEnv)
assign("lambda2", lambda2, envir = .GlobalEnv)
assign("tau2", tau2, envir = .GlobalEnv)
assign("lambda2", lambda2, envir = .GlobalEnv)
assign("sig2.dat", sig2.dat, envir = .GlobalEnv)

}
