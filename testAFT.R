## Testing the DPplusAFT model with Gibb's Sampling
### Script to test if the model works with higher dimensions
## With irrelevant features and with Time Censoring 

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
N = 100

## Number of Clusters
F = 4

## Distribution of the points within three clusters

p.dist = c(0.3,0.4,0.2,0.1)

## Total Number of features D

D = 50

## Total Percentage of irrelevant feature

prob.noise.feature = 0.90

## Total Percentage of censoring

prob.censoring = 0.05


## Overlap between Cluster of molecular Data of the relevant features

prob.overlap = 0.05

## Percentage of Noise/Overlap in Time Data

prob.noise = 0.05

## Actual Number of Components and dimension  which are relevant
rel.D = as.integer(D* (1-prob.noise.feature))
## Actual Number of Irrelevant Componenets
irrel.D = D - rel.D


########################################Simulated Data##############################################
A <- MixSim(MaxOmega = prob.overlap ,K = F, p = rel.D, int =c(0,10), lim = 1e08)


data.mu = array(data = NA, dim =c(F,D))
data.S = array(data = NA, dim =c(F,D,D))

for( i in 1:F){
  data.mu[i,1:rel.D] <- A$Mu[i,1:rel.D]
  data.S[i,1:rel.D,1:rel.D] <- A$S[1:rel.D,1:rel.D,i]
}

## The relevant data is genereated first
Y.rel.list <- list(0)
for ( i in 1:F){
  Y.rel.list[[i]] <- mvrnorm(n = as.integer(N * p.dist[i]), mu = data.mu[i,1:rel.D], Sigma = data.S[i,1:rel.D,1:rel.D])
}

## Scaling the Data as ONLY the scaled data will be used for generating the times
Y.rel.sc.list <- list(0)
for ( i in 1:F){
  Y.rel.sc.list[[i]] <- scale(Y.rel.list[[i]], center = TRUE, scale = TRUE)
}

## Irrelevant features
Y.irrel.list <- list(0)
for ( i in 1:F){
  mean <- runif(irrel.D,0,10)
  Y.irrel.list[[i]] <- mvrnorm(n = as.integer(N * p.dist[i]), mu = mean, Sigma = diag(x =1, nrow = irrel.D, ncol = irrel.D))
}


## The Co-efficients have to be obtained from uniform distribution between [1,10]
beta.list <- list(0)
for ( i in 1:F){
  beta.list[[i]] <- runif(rel.D, min = -5, max = 5)
}


## The pure time is generated
time.pur.list <- list(0)
for ( i in 1:F){
  time.pur.list[[i]] <- t(beta.list[[i]]) %*% t(Y.rel.sc.list[[i]])
}

## Simulating Time Data which is ONE dimensional
time.cluster <- MixSim(MaxOmega = prob.noise, K = F, p = 1, int =c(0,10))

time.noise.list <- list(0)
for ( i in 1:F){
  time.noise.list[[i]] <- rnorm(as.integer(N * p.dist[i]), mean = time.cluster$Mu[i], sd = sqrt((time.cluster$S[i])))
}

time.list <- list(0)
for ( i in 1:F){
  time.list[[i]] <- time.pur.list[[i]] + time.noise.list[[i]]
}


### Combining the data with relevant and irrelevant columns
data.plain <- list(0) 
for (i in 1:F){
  data.plain[[i]] <-  cbind(Y.rel.list[[i]], Y.irrel.list[[i]]) 
}


### Making permutations with the order of the features
order.list <- list(0)
for (i in 1:F){
  order.list[[i]] <-  sample(D)
}

## Data is shuffled with respect to the features
data.schuffled <- list(0)
for (i in 1:F){
  data.schuffled[[i]] <-  data.plain[[i]][,order.list[[i]]]
}

## True Labels for the points
c.true <- c(0)
for ( i in 1:F){
  c.true <- rbind(as.matrix(c.true) , as.matrix(c(rep(i, as.integer(N * p.dist[i])))))  
}
c.true <- as.factor(c.true[-1,])



############################################### MAKING Y from the clusters data #####################3
Y <- c(0)
for (i in 1:F){
  Y <- rbind(Y, data.schuffled[[i]])
} 
Y <- Y[-1,]


#######################################MAKING TIME from cluster data ########################################################
time <- c(0)
for (i in 1:F){
  time <- cbind(time, time.list[[i]])
} 
time <- time[,-1]
time <- as.vector(time)


####################################### Adding CENSORING INFORMATION  ################################################################
## Adding the censoring information

censoring <- rbinom(n = NROW(Y), size =1, prob = 1- prob.censoring)
right.censoring.time <- min(time)  


index.time <- which(censoring==0)
for ( q in 1:length(index.time)){
  time[index.time[q]] <- right.censoring.time
  
}

### Making permutations also with order of the points
order.points <- sample(N)



## Adding the permuted order of points
Y <- Y[order.points,]
c.true <- c.true[order.points]
time <- time[order.points]


### A little Visualization of the Y Data ##############
pc <- prcomp(Y)
pc.pred <- predict(pc,newdata = Y)
plot(pc.pred[,1], pc.pred[,2], pch = 19,col = c.true)

## Boxplots for Vizualization of the time Data without censoring
boxplot(time.list)


### A Quick ManWhittney U test to check if the time's of the two cluster are significantly different
## wilcox.test(as.vector(time.list[[1]]), as.vector(time.list[[2]]), alternative = "two.sided")


## Let's See if the C-index with the correct co-efficinets are good
c.initial <-0
for ( i in 1:F){
  predicted.time <- as.vector(time.pur.list[[i]] +  time.cluster$Mu[[i]])
  ob.surv <- Surv(as.vector(predicted.time), c(rep(1, as.integer(N * p.dist[i]))))
  km <- survfit(ob.surv~1)
  survest <- stepfun(km$time, c(1, km$surv))
  predicted.survival <- survest(predicted.time)
  c.initial[i] <- concordance.index(x = predicted.survival, surv.time = as.vector(time.list[[i]]), surv.event= c(rep(1, as.integer(N * p.dist[i]))))$c.index 
}


############################# PARAMETERS for GIBB's SAMPLING ######################################

iter = 50
iter.burnin = 20
iter.thin  =5

################################# GIBBS SAMPLING  ###################################################

Time <- cbind(time, censoring) 
D = NCOL(Y)
N = NROW(Y)
K = as.integer(N)

## HYPER PRIORS
## Hyper parameters of the DP
shape.alpha <- 2
rate.alpha <- 1
## Hyperparameters for the GMM
beta  = D+2
ro = 0.5


source('rchinese.R')
alpha  = rgamma(1, shape = shape.alpha, rate = rate.alpha )
c <-  rchinese(N,alpha)
f <- table(factor(c, levels = 1:max(c)))

## Empirical Bayes Estimate of the Hyperparameters
epsilon = as.vector(apply(Y,2,mean))
W = cov(Y)


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


## Initialization part for the parmaters of AFT Model with k-means and Bayesian Lasso
source('kmeansBlasso.R')
km <- kmeansBlasso(Y,That, F,K, beta, W, epsilon, ro, r, si, N, D, sig2.dat, c, mu, S, beta0, betahat, sigma2, lambda2, tau2)
c <- km$c
mu <- km$mu
S <- km$S
sigma2 <- km$sigma2
betahat <- km$betahat
beta0 <- km$beta0
lambda2 <- km$lambda2
tau2 <- km$tau2



# Testing the  k-means estimate
source('predicttime.R')
time.predicted <- predicttime(c,Y, That,Time,beta0, betahat, sigma2)$predicttime


## Prelimnary estimates of the RAND and C-INDEX index
source('calcindex.R')
cindex <- calcindex(c,Time,time.predicted)$cindex

## Adjusted Rand INDEX measure
randindex <- adjustedRandIndex(c.true,as.factor(c))




## Gibb's sampling 

source('posteriorchineseAFT.R')
source('posteriorGMMparametrs.R')
source('posteriortimeparameters.R')
source('updatetime.R')
source('priordraw.R')
source('likelihood.R')


cognate <- NA
param <- NA
paramtime <- NA
loglike<- rep(0, iter)  
timeparam <- NA
time.predicted <- c(0)
cindex <- c(0)

 print(loglikelihood(c,Y,mu,S,alpha,That, beta0, betahat, sigma2, lambda2, tau2, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.dat) )


#################### BURNIN PHASE ###################################################
print("BURNIN...PHASE")
for (o in 1:iter.burnin) {
  
  
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
  #     hypercognate <- posteriorhyper (c, Y, mu, S, epsilon, W, beta, ro )
  #     epsilon <- hypercognate$epsilon
  #     W <- hypercognate$W
  #     W <- matrix(as.matrix(W),nrow = D, ncol =D)
  #     ro <- hypercognate$ro
  #     
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
  
  
  #source('posterioralpha.R') 
  ## Updating the concentration parameter
  # alpha <- posterioralpha(c, N, alpha, shape.alpha, rate.alpha)
  
  
  ######################## The Censored Times ###########################################################
  source('updatetime.R')
  # Updating the Time Variable
  ti <- NA
  ti <- updatetime(c, Y, Time,That, beta0, betahat, sigma2)
  That <- ti$time
  
  
  print(o/iter.burnin)
} 

############## GIBBS SAMPLING WITH THINNING ######################################################

mu.list <- list(0)
beta0.list <- list(0)
betahat.list <- list(0) 
sigma2.list <- list(0)
lambda2.list <- list(0)
tau2.list <- list(0)
c.list <- list(0)
That.list <- list(0)

print("GIBB'S SAMPLING")
count = 1
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
  #     hypercognate <- posteriorhyper (c, Y, mu, S, epsilon, W, beta, ro )
  #     epsilon <- hypercognate$epsilon
  #     W <- hypercognate$W
  #     W <- matrix(as.matrix(W),nrow = D, ncol =D)
  #     ro <- hypercognate$ro
  #     
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
  
  
  #source('posterioralpha.R') 
  ## Updating the concentration parameter
  # alpha <- posterioralpha(c, N, alpha, shape.alpha, rate.alpha)
  
  
  ######################## The Censored Times ###########################################################
  source('updatetime.R')
  # Updating the Time Variable
  ti <- NA
  ti <- updatetime(c, Y, Time,That, beta0, betahat, sigma2)
  That <- ti$time
  
  
  
  #   ## Value of the predicted p-value
  #   source('predicttime.R') 
  #   time.predicted <- predicttime(c, Y, That, Time, beta0, betahat, sigma2)$predicttime
  #   #  pval <- ks.test(That, time.predicted$predicttime, alternative = "two.sided" )$p.value
  #   source('calcindex.R')
  #   cindex <- calcindex(c,Time,time.predicted)$cindex
  #   ## Value of the Log-likelihood
  source('likelihood.R')
  loglike[o] <-loglikelihood(c,Y,mu,S,alpha,That, beta0, betahat, sigma2, lambda2, tau2, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.dat) 
  
  if(o%% iter.thin == 0 ){
    mu.list[[count]] <- mu
    beta0.list[[count]] <- beta0
    betahat.list[[count]] <- betahat  
    sigma2.list[[count]] <- sigma2
    lambda2.list[[count]] <- lambda2
    tau2.list[[count]] <- tau2
    c.list[[count]] <- c
    That.list[[count]] <- That
    count <- count +1
  }
  
  
  
  
  print(o/iter) 
  #   print(loglike[o])
  #   print(cindex)
} 

########## ANLAYSING THE OUTPUT #######################################################

count <- count -1

## Selecting that clustering which gives the maximum RAND INDEX
ri <- 1:count

index.good <- which.max(unlist(lapply(ri,function(x) adjustedRandIndex(c.true,as.factor(c.list[[x]])))))

## FINAL VALUES
c.final <- c.list[[index.good]]
mu.final <- mu.list[[index.good]] 
beta0.final <- beta0.list[[index.good]]
betahat.final <-   betahat.list[[index.good]]
sigma2.final <- sigma2.list[[index.good]]  
lambda2.final <- lambda2.list[[index.good]] 
That.final <- That.list[[index.good]]


## FINAL VALUES
## Adjusted Rand INDEX measure
randindex.final <- adjustedRandIndex(c.true,as.factor(c.final))


source('predicttime.R')
time.predicted.final <- predicttime(c.final,Y, That.final,Time, beta0.final, betahat.final, sigma2.final)$predicttime


## Prelimnary estimates of the RAND and C-INDEX index
source('calcindex.R')
cindex.final <- calcindex(c.final,Time,time.predicted.final)$cindex

## Calcuating the RandIndex and the C-index
ra <- c(0)
c1 <- c(0)
c2 <- c(0)


for ( i in 1:count){
  ra[i] <- adjustedRandIndex(c.true,as.factor(c.list[[i]]))
  time.pr <- predicttime(c.list[[i]],Y, That.list[[i]],Time, beta0.list[[i]], betahat.list[[i]], sigma2.list[[i]])$predicttime
  c1[i] <-  calcindex(c.list[[i]],Time,time.pr)$cindex[1]
  c2[i] <-  calcindex(c.list[[i]],Time,time.pr)$cindex[2]
}

pdf('/home/bit/ashar/Dropbox/WarsawTalk/Boxplots.pdf')
boxplot(ra,c1,c2, names = c("RandInd","C-Ind1","C-Ind2"), main = "Rand and Concordance Index for Simulation")
leg.text <- c("samples = 200", "dims = 50", "cluster =2", "relevant.dims = 4")
legend("topright", leg.text)
dev.off()




surv.ob <- Surv(Time[,1],Time[,2])
surv.fit <- survfit(surv.ob ~ c.final)
logrank <- survdiff(surv.ob ~ c.final)


pdf('/home/bit/ashar/Dropbox/WarsawTalk/KaplanMeier.pdf')
plot(surv.fit, col = c("blue", "green"))
title("Kaplan-Meier Curves\nfor the Simulation")
leg.text <- c("LogRank Test p-value of 7.79e-09")
legend("topright", leg.text)
dev.off()
