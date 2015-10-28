## This script selects the seed which g

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
library(mixAK)
library(mclust)
library(monomvn)

randindexi <- c(0)
for ( l in 1:30){
  
  set.seed(l)

#################################### SIMULATED DATA PROPERTIES ####################################################
## Number of points
N = 100

## Number of Clusters
F = 2

## Distribution of the points within three clusters

p.dist = c(0.7,0.3)

## Total Number of features D

D = 50

## Total Percentage of irrelevant feature

prob.noise.feature = 0.90

## Total Percentage of censoring

prob.censoring = 0.10

## Overlap between Cluster of molecular Data of the relevant features

prob.overlap = 0.05

## Percentage of Noise/Overlap in Time Data

prob.noise = 0.05

## Actual Number of Components and dimension  which are relevant
rel.D = as.integer(D* (1-prob.noise.feature))
## Actual Number of Irrelevant Componenets
irrel.D = D - rel.D


## Generating Data with overlap only with the relevant features
A <- MixSim(MaxOmega = prob.overlap ,K = F, p = rel.D, int =c(-1.5,1.5), lim = 1e08)


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
  mean <- runif(irrel.D,-1.5,1.5)
  Y.irrel.list[[i]] <- mvrnorm(n = as.integer(N * p.dist[i]), mu = mean, Sigma = diag(x =1, nrow = irrel.D, ncol = irrel.D))
}

### Combining the data with relevant and irrelevant columns
data.old <- list(0) 
for (i in 1:F){
  data.old[[i]] <-  cbind(Y.rel.list[[i]], Y.irrel.list[[i]]) 
}

############################################### MAKING Y from the clusters data #####################3
Y.old <- c(0)
for (i in 1:F){
  Y.old <- rbind(Y.old, data.old[[i]])
} 
Y.old <- Y.old[-1,]

#########################################################################################
X <- Y.old

rel.X <- as.matrix(X[,1:rel.D])

obj.qr <- qr(X)

rk <- obj.qr$rank

alpha <- qr.Q(obj.qr)[,1:rel.D]

gamma <- qr.Q(obj.qr)[,(1+rel.D):rk]

matT <- matrix(runif(n = rel.D*(rk -rel.D), min = -0.005, max= 0.005), nrow = rel.D, ncol = (rk -rel.D))

matP <- t(matT) %*% matT

max.eig <- eigen(matP)$values[1]

max.corr <- sqrt(max.eig)/sqrt(1 + max.eig)

linear.space <- gamma + alpha %*% matT

irrel.X <- matrix(NA, nrow = N, ncol = irrel.D)

for ( i in 1: irrel.D){
  
  matTemp <- matrix(runif(n = (rk -rel.D), min = -1.5, max= 1.5),  nrow = (rk-rel.D), ncol =1)
  irrel.X[,i] <- as.vector(linear.space %*% matTemp)
  
}

## Checking if the covariance is indeed small

cov.mat <- cov(rel.X,irrel.X)

boxplot(cov.mat)

## Building the full data matrix

X.full <- cbind(rel.X, irrel.X)


levelplot(cov(X.full[,1:20]))

Y <- X.full

#########################################################################################
##########################################################################################
##### Now WE DEAL WITH CLUSTERED DATA AND GENERATE NON RELEVANT FEATURES INDEPENENTLY ###


## Selcting the beta co-efficients

## The Co-efficients have to be obtained from uniform distribution between [-3,3]
beta.list <- list(0)
half <- rel.D/2
ohalf <- rel.D - half
for ( i in 1:F){
  beta.list[[i]] <- as.vector(rbind(runif(half, min = -3, max = -0.1), runif(ohalf, min = 0.1, max = 3)))
}

## The pure time is generated
time.pur.list <- list(0)
for ( i in 1:F){
  time.pur.list[[i]] <- t(beta.list[[i]]) %*% t(Y.rel.sc.list[[i]])
}

## Simulating Time Data which is ONE dimensional
time.cluster <- MixSim(MaxOmega = prob.noise, K = F, p = 1, int =c(0,3))

time.noise.list <- list(0)
for ( i in 1:F){
  time.noise.list[[i]] <- rnorm(as.integer(N * p.dist[i]), mean = time.cluster$Mu[i], sd = sqrt((time.cluster$S[i])))
}

time.list <- list(0)
for ( i in 1:F){
  time.list[[i]] <- time.pur.list[[i]] + time.noise.list[[i]]
}



## True Labels for the points
c.true <- c(0)
for ( i in 1:F){
  c.true <- rbind(as.matrix(c.true) , as.matrix(c(rep(i, as.integer(N * p.dist[i])))))  
}
c.true <- as.factor(c.true[-1,])


#######################################MAKING TIME from cluster data ########################################################
## Rela time without Censoring
time.real <- c(0)
for (i in 1:F){
  time.real <- cbind(time.real, time.list[[i]])
} 
time.real <- time.real[,-1]
time.real <- as.vector(unlist(time.real))


####################################### Adding CENSORING INFORMATION  ################################################################
## Adding the censoring information to the TIME

censoring <- rbinom(n = NROW(Y), size =1, prob = 1- prob.censoring)
right.censoring.time <- min(time.real)  

time <- time.real

index.time <- which(censoring==0)
for ( q in 1:length(index.time)){
  time[index.time[q]] <- right.censoring.time
  
}


## Boxplots for Vizualization of the time Data without censoring
boxplot(time.list)


### A Quick ManWhittney U / Kruskal test  test to check if the time's of the two cluster are significantly different
## wilcox.test(as.vector(time.list[[1]]), as.vector(time.list[[2]]), alternative = "two.sided")

kruskal.test(time, c.true)


## Let's calculate MSE between predicted times with the correct co-efficients and Real Time i.e. without censoring
Ytemp <- matrix(NA, nrow = N, ncol = D)
numclust <- table(factor(c.true, levels = 1:F))
activeclass<- which(numclust!=0)
for ( i in 1:length(activeclass)) {
  
  clust <- which(c.true == activeclass[i])
  
  Ytemp[clust,1:D] <- scale(Y[clust,1:D], center = TRUE, scale = TRUE)
}
time.predicted.correct <- c(0)
for ( i in 1:N){  
  time.predicted.correct[i] <-  rnorm(1, mean = time.cluster$Mu[c.true[i]] + t(beta.list[[c.true[i]]]) %*% Ytemp[i,1:rel.D], sd = sqrt(time.cluster$S[c.true[i]]))
}

source('calcrmse.R')
er_min <- calcrmse(time.real,time.predicted.correct)$rmse

## Let's See if the Clusters are separate WITH AND WITHOUT THE RELEVANT FEATURES
Y.rel <- c(0)
for (i in 1:F){
  Y.rel <- rbind(Y.rel, Y.rel.list[[i]])
} 
Y.rel <- Y.rel[-1,]

## Just the Relevant Features Data
pc <- prcomp(Y.rel)
pc.pred <- predict(pc,newdata = Y.rel)
plot(pc.pred[,1], pc.pred[,2], pch = 19,col = c.true)

## With all the features but CORRELATED DATA
pc <- prcomp(Y.old)
pc.pred <- predict(pc,newdata = Y.old)
plot(pc.pred[,1], pc.pred[,2], pch = 19,col = c.true)

## All features but uncorrelated data
pc <- prcomp(Y)
pc.pred <- predict(pc,newdata = Y)
plot(pc.pred[,1], pc.pred[,2], pch = 19,col = c.true)

############################# PARAMETERS for GIBB's SAMPLING ######################################

iter = 300
iter.burnin = 100
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

## Check What is the RMSE
er_random <- calcrmse(time.real,That)$rmse


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

## Check RMSE
er_kmeansblasso <- calcrmse(time.real,time.predicted)$rmse


##See How close is the predicted time with the real time
wilcox.test(as.vector(time.predicted)[index.time], as.vector(time.real)[index.time], paired = TRUE)


## Adjusted Initial Rand INDEX measure
randindexi[l] <- adjustedRandIndex(c.true,as.factor(c))

print(l)
}



## MAKING SOME BOXPLOT ##################
 load(file = '/home/bit/ashar/ExpressionSets/Simulations/Oneminuserror.RData')
load(file = '/home/bit/ashar/ExpressionSets/Simulations/TopTenRandIndex.RData')


pdf('/home/bit/ashar/Boxplots.pdf')
boxplot(toptenrand,oneminuserror, names = c("AdjustedRandIndex","1 - Aveg.NMSE"), main = "Statistics on 10 Simulation Runs", ylim =c(0.5,1))
leg.text <- c("samples = 200", "dims = 50", "cluster =2", "relevant.dims = 4")
legend("bottomright", leg.text)
dev.off()


