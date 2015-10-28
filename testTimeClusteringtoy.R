### This script tests whether the model can recover clustering structure in the absence of clustering in the the molecular data (or very high overlap in the time clusters)
### This script uses Clustering from time to initialize the cluster assignments


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
library(cluster)
library(flexmix)

randindexi <- c(0)

##### THE PART WITH THE CORRECT SEED ###########################
#################################### SIMULATED DATA PROPERTIES ####################################################
## Number of points
N = 150

## Number of Clusters
F = 2

## Distribution of the points within three clusters

p.dist = c(0.7,0.3)

## Total Number of features D

D = 4

## Total Percentage of irrelevant feature

prob.noise.feature = 0.50

## Total Percentage of censoring

prob.censoring = 0.10

## Overlap between Cluster of molecular Data of the relevant features

prob.overlap = 0.05

## Percentage of Noise/Overlap in Time Data

prob.noise = 0.01

## Actual Number of Components and dimension  which are relevant
rel.D = as.integer(D* (1-prob.noise.feature))
## Actual Number of Irrelevant Componenets
irrel.D = D - rel.D

## The data
A <- MixSim(BarOmega = prob.overlap ,K = F, p = rel.D, int =c(-1.0,1.0), lim = 1e09)


data.mu = array(data = NA, dim =c(F,D))
data.S = array(data = NA, dim =c(F,D,D))

for( i in 1:F){
  data.mu[i,1:rel.D] <- A$Mu[i,1:rel.D]
  data.S[i,1:rel.D,1:rel.D] <- A$S[1:rel.D,1:rel.D,i]
}

## The relevant data is genereated first ASSUMING DATA COMING ONLY FROM FIRST CLUSTER
Y.rel.list <- list(0)
for ( i in 1:F){
  Y.rel.list[[i]] <- mvrnorm(n = as.integer(N * p.dist[i]), mu = data.mu[1,1:rel.D], Sigma = data.S[1,1:rel.D,1:rel.D])
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

######################################################################################### Making the irrelevant features independent from the dependent features #############
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


levelplot(cov(X.full[,1:D]))

Y <- X.full

#########################################################################################
##########################################################################################
##### Now WE DEAL WITH CLUSTERED DATA AND GENERATE NON RELEVANT FEATURES INDEPENENTLY ###


## Selcting the beta co-efficients

## The Co-efficients have to be obtained from uniform distribution between [-3,3]
beta.list <- list(0)
for ( i in 1:F){
  beta.list[[i]] <- as.vector(rbind(2*i,i-6))
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
  time.noise.list[[i]] <- rnorm(as.integer(N * p.dist[i]), mean = 3*i, sd = 0.2 *i)
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
## Real time without Censoring
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
beta  = D+1
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
source('flexmixBLASSO.R')
fm <- flexmixBlasso(Y,That, F,K, beta, W, epsilon, ro, r, si, N, D, sig2.dat, c, mu, S, beta0, betahat, sigma2, lambda2, tau2)
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
er_kmeansblasso <- calcrmse(time.real,time.predicted)$rmse


##See How close is the predicted time with the real time
wilcox.test(as.vector(time.predicted)[index.time], as.vector(time.real)[index.time], paired = TRUE)


## Adjusted Initial Rand INDEX measure
 randindexi <- adjustedRandIndex(c.true,as.factor(c))



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

init.likli <- loglikelihood(c,Y,mu,S,alpha,That, beta0, betahat, sigma2, lambda2, tau2, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.dat) 
rmse <- c(0)
randy <- c(0)
likli <- c(0)

o =1
#################### BURNIN PHASE ###################################################
o.iter = o
print("BURNIN...PHASE")
for (o in o.iter:iter.burnin) {
  
  
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
  
  source('posteriorbeta.R')
  beta <- posteriorbeta(c, beta, D, S, W)
  
  
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
  
  
  ##################### Print SOME Statistics #####################################################
  randy[o] <- adjustedRandIndex(c.true,as.factor(c))
  print(randy[o])
  rmse[o] <- calcrmse(time.real,That)$rmse
  print(rmse[o])
  likli[o] <- loglikelihood(c,Y,mu,S,alpha,That, beta0, betahat, sigma2, lambda2, tau2, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.dat)
  print(likli[o])
  print(o/iter.burnin)
} 

############## GIBBS SAMPLING WITH THINNING ######################################################
nmrse <- c(0)
mu.list <- list(0)
beta0.list <- list(0)
betahat.list <- list(0) 
sigma2.list <- list(0)
lambda2.list <- list(0)
tau2.list <- list(0)
c.list <- list(0)
That.list <- list(0)
iter = 100
iter.thin =2
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
  
  
  
  
  
  ######################## The Censored Times ###########################################################
  source('updatetime.R')
  # Updating the Time Variable
  ti <- NA
  ti <- updatetime(c, Y, Time,That, beta0, betahat, sigma2)
  That <- ti$time
  
  
  ################## PARAMETERS OF THE DP Mixture Model ######################################################
  ## Updating the parameters based on the observations 
  
  if(o%% iter.thin == 0 ){
    mu.list[[count]] <- mu
    beta0.list[[count]] <- beta0
    betahat.list[[count]] <- betahat  
    sigma2.list[[count]] <- sigma2
    lambda2.list[[count]] <- lambda2
    tau2.list[[count]] <- tau2
    c.list[[count]] <- c
    That.list[[count]] <- That
    time.predicted <- predicttime(c,Y, That,Time,beta0, betahat, sigma2)$predicttime
    nmrse[count] <- calcrmse(time.real,time.predicted)$rmse
    count <- count +1
  }
  
  
  
  
  print(o/iter) 
  #   print(loglike[o])
  #   print(cindex)
} 

########## ANLAYSING THE OUTPUT #######################################################

count <- count -1

c.matrix <- matrix(NA, nrow = N, ncol = count)
for ( i in 1:count){
  c.matrix[,i] <- c.list[[i]]
}
c.final <- apply(c.matrix,1,median)

surv.ob <- Surv(time.real,censoring)
logrank <- survdiff(surv.ob ~ c.final)

list.betahat <- list(0)

for ( i in 1:count){
  list.betahat[[i]] <- (betahat.list[[i]][1:2,] != 0) +0
}

betahat1.final <- matrix(NA, nrow = count, ncol = D)
for ( i in 1:count){
  betahat1.final[i,] <- list.betahat[[i]][1,]
}
betahat2.final <- matrix(NA, nrow = count, ncol = D)
for ( i in 1:count){
  betahat2.final[i,] <- list.betahat[[i]][2,]
}
### Final bethats
final.betahat1 <- apply(betahat1.final,2,mean)
final.betahat2 <- apply(betahat2.final,2,mean)


### Probability of betahat of genes
final.betahat <- rbind(final.betahat1, final.betahat2)
rownames(final.betahat) = c("cluster_1","cluster_2")

colnames(final.betahat) =  c(rep("relevant",rel.D),rep("irrelevant",irrel.D))

pc <- prcomp(Y)
pc.pred <- predict(pc,newdata = Y)
plot(pc.pred[,1], pc.pred[,2], pch = 19,col = c.final)


heatmapdata <- as.data.frame(final.betahat)
pdf("/home/bit/ashar/ExpressionSets/Simulations/featureselection.pdf")
heatmap.2(t(as.matrix(heatmapdata)),dendrogram="none", col =cm.colors(180), margins=c(6,10), main = "Posterior probabilities \n for Selection \n in 1 Simulation ", cexCol = 0.85, cexRow = 0.7, Rowv = FALSE)
dev.off()
