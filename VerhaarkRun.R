### Script to test if the model works with the Verhaak Data


rm(list = ls())
load('/home/bit/ashar/ExpressionSets/Verhark/30pathwayVerhaark.RData')
load('/home/bit/ashar/ExpressionSets/Verhark/phenoVerhaark.RData')
setwd('/home/bit/ashar/Dropbox/Code/DPmixturemodel/DPplusAFT')

### Load Results from the Last time
##load('/home/bit/ashar/ExpressionSets/Verhark/Verhark.RData')

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


Y <- data.combined
time <- pheno[,3]
censoring <- pheno[,2]
kruskal.test(time, pheno[,4])
c.true <- pheno[,4]
## p-value 0.851 This means that survival times are not significantly different YEEEEY!!!!


## Let's See if the Clusters are separate on PCA

## All features but uncorrelated data
pc <- prcomp(Y)
pc.pred <- predict(pc,newdata = Y)
plot(pc.pred[,1], pc.pred[,2], pch = 19,col = c.true)

############################# PARAMETERS for GIBB's SAMPLING ######################################

iter = 10000
iter.burnin = 1000
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

# ## Check What is the RMSE
# calcrmse(time.real,That)$rmse
F = 4

surv.ob <- Surv(time,censoring)

## Initialization part for the parmaters of AFT Model with k-means and Bayesian Lasso
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

### Checking the Plot with new Clustering
plot(pc.pred[,1], pc.pred[,2], pch = 19,col = c.init)
## Checking the rand index
adjustedRandIndex(c.true,as.factor(c.init))
## So only 0.45 matching

## Checking the significance of the different time curves and silhouette indices
dis <- dist(Y)
sil <- silhouette(c.init,dis)
## Checking the average silhouette index
unlist(summary(sil))$avg.width
logrank <- survdiff(surv.ob ~ c.init)
## Checkingthe P-value
1 - pchisq(unlist(logrank)$chisq, df=3)


# Testing the  k-means estimate
source('predicttime.R')
time.predicted <- predicttime(c,Y, That,Time,beta0, betahat, sigma2)$predicttime



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
randy <- c(0)
likli <- c(0)
sili <- c(0)
plogi <- c(0)

#################### BURNIN PHASE ###################################################
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
  
  
  #source('posterioralpha.R') 
  ## Updating the concentration parameter
  # alpha <- posterioralpha(c, N, alpha, shape.alpha, rate.alpha)
  
  
  
  
  
  ######################## The Censored Times ###########################################################
  source('updatetime.R')
  # Updating the Time Variable
  ti <- NA
  ti <- updatetime(c, Y, Time,That, beta0, betahat, sigma2)
  That <- ti$time
  
  
  ##################### Print SOME Statistics #####################################################
  randy[o] <- adjustedRandIndex(c.true,as.factor(c))
  print(randy[o])
  likli[o] <- loglikelihood(c,Y,mu,S,alpha,That, beta0, betahat, sigma2, lambda2, tau2, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.dat)
  print(likli[o])
  ## silhouette index
  ## Checking the average silhouette index
  sili[o] <- unlist(summary(silhouette(c,dis)))$avg.width
  print(sili[o])
  ## Checkingthe P-value
  numclust <- table(factor(c, levels = 1:K))
  activeclust <- which(numclust!=0)
  deg <- length(activeclust) -1
  plogi[o] <- 1 - pchisq(unlist(survdiff(surv.ob ~ c))$chisq, df=deg)
  print(plogi[o])
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

############## GIBBS SAMPLING WITH THINNING ######################################################
iter =100
thin = 2
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
  
  
  source('posterioralpha.R') 
  ## Updating the concentration parameter
  alpha <- posterioralpha(c, N, alpha, shape.alpha, rate.alpha)
  
  
  ######################## The Censored Times ###########################################################
  source('updatetime.R')
  # Updating the Time Variable
  ti <- NA
  ti <- updatetime(c, Y, Time,That, beta0, betahat, sigma2)
  That <- ti$time
  
  
  
  
  #   ## Value of the Log-likelihood
  source('likelihood.R')
  loglike[o] <-loglikelihood(c,Y,mu,S,alpha,That, beta0, betahat, sigma2, lambda2, tau2, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.dat) 
  
  if(o%% thin == 0 ){
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

#### Taking the average values for the results #########################################

c.matrix <- matrix(NA, nrow = N, ncol = count)
for ( i in 1:count){
  c.matrix[,i] <- c.list[[i]]
}
c.final <- apply(c.matrix,1,median)

surv.fit <- survfit(surv.ob ~ c.final)
logrank <- survdiff(surv.ob ~ c.final)


#### SURVIVAL OBJECTS BASED ON DAYS ########################################

surv.days <- Surv(exp(time),censoring)
surv.days.fit <- survfit(surv.days ~ c.final)
surv.days.logrank <- survdiff(surv.days ~ c.final)

ver.surv.days.fit <- survfit(surv.days ~ c.true)
ver.surv.days.logrank <- survdiff(surv.days ~ c.true)


## Calculate Plots for our Clustering
pdf('/home/bit/ashar/KaplanMeierVerhaakDPMODEL.pdf')
plot(surv.days.fit, col =c("red","blue","green","black"), xlab = "Days")
title("Kaplan-Meier Curves\n with clusters \n from our DP model ")
leg.text <- c("LogRank Test p-value of 0.129")
legend("topright", leg.text)
dev.off()

pdf('/home/bit/ashar/MolecularVerhaakOURCLUSTERING.pdf')
pc <- prcomp(Y)
pc.pred <- predict(pc,newdata = Y)
plot(pc.pred[,1], pc.pred[,2], pch = 19,col = c.final)
title("PCA plot for Verhaak Data \n from our DP model \n with Pathway features")
dev.off()

pdf('/home/bit/ashar/KaplanMeierVerhaakOriginalClustering.pdf')
plot(ver.surv.days.fit, col =c("red","blue","green","black"), xlab = "Days")
title("Kaplan-Meier Curves\n with clusters \n from original Verhaak classification ")
leg.text <- c("LogRank Test p-value of 0.952")
legend("topright", leg.text)
dev.off()


load(file = "/home/bit/ashar/ExpressionSets/Verhark/OriginalVerhaakData.RData")
class(exprs.norm) <- "numeric"

pdf('/home/bit/ashar/MolecularVerhaakOriginalCLustering.pdf')
pc <- prcomp(exprs.norm)
pc.pred <- predict(pc,newdata = as.matrix(exprs.norm))
plot(pc.pred[,1], pc.pred[,2], pch = 19,col = c.true)
title("PCA plot for Verhaak Data \n from original Verhaak classification \n with Genes as features")
dev.off()





list.betahat <- list(0)

for ( i in 1:count){
  list.betahat[[i]] <- (betahat.list[[i]][1:4,] != 0) +0
}

betahat1.final <- matrix(NA, nrow = count, ncol = D)
for ( i in 1:count){
  betahat1.final[i,] <- list.betahat[[i]][1,]
}
betahat2.final <- matrix(NA, nrow = count, ncol = D)
for ( i in 1:count){
  betahat2.final[i,] <- list.betahat[[i]][2,]
}
betahat3.final <- matrix(NA, nrow = count, ncol = D)
for ( i in 1:count){
  betahat3.final[i,] <- list.betahat[[i]][3,]
}
betahat4.final <- matrix(NA, nrow = count, ncol = D)
for ( i in 1:count){
  betahat4.final[i,] <- list.betahat[[i]][4,]
}

### Final bethats
final.betahat1 <- apply(betahat1.final,2,mean)
final.betahat2 <- apply(betahat2.final,2,mean)
final.betahat3 <- apply(betahat3.final, 2,mean)
final.betahat4 <- apply(betahat4.final, 2,mean)


### Probability of betahat of genes
final.betahat <- rbind(final.betahat1, final.betahat2, final.betahat3, final.betahat4)
rownames(final.betahat) = c("cluster_1","cluster_2","cluster_3","cluster_4")
yy <- as.list(KEGGPATHID2NAME)
colnames(final.betahat) <- as.character(unlist(yy[colnames(Y)])) 






#### Subset the betahat such that only non zero elements

subset.betahat <- which(apply(final.betahat,2,mean)!=0 )

## A new Data frame with only these pathways
Y.new <- Y[,subset.betahat]
colnames(Y.new) <- names(subset.betahat)
rownames(Y.new) <- c.final

ind1 <- which(c.final ==1)
ind2 <- which(c.final == 2)
ind3 <- which(c.final ==3)
ind4 <- which(c.final ==4)

Y.order <- matrix(0, nrow = N, ncol = ncol(Y.new))
Y.order <- rbind(Y.new[ind1,],Y.new[ind2,],Y.new[ind3,],Y.new[ind4,])


pdf("/home/bit/ashar/ExpressionSets/Verhark/PathwayHeatmap.pdf")
heatmap.2(t(Y.order), dendrogram = "none",col = , Rowv = FALSE, Colv = FALSE, ColSideColors = rownames(Y.order), cexRow =0.6, margins = c(5,12), trace = "none", main = "Pathway expression \n Verhaak data") 
dev.off()



heatmapdata <- as.data.frame(final.betahat)
pdf("/home/bit/ashar/ExpressionSets/Verhark/heatmap.pdf")
heatmap.2(t(as.matrix(heatmapdata)),dendrogram="none", col =cm.colors(90), margins=c(6,10), main = "Marginal Posterior prob. ", cexCol = 0.85, cexRow = 0.7, ylab ="KEGG pathways", trace = "none")
dev.off()


#### BAR PLOTS ##############
pdf("/home/bit/ashar/ExpressionSets/Verhark/Barplot.pdf")
par(mfrow=c(2,2))
barplot(final.betahat[1,], main = "Cluster1")
barplot(final.betahat[2,], main = "Cluster2")
barplot(final.betahat[3,], main = "Cluster3")
barplot(final.betahat[4,], main = "Cluster4")
dev.off()


surv.fit <- survfit(surv.ob ~ c.final)
logrank <- survdiff(surv.ob ~ c.final)
logrank_verhaak <- survdiff(surv.ob ~ c.true)



pdf('/home/bit/ashar/Dropbox/WarsawTalk/KaplanMeier.pdf')
plot(surv.fit, col = c("blue", "green"))
title("Kaplan-Meier Curves\nfor the Simulation")
leg.text <- c("LogRank Test p-value of 7.79e-09")
legend("topright", leg.text)
dev.off()
