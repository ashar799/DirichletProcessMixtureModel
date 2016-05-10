##### This script CLUSTERS Periphery Cells ON THE PvsZ signature and PFS survival  ####################
library('xlsx')
library('limma')
rm(list =ls())
setwd("/home/abidata/Dinis/VIP/signatureANNO")

## load data ##

load('/home/abidata/Dinis/VIP/signatureANNO/Data/ExpressionConsole.normalised.RData')

## The phenoData
tab <- read.xlsx(file  = '/home/abidata/Dinis/VIP/Sample_info/14-12-16 Samples for classifier.xlsx', sheetIndex =1)

list.patients <- tab[2:186,4]


## Include ONLY those features that have correspoding annotation ~ 27148 out of 70000
## Include ONLY those Samples which are needed for classification 

exprs <- eset.ec.norm[featureNames(eset.ec.norm) %in% rownames(anno.GENENAME),sampleNames(eset.ec.norm) %in% list.patients]


##########################################################################################
######## SIGNATURE ####################################################################
##########################################################################################


## LIMMA PREFILTERING #################################################

pheno  <- pData(exprs) 

pheno$Case <- as.factor(as.matrix(pheno$Case))
pheno$Topo <- as.factor(as.matrix(pheno$Topo))
pheno$Class <- as.factor(as.matrix(pheno$Class))


mm <- model.matrix( ~ 0 + Topo + Case  + Class , pheno)

cc <- colnames(mm)

cc <-  gsub("_","", cc)

cc <- gsub("-","", cc)

colnames(mm) <- cc

cont.matrix <- makeContrasts(ZvsP = Topocenter - Topoperiphery, levels= mm)

fit <- lmFit(exprs, design = mm)
fit2 <- contrasts.fit(fit, cont.matrix)
fit3 <- eBayes(fit2)

DEG_table <- topTable(fit3, adjust="BH", coef='ZvsP', number = Inf)

## Only Those genes which have a high FC 1.4 and significance


list <-  rownames(DEG_table)[abs(DEG_table$logFC) > 1.2 & DEG_table$adj.P.Val < 0.01]
##############################################################################################
###############################################################################################
###### OUR SIGNATURE ###########################################################################

signature <- list

library('survival')
library('stats')


## load Phenodata ######
phen <- read.csv(file = '/home/abidata/Dinis/VIP/Sample_info/15-02-10 VIP Sample Collection Data_v4.3.csv', sep = ",", quote = '\"', header = TRUE)


### Getting those samples which are periphery
pheno.per <- phen[(phen$Topo == "periphery") & (!is.na(phen$Topo)),]

######################################################################################################
################# PROGRESSION FREE SURVIVAL ##########################################################
######################################################################################################



### Only those patients who have a NON NA is PFS

pheno.ready <- pheno.per[!is.na(pheno.per$PFS),]





### Load conversion matrix
anno <- read.csv('Data/final.ANNO.v34.csv', sep = ',')

#########################################################################################
########################################################################################
### EVEN THOUGH AFFYIDS MAY HAVE SOME ANNOTATIONS, 
### signature 
sig.anno <- anno[anno$ENTREZ.ID  %in% signature,]    

## Defining the data set
Xsig <- exprs(eset.ec.norm[featureNames(eset.ec.norm) %in% (sig.anno[,1]),])

### Making the data matrices ready for the COX PH
Xsig.ready <- subset(t(Xsig), subset = rownames(t(Xsig)) %in% pheno.ready$USI)


## Order the rows in the phenodata and the expression to be the same
ind1 <- match(rownames(Xsig.ready), pheno.ready$USI)

pheno.sig <- pheno.ready[ind1,]

## JUST CHECKING
pheno.sig$USI == rownames(Xsig.ready)


######## Getting survival times (PFS) and status #####################

time <- as.numeric(as.matrix(pheno.sig$PFS))
censoring <- 1- pheno.sig$Alive.


##### Applying the MODEL on CLUSTERING #####################################

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

###### Just Doing some changes to the data ##########################################3
time <- log(time*30)
Y <- Xsig.ready
index.inf <- which(time == -Inf)
time <- time[-index.inf]
Y <- Y[-index.inf,]
censoring <- censoring[-index.inf]
That <- time
## Let's See if the Clusters are separate on PCA

## All features but uncorrelated data
pc <- prcomp(Y)
pc.pred <- predict(pc,newdata = Y)
plot(pc.pred[,1], pc.pred[,2], pch = 19)



############################# PARAMETERS for GIBB's SAMPLING ######################################
iter.burnin = 100
iter.thin  =5

################################# GIBBS SAMPLING  ###################################################
Time <- c(0)
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

setwd('/home/bit/ashar/Dropbox/Code/DPmixturemodel/DPplusAFT')
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
F = 2

surv.ob <- Surv(exp(time),censoring)

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

library(cluster)
## Checking the significance of the different time curves and silhouette indices
dis <- dist(Y)
sil <- silhouette(c.init,dis)
## Checking the average silhouette index
sil.width0 <- unlist(summary(sil))$avg.width
logrank <- survdiff(surv.ob ~ c.init)
## Checkingthe P-value
logrank0 <- 1 - pchisq(unlist(logrank)$chisq, df=F-1)


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

timeparam <- NA
time.predicted <- c(0)


likli0 <- print(loglikelihood(c,Y,mu,S,alpha,That, beta0, betahat, sigma2, lambda2, tau2, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.dat) )
sili <- c(0)
plogi <- c(0)
o.initi <- o
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

surv.days <- Surv((30.5*time),censoring)
surv.days.fit <- survfit(surv.days ~ c.final)
surv.days.logrank <- survdiff(surv.days ~ c.final)






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
rownames(final.betahat) = c("cluster_red","cluster_black")

index.match <- match(colnames(Y),sig.anno[,1])
colnames(final.betahat) <- sig.anno[,3][index.match]



heatmapdata <- as.data.frame(final.betahat)
pdf("/home/bit/ashar/ExpressionSets/VIPdataset/SurvialRelatedSignature.pdf")
heatmap.2(t(as.matrix(heatmapdata)),dendrogram="none", col =cm.colors(180), margins=c(6,10), main = "Cluster \n specific \n Biomarkers ", cexCol = 0.85, cexRow = 0.7, trace = "none")
dev.off()

write.csv(final.betahat, file = 'betahatVIP.csv')


####### Generate HeatMasps ################################################3
subset.betahat <- which(apply(final.betahat,2,mean) != 0 )

## A new Data frame with only these pathways
Y.new <- Y[,subset.betahat]
colnames(Y.new) <- names(subset.betahat)
rownames(Y.new) <- c.final

ind1 <- which(c.final ==1)
ind2 <- which(c.final == 2)


Y.order <- matrix(0, nrow = N, ncol = ncol(Y.new))
Y.order <- rbind(Y.new[ind1,],Y.new[ind2,])


pdf("/home/bit/ashar/ExpressionSets/VIPdataset/ExpressionHeatmap.pdf")
heatmap.2(t(Y.order), dendrogram = "none",col = , Rowv = FALSE, Colv = FALSE, ColSideColors = rownames(Y.order), cexRow =0.6, margins = c(5,12), trace = "none", main = "Periphery Cells \n Expression level  with \n survival related genes") 
dev.off()

######################### CHANGES TO THE PLOT AFTER MEETING OF 8th September #############################################################
######################## I was asked to redo the survival plots WITHOUT LOGRANK TEST BUT JUST DIFFERENCES IN THE MEDIAN SURVIVAL ##########
## Let's clear up the memory first#########
rm(list = ls())
load("/home/bit/ashar/ExpressionSets/VIPdataset/VIPDP.RData")
### The pheno.sig has information about the patient pairing
patient.cases <- as.numeric(levels(as.factor(as.numeric(as.matrix(pheno.sig$Case)))))


### Let's create a vector which will store vectors which store indices which I will include in my survival plot

## Inde just has indices for each unique patient 
inde <- list(0)
for ( i in 1:length(patient.cases)){
  inde[[i]] <- which(pheno.sig$Case == patient.cases[i])
 }
## check will have indices for the patients that are good
check <- list(0)
for ( i in 1:length(inde)){
  if (length(unique(c.final[inde[[i]]])) == 1)
    {
    check[[i]] <- inde[[i]] 
    }else {
    check[[i]] <- NA
  }
}
## saving good indices
d<- unlist(check)
good.indices <- d[!is.na(d)]

### Subsetting all the necessary objects
time.subset <- time[good.indices]
censoring.subset <- censoring[good.indices]
sub.c.final <- c.final[good.indices]

sub.surv.days <- Surv((30.5*time.subset),censoring.subset)
sub.surv.days.fit <- survfit(sub.surv.days ~ sub.c.final)

pdf('/home/bit/ashar/ExpressionSets/VIPdataset/VIP_PCells_clusters.pdf')
par(mfrow=c(1,2))
plot(pc.pred[,1], pc.pred[,2], pch = 19,col = c.final, main = "P-Cells", xlab = "PC1", ylab = "PC2")
plot(sub.surv.days.fit, col = c("black", "red"), main = "Kaplan Meier Curves", xlab = " Median PFS black 262d \n Median PFS Red 204d")
dev.off()

