## Testing the DPplusAFT model with Gibb's Sampling
### Script to test if the model works with higher dimensions
## With irrelevant features and with Time Censoring 

rm(list = ls())
setwd('/home/bit/ashar/Dropbox/Code/DPmixturemodel/DPplusAFT')
load(file = '/home/bit/ashar/ExpressionSets/Simulations/Simulations.RData')
library(rms)
library(pec)
library(ipred)
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
library(survival)
############################################################################################
######### LOG LIKELIHOOD ###################################################################

### Getting cluster assignments based on the K-means clustering
gr.km <- kmeans(Y, 2, nstart =10)
adjustedRandIndex(c.true,as.factor(gr.km$cluster))

brier.km <- c(0)
brier.aft <- c(0)
brier.cox <- c(0)

for (i in 1:F){
  
  ### Cluster Specific Data
  
  ind <- which((gr.km$cluster) == i)
  time.tmp <- time[ind]
  censoring.tmp <- censoring[ind]
  Y.tmp <- Y[ind,]
  rownames(Y.tmp) <- as.character(c(1:nrow(Y.tmp)))
  
  ########## PART FOR THE KAPLAN MEIER CURVE ESTIMATE ##################
  
  smod <-  Surv(exp(time.tmp), censoring.tmp)
  fit.km <- survfit(smod ~ 1)
  brier.km[i] <- crps(pec(fit.km, formula =  Surv(exp(time.tmp), censoring.tmp) ~ 1, data = as.data.frame(Y.tmp), times = exp(time.tmp)))[1]
  
  
  ########## PART FOR THE AFT ###########################################
  L = length(ind)
  f1 <- psm(smod ~ Y.tmp , dist="lognormal")
  S1 <- Survival(f1)
  mat.tmp <- matrix(NA, nrow = L, ncol = L)
  for (j in 1:L){
    mat.tmp[,j] <- S1(exp(time.tmp[j]),f1$linear.predictors)
    
  }
  
  #### BRIER SCORE FOR THE AFT ######################################
  brier.aft[i] <- sbrier(smod,mat.tmp,exp(time.tmp))[1]
  
  
  ########### PART FOR THE COXPH ##########################################
  f2 <-   coxph(smod ~ Y.tmp, data = as.data.frame(Y.tmp))
  fit.coxph <- survfit(f2, newdata = as.data.frame(Y.tmp))
  
  
  #### BRIER SCORE FOR COX PH
  brier.cox[i] <- crps(pec(list("CoxPH"=f2),data= as.data.frame(Y.tmp),formula=Surv(exp(time.tmp), censoring.tmp) ~ Y.tmp))[[2]]
  
  ############# Holger's code ############################################
  
  
  
}

### The above algorithm do not necessarily work when High Dimensional Data is concerned
### Changes the code to work with high dimensional data 
brier.reg.aft <- c(0)
brier.reg.cox <- c(0)

for ( i in 1:F){
  
ind <- which((gr.km$cluster) == i)
time.tmp <- time[ind]
censoring.tmp <- censoring[ind]
Y.tmp <- Y[ind,]
rownames(Y.tmp) <- as.character(c(1:nrow(Y.tmp)))

########## PART WITH DECLARING THE .... ##################

smod <-  Surv(exp(time.tmp), censoring.tmp)
L = length(ind)

########## PART FOR THE AFT ###########################################



#### Using GLMpath to calculate regularized Cox Model
library(glmpath)
coxreg <- list(0)
coxreg$x <- Y.tmp
coxreg$time <- exp(time.tmp)
coxreg$status <- censoring.tmp
cv.path <- cv.coxpath(data = coxreg, mode = "lambda" )
f.reg <- predict(object = path, data = coxreg, s = 10, type =  "coxph", mode = "lambda")
fit.coxregph <-  survfit(f.reg)
brier.reg.cox[i] <- crps(pec(list("CoxPH"=f.reg),data= as.data.frame(Y.tmp),formula=Surv(exp(time.tmp), censoring.tmp) ~ Y.tmp))[[2]]

### For the AFT Model I should use glmnet
library(glmnet)
reg <- cv.glmnet(x = Y.tmp, y = time.tmp, family = "gaussian")
linear.pred <- predict(object =reg, newx = Y.tmp, s= "lambda.min")
coeff.pred <- coef(object =reg, newx = Y.tmp, s= "lambda.min")
rel.coeff <- coeff.pred[2:(D+1)] 
ind.rel <- which(rel.coeff !=0)
### Calculate the error variance as it is required to calcuate the Survival Probabilities for AFT model
predicted.tmp <- predict(object = reg, newx = Y.tmp, s = "lambda.min") 
sum.err <-c(0)
for ( j in 1:L){
 sum.err[j] <- (time.tmp[j] - predicted.tmp[j])^2 
}
variance = sqrt(median(sum.err))

## This is just to get a survival function1 - pnorm((t.trans - lp)/exp(parms))
aft_survival = function (times = NULL, lp = NULL, parms = variance) 
{
  t.trans <- logb(times)
  names(t.trans) <- format(times)
  1 - pnorm((t.trans - lp)/(parms))
}

mat.reg.tmp <- matrix(NA, nrow = L, ncol = L)
for (j in 1:L){
  mat.reg.tmp[,j] <- aft_survival(exp(time.tmp[j]),predicted.tmp)
}
brier.reg.aft[i] <- sbrier(smod,mat.reg.tmp,exp(time.tmp))[1]
}


