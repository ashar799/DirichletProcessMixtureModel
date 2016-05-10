######### This file Runs the Model for Lung Adenocarcinoma ##########################
rm(list = ls())
################################################
######## Load the Data from the folder #########
load("/home/bit/ashar/ExpressionSets/ONE_VIEW/lungCancer/DataLungAdeno-only1Datset.RData")
obt <- scale(relev$Y.train, center = TRUE, scale = TRUE)
Y <- scale(relev$Y.train, center = attr(obt,"scaled:center"), scale = (attr(obt,"scaled:scale")))
Y.new <- scale(relev$Y.test, center = attr(obt,"scaled:center"), scale = attr(obt,"scaled:scale"))

Y <- impute.knn(Y)$data
Y.new <- impute.knn(Y.new)$data




time <- log(relev$time.train)
censoring <- relev$censoring.train

Y.new <- relev$Y.test
time.new <- log(relev$time.test)
censoring.new <- relev$censoring.test

c.true <-  relev$c.true
c.true.new <- relev$c.true.new


########## Fitting Of Univariate Cox Models to all the features
N <- nrow(Y)
D <- ncol(Y)
N.new <- nrow(Y.new)

############################# PARAMETERS for GIBB's SAMPLING ####
iter = 100
iter.burnin = 100
iter.thin  = 5
k = 3
Nps = 20

######################### Initialize the Parameters ################
source('initializeDPMM.R')
initializeDPMM()


########### Train the Model #########################################
source('burninDPMM.R')
burninDPMM()

source('gibbsDPMM.R')
gibbsDPMM()

########## Analyze the fit ##########################################
### Good feature selection from heatmap plus cindex plus randindex
source('SIManalyzeDPMM.R')
SIManalyzeDPMM()



##########Compare the Kaplan Meier Curves #############
####Fitting Survival Curves
surv.ob <- Surv(time,censoring)
surv.fit1 <- survfit(surv.ob ~ c.true)
surv.fit2 <- survfit(surv.ob ~ c.final)
logrank1 <- survdiff(surv.ob ~ c.true)
logrank2 <- survdiff(surv.ob ~ c.final)


######## Predict on New Data Set  BASED ON JUST THE MOLECULAR DATA #####################################
source('predictCLUSTER.R')
predictCLUSTER(Y.new)
## Check the predicted Rand Index 
print(posteriorprobMOL)
c.test <- apply(posteriorprobMOL,1,which.max)
i.test.randindex <- adjustedRandIndex(c.test,c.true.new)


source('predictchineseAFTtime.R')
predictchineseAFTtime(Y.new)
### Check of the Predicted C-index 
predicted.cindex <- survConcordance(Surv(exp(time.new),censoring.new) ~ exp(-post.time.avg))[1]


