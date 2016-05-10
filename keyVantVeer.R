rm(list =ls())
############################################
library(penalized)
data(nki70)

Y <- as.matrix(nki70[,8:77])
time <- log(nki70[,1]* 30)
event <- nki70[,2]
censoring <- event

N <- nrow(Y)
D <- ncol(Y)
############################# PARAMETERS for GIBB's SAMPLING ####
iter = 100
iter.burnin = 100
iter.thin  = 5
k = 2


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


######## Predict on New Data Set #####################################
source('predictCLASS.R')
predictCLASS(Y.new, time.new)
## Check the predicted Rand Index 
print(posteriorprob)
test.randindex <- adjustedRandIndex(apply(posteriorprob,1,which.max),c.true.new)

######## Predict on New Data Set  BASED ON JUST THE MOLECULAR DATA #####################################
source('predictCLUSTER.R')
predictCLUSTER(Y.new)
## Check the predicted Rand Index 
print(posteriorprobMOL)
i.test.randindex <- adjustedRandIndex(apply(posteriorprobMOL,1,which.max),c.true.new)


source('predictchineseAFTtime.R')
predictchineseAFTtime(Y.new)
### Check of the Predicted C-index 
predicted.cindex <- survConcordance(Surv(exp(time.new),censoring.new) ~ exp(-post.time.avg))[1]

########### Check Prediction Ground Truth
source('predictionGroundTruth.R')
predictionGroundTruth()
