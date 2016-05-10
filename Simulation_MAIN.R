### This simulation Checks if the model can predict well with 10 Dimensions
### Also CHECK THE TIME REQUIRED FOR THE MODEL
### Remove The Past 

rm(list = ls())
#################################### SIMULATED DATA PROPERTIES ####################################################
## Number of points
N.test = 100
N.train = 100

## Number of Clusters
F = 2

## Distribution of the points within three clusters

p.dist = c(0.5,0.5)

## Total Number of features D

D = 20

## Total Percentage of irrelevant feature
prob.noise.feature = 0.2


## Overlap between Cluster of molecular Data of the relevant features
prob.overlap = 0.01

###### Get the Data #####################################

## Initialize the Training Data
source('simulateDPMM.R')
simulateDPMM()

####### Pre process the data ###############
#source('preprocessDPMM.R')
#preprocessDPMM()

###### Use a Penalized LASSO to get to know which features are relevant##############
# smod <-  Surv(exp(time), censoring)
# reg.pcox <- cv.glmnet(x = Y.dat, y = time)
# ind.rel <- unlist(predict(object =reg.pcox, newx = Y, s = mean(reg.pcox$lambda[28:30]),type = "nonzero"))
# Y <- Y.dat[,ind.rel]
# Y.new <- Y.new.dat[,ind.rel]
# D <- ncol(Y) 
# pc <- prcomp(Y)
# pc.pred <- predict(pc,newdata = Y)
# plot(pc.pred[,1], pc.pred[,2], pch = 19, col =c.true, main = 'Main Training Data After Preprocessing')

######### OR Don't preprocess use this ########3
Y <- Y.dat
Y.new <- Y.new.dat

############################# PARAMETERS for GIBB's SAMPLING ####
iter = 100
iter.burnin = 50
iter.thin  = 5
k = 2


######################### Initialize the Parameters ################
source('initializeDPMM.R')
initializeDPMM()


######### Ground Truth ##############################
source('SIMgroundtruth.R')
SIMgroundtruth()

source('iSIMgroundtruth.R')
iSIMgroundtruth()

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

######### Use PreMiuM ##############3
source('premium.R')
direc <- as.character('/home/bit/ashar/ownCloud/Research/DPMMSIMULATIONS/OneView/D20Noise20perOverlap01per/premium')
premium(direc)


