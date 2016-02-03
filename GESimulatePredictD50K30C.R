### This simulation Checks if the model can predict well with 10 Dimensions
### Also CHECK THE TIME REQUIRED FOR THE MODEL


### Remove The Past 
rm(list = ls())
#################################### SIMULATED DATA PROPERTIES ####################################################
## Number of points
N.test = 100
N.train =100

## Number of Clusters
F = 3

## Distribution of the points within three clusters

p.dist = c(0.4,0.3,0.3)

## Total Number of features D

D = 50

## Total Percentage of irrelevant feature

prob.noise.feature = 0.6


## Overlap between Cluster of molecular Data of the relevant features

prob.overlap = 0.01

###### Get the Data #####################################

## Initialize the Training Data
source('simulateDPMM.R')
simulateDPMM()



############################# PARAMETERS for GIBB's SAMPLING ####
iter = 100
iter.burnin = 5*D
iter.thin  =5
k =3

########################### Initialize the functions ############



######################### Initialize the Parameters ################
source('initializeDPMM.R')
initializeDPMM()


######### Ground Truth ##############################
source('SIMgroundtruth.R')
#SIMgroundtruth()



########### Train the Model #########################################
source('burninDPMM.R')
burninDPMM()


source('gibbsDPMM.R')
gibbsDPMM()

########## Analyze the fit ##########################################
source('SIManalyzeDPMM.R')
SIManalyzeDPMM()


######## Predict on New Data Set #####################################
source('predictCLASS.R')
predictCLASS(Y.new, time.new)
## Check how much concordance is there
test.randindex <- adjustedRandIndex(apply(posteriorprob,1,which.max),c.true.new)

source('predictchineseAFTtime.R')
predictchineseAFTtime(Y.new)

### Check of the Predicted C-index ## 91% concordance
predicted.cindex <- survConcordance(Surv(time.new,censoring.new) ~ exp(-post.time.avg))[1]


