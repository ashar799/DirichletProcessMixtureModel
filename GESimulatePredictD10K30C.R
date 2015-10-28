### This simulation Checks if the model can predict well with 10 Dimensions
### Also CHECK THE TIME REQUIRED FOR THE MODEL


source('importDP.R')

### Remove The Past 
rm(list = ls())
#################################### SIMULATED DATA PROPERTIES ####################################################
## Number of points
N.test = 100
N.train =100

N = N.test 

## Number of Clusters
F = 3

## Distribution of the points within three clusters

p.dist = c(0.4,0.3,0.3)

## Total Number of features D

D = 10

## Total Percentage of irrelevant feature

prob.noise.feature = 0.6

## Total Percentage of censoring

prob.censoring = 0.10

## Overlap between Cluster of molecular Data of the relevant features

prob.overlap = 0.01

## Percentage of Noise/Overlap in Time Data

prob.noise = 0.1


############################# PARAMETERS for GIBB's SAMPLING ####
iter = 100
iter.burnin = 60
iter.thin  =5
k =3

########################### Initialize the functions ############

source('simulateDPMM.R')
simulateDPMM()


source('initializeDPMM.R')
initializeDPMM()



######### Ground Truth ##############################
source('SIMgroundtruth.R')
SIMgroundtruth()



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

source('predictchineseAFTtime.R')
predictchineseAFTtime(Y.new)


