### This simulation Checks if the model can predict well with 10 Dimensions
### Also CHECK THE TIME REQUIRED FOR THE MODEL
### Remove The Past 
rm(list = ls())

RIground <- c(0)
CIground <- c(0)


for ( v in 1:20){
  

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
prob.noise.feature = 0.5
  
  
## Overlap between Cluster of molecular Data of the relevant features
prob.overlap = 0.10
#################################### SIMULATED DATA PROPERTIES ####################################################

############################# PARAMETERS for GIBB's SAMPLING ####
iter = 100
iter.burnin = 100
iter.thin  = 5
k = 2
F = k


## Initialize the Training Data
source('simulateDPMM.R')
simulateDPMM()

######### OR Don't preprocess use this ########3
Y <- Y.dat
Y.new <- Y.new.dat


########### Check Prediction Ground Truth
source('predictionGroundTruth.R')
predictionGroundTruth()

RIground[v] <- predRandIndex
CIground[v] <- predCIndex

}

