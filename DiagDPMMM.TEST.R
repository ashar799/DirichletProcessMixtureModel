## Script for a 2D visualization using a 4D example
### Remove The Past 
rm(list = ls())
#################################### SIMULATED DATA PROPERTIES ####################################################
## Number of points
N.test = 100
N.train = 100

## Number of Clusters
F = 2

## Distribution of the points within three clusters

p.dist = c(0.4,0.6)

## Total Number of features D
D = 4

## Total Percentage of irrelevant feature FOR THE SURVIVAL TIMES
prob.noise.feature = 0.5


## Overlap between Cluster of molecular Data of the relevant features
prob.overlap = 0.01

###### Get the Data #####################################
## Initialize the Training Data
source('simulateDPMM.R')
simulateDPMM()



############################# PARAMETERS for GIBB's SAMPLING ####
iter = D*50
iter.burnin = D*50
iter.thin  =5
k =2


######################### Initialize the Parameters ################
source('initializeDPMM.R')
initializeDPMM()



########### Train the Model #########################################
source('burninDPMM.R')
burninDPMM()


source('gibbsDPMM.R')
gibbsDPMM()

########## Analyze the fit ##########################################
source('SIManalyzeDPMM.R')
SIManalyzeDPMM()

pdf('/home/bit/ashar/Dropbox/WarsawTalk/2DGeneExpression.pdf')
plot(Y[,3],Y[,4],col =c("blue", "green","orange")[c.final], main = "2D Gene Expression", pch =19)
leg.text <- c("Rel.Vertical ","Blue 60","Green 40")
legend("bottomleft", leg.text, pt.cex = 0.5)
dev.off()



## Draw Corresponding Kaplan Meier Curves
surv.ob <- Surv(time,censoring)
surv.fit <- survfit(surv.ob ~ c.final)
logrank <- survdiff(surv.ob ~ c.final)


pdf('/home/bit/ashar/Dropbox/WarsawTalk/2DKaplanMeierCurve.pdf')
plot(surv.fit, col = c("blue", "green","orange"))
title("Kaplan-Meier Curves\nfor the 2 Clusters")
leg.text <- c("p-value 1.63e-02")
legend("topright", leg.text)
dev.off()

