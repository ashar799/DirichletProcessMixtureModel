## Testing the DPplusAFT model with Gibb's Sampling
### Script to test if the model works with higher dimensions
## D = 100
## K =3 (No of Cluster)
## With irrelevant features and with Time Censoring 

rm(list = ls())
setwd('/home/bit/ashar/Dropbox/Code/DPmixturemodel/DPplusAFT')
library(MASS)
library(mvtnorm)
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


#################################### SIMULATED DATA PROPERTIES ####################################################
## Number of points
N = 120

## Number of Clusters
F =3

## Distribution of the points within three clusters

p.dist = c(0.3,0.4,0.3)

## Total Number of features

D =30

## Total Percentage of irrelevant feature

prob.noise.feature = 0.80

## Total Percentage of censoring

prob.censoring = 0.05


## Overlap between Cluster

prob.overlap = 0.10

## Percentage of Noise/Overlap in Time Data

prob.noise = 0.10

## Actual Number of Components and dimension and dimension which are relevant
rel.D = as.integer(D* (1-prob.noise.feature))
## Actual Number of Irrelevant Componenets
irrel.D = D - rel.D


########################################Simulated Data##############################################
A <- MixSim(MaxOmega = prob.overlap ,K = F, p = rel.D, int =c(0,100))


data.mu = array(data = NA, dim =c(F,D))
data.S = array(data = NA, dim =c(F,D,D))

for( i in 1:F){
  data.mu[i,1:rel.D] <- A$Mu[i,1:rel.D]
  data.S[i,1:rel.D,1:rel.D] <- A$S[1:rel.D,1:rel.D,i]
}

## Initializing the points
Y.list <- list(0)
for ( i in 1:F){
Y.list[[i]] <- matrix(0, nrow =  as.integer(N * p.dist[i]), ncol = D)
}

## The relevant data is genereated first
Y.rel.list <- list(0)
for ( i in 1:F){
  Y.rel.list[[i]] <- mvrnorm(n = as.integer(N * p.dist[i]), mu = data.mu[i,1:rel.D], Sigma = data.S[i,1:rel.D,1:rel.D])
  Y.rel.list[[i]] <- scale(Y.rel.list[[i]])
}

## Generating irrelevant data
Y.irrel.list <- list(0)
for ( i in 1:F){
  Y.irrel.list[[i]] <- mvrnorm(n = as.integer(N * p.dist[i]), mu = rep(0, irrel.D), Sigma = diag(irrel.D))
}

## The Co-efficients have to be obtained from uniform distribution between [1,10]
beta.list <- list(0)
for ( i in 1:F){
beta.list[[i]] <- runif(rel.D, min =1, max = 10)
}


## The pure time is generated
time.pur.list <- list(0)
for ( i in 1:F){
  time.pur.list[[i]] <- t(beta.list[[i]]) %*% t(Y.rel.list[[i]])
}

## Simulating Time Data which is ONE dimensional
time.cluster <- MixSim(MaxOmega = prob.noise, K = F, p = 1, int =c(0,10))

time.noise.list <- list(0)
for ( i in 1:F){
  time.noise.list[[i]] <- rnorm(as.integer(N * p.dist[i]), mean = time.cluster$Mu[i], sd = sqrt(time.cluster$S[i]))
}

time.list <- list(0)
for ( i in 1:F){
  time.list[[i]] <- time.pur.list[[i]] + time.noise.list[[i]]
}

### Combining the data with relevant and irrelevant columns

data.plain <- list(0) 
for (i in 1:F){
 data.plain[[i]] <-  cbind(Y.rel.list[[i]], Y.irrel.list[[i]]) 
}


### Making permutations with the order of the features
order.list <- list(0)
for (i in 1:F){
  order.list[[i]] <-  sample(D)
}

## Data is shuffled
data.schuffled <- list(0)
for (i in 1:F){
  data.schuffled[[i]] <-  data.plain[[i]][,order.list[[i]]]
}



############################################### MAKING Y from the clusters data #####################3
Y <- c(0)
for (i in 1:F){
  Y <- rbind(Y, data.schuffled[[i]])
} 
Y <- Y[-1,]

#######################################MAKING TIME from cluster data ########################################################
time <- c(0)
for (i in 1:F){
  time <- cbind(time, time.list[[i]])
} 
time <- time[,-1]
time <- as.vector(time)


####################################### Adding CENSORING INFORMATION  ################################################################
## Adding the censoring information

censoring <- rbinom(n = NROW(Y), size =1, prob = 1- prob.censoring)
right.censoring.time <- min(time)  


index.time <- which(censoring==0)
for ( q in 1:length(index.time)){
  time[index.time[q]] <- right.censoring.time
  
}

############################# PARAMETERS for GIBB's SAMPLING ######################################

iter = 500
iter.samples = 100

################################# CALLING MAIN FUNCTION ###################################################

result <- NA

source('dirichletmixture.R')
result <- dirichletmixture(Y, time, censoring, iter,F) 
