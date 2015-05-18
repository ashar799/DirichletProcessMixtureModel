## Script for a 2 dimensional case Vizualization

rm(list = ls())
setwd('/home/bit/ashar/Dropbox/Code/DPmixturemodel/DPplusAFT')

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


#################################### SIMULATED DATA PROPERTIES ####################################################
## Number of points
N = 100

## Number of Clusters
F = 2

## Distribution of the points within three clusters

p.dist = c(0.6,0.4)

## Total Number of features D

D = 2

## Total Percentage of irrelevant feature

prob.noise.feature = 0.90

## Total Percentage of censoring

prob.censoring = 0.05


## Overlap between Cluster of molecular Data of the relevant features

prob.overlap = 0.05

## Percentage of Noise/Overlap in Time Data

prob.noise = 0.05

## Actual Number of Components and dimension  which are relevant
rel.D = as.integer(D* (1-prob.noise.feature))
## Actual Number of Irrelevant Componenets
irrel.D = D - rel.D


########################################Simulated Data##############################################
A <- MixSim(MaxOmega = prob.overlap ,K = F, p = rel.D, int =c(0,10), lim = 1e08)


data.mu = array(data = NA, dim =c(F,D))
data.S = array(data = NA, dim =c(F,D,D))

for( i in 1:F){
  data.mu[i,1:rel.D] <- A$Mu[i,1:rel.D]
  data.S[i,1:rel.D,1:rel.D] <- A$S[1:rel.D,1:rel.D,i]
}

## Initializing the points
# Y.list <- list(0)
# for ( i in 1:F){
# Y.list[[i]] <- matrix(0, nrow =  as.integer(N * p.dist[i]), ncol = D)
# }

# ## The relevant data is genereated first
# Y.list <- list(0)
# for ( i in 1:F){
#   Y.list[[i]] <- mvrnorm(n = as.integer(N * p.dist[i]), mu = data.mu[i,1:D], Sigma = diag(x =1, nrow = D, ncol = D))
# }
# 

## The relevant data is genereated first
Y.rel.list <- list(0)
for ( i in 1:F){
  Y.rel.list[[i]] <- mvrnorm(n = as.integer(N * p.dist[i]), mu = data.mu[i,1:rel.D], Sigma = data.S[i,1:rel.D,1:rel.D])
}

## Scaling the Data as ONLY the scaled data will be used for generating the times
Y.rel.sc.list <- list(0)
for ( i in 1:F){
  Y.rel.sc.list[[i]] <- scale(Y.rel.list[[i]], center = TRUE, scale = TRUE)
}

## Irrelevant features
Y.irrel.list <- list(0)
for ( i in 1:F){
  mean <- runif(irrel.D,0,10)
  Y.irrel.list[[i]] <- mvrnorm(n = as.integer(N * p.dist[i]), mu = mean, Sigma = diag(x =1, nrow = irrel.D, ncol = irrel.D))
}


## The Co-efficients have to be obtained from uniform distribution between [1,10]
beta.list <- list(0)
for ( i in 1:F){
  beta.list[[i]] <- runif(rel.D, min = -5, max = 5)
}

# order.rel <- list(0)
# for (i in 1:F){
# order.rel[[i]] <- sample(D, rel.D)  
# }


## The relevant data is 
#  Y.rel.list <- list(0)
# for ( i in 1:F){
#   Y.rel.list[[i]] <- Y.list[[i]][,order.rel[[i]]]
# }

# ## Scaling the Data as ONLY the scaled data will be used for generating the times
# Y.rel.sc.list <- list(0)
# for ( i in 1:F){
#   Y.rel.sc.list[[i]] <- scale(Y.rel.list[[i]], center = TRUE, scale = TRUE)
# }
# 
# ## Irrelevant feature
# Y.irrel.list <- list(0)
# for ( i in 1:F){
#   Y.irrel.list[[i]] <- Y.list[[i]][,-order.rel[[i]]]
# }


## The pure time is generated
time.pur.list <- list(0)
for ( i in 1:F){
  time.pur.list[[i]] <- t(beta.list[[i]]) %*% t(Y.rel.sc.list[[i]])
}

## Simulating Time Data which is ONE dimensional
time.cluster <- MixSim(MaxOmega = prob.noise, K = F, p = 1, int =c(0,10))

time.noise.list <- list(0)
for ( i in 1:F){
  time.noise.list[[i]] <- rnorm(as.integer(N * p.dist[i]), mean = time.cluster$Mu[i], sd = sqrt((time.cluster$S[i])))
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

## Data is shuffled with respect to the features
data.schuffled <- list(0)
for (i in 1:F){
  data.schuffled[[i]] <-  data.plain[[i]][,order.list[[i]]]
}



## True Labels for the points
c.true <- c(0)
for ( i in 1:F){
  c.true <- rbind(as.matrix(c.true) , as.matrix(c(rep(i, as.integer(N * p.dist[i])))))  
}
c.true <- as.factor(c.true[-1,])



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

### Making permutations also with order of the points
order.points <- sample(N)



## Adding the permuted order of points
Y <- Y[order.points,]
c.true <- c.true[order.points]
time <- time[order.points]

pdf('/home/bit/ashar/Dropbox/WarsawTalk/2DGeneExpression.pdf')
plot(Y[,1],Y[,2],col =c("blue", "green")[c.true], main = "2D Gene Expression", pch =19)
leg.text <- c("Rel.Vertical ","Blue 60","Green 40")
legend("bottomleft", leg.text, pt.cex = 0.5)
dev.off()



## Draw Corresponding Kaplan Meier Curves
surv.ob <- Surv(time,censoring)
surv.fit <- survfit(surv.ob ~ c.true)
logrank <- survdiff(surv.ob ~ c.true)


pdf('/home/bit/ashar/Dropbox/WarsawTalk/2DKaplanMeierCurve.pdf')
plot(surv.fit, col = c("blue", "green"))
title("Kaplan-Meier Curves\nfor the 2 Clusters")
leg.text <- c("p-value 1.63e-02")
legend("topright", leg.text)
dev.off()




