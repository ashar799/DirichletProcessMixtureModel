#### A script to simulate Gene Expression with some part of the columns being completely unrelated

rm(list = ls())

N = 100

D  = 500

rel.D = 20

irrel.D = 480

X <- matrix(runif(n =N*D, min = -1.5, max = 1.5), nrow = N, ncol = D)

rel.X <- as.matrix(X[,1:rel.D])

obj.qr <- qr(X)


alpha <- qr.Q(obj.qr)[,1:rel.D]

gamma <- qr.Q(obj.qr)[,(1+rel.D):N]

matT <- matrix(runif(n = rel.D*(N-rel.D), min = -0.005, max= 0.005), nrow = rel.D, ncol = (N-rel.D))

matP <- t(matT) %*% matT

max.eig <- eigen(matP)$values[1]

max.corr <- sqrt(max.eig)/sqrt(1 + max.eig)

linear.space <- gamma + alpha %*% matT

irrel.X <- matrix(NA, nrow = N, ncol = irrel.D)

for ( i in 1: irrel.D){
  
  matTemp <- matrix(runif(n = (N-rel.D), min = -1.5, max= 1.5),  nrow = (N-rel.D), ncol =1)
  irrel.X[,i] <- as.vector(linear.space %*% matTemp)

}

## Checking if the covariance is indeed small

cov.mat <- cov(rel.X,irrel.X)

boxplot(cov.mat)

## Building the full data matrix

X.full <- cbind(rel.X, irrel.X)


levelplot(cov(X.full[,1:100]))
#########################################################################################
##########################################################################################
##### Now WE DEAL WITH CLUSTERED DATA AND GENERATE NON RELEVANT FEATURES INDEPENENTLY ###

library(MixSim)
########################################Simulated Data##############################################
F =3

p.dist = c(0.3,0.4,0.3)


prob.overlap = 0.05

A <- MixSim(MaxOmega = prob.overlap ,K = F, p = rel.D, int =c(-1.5,1.5), lim = 1e08)


data.mu = array(data = NA, dim =c(F,D))
data.S = array(data = NA, dim =c(F,D,D))

for( i in 1:F){
  data.mu[i,1:rel.D] <- A$Mu[i,1:rel.D]
  data.S[i,1:rel.D,1:rel.D] <- A$S[1:rel.D,1:rel.D,i]
}

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

Y.rel <- c(0)
for (i in 1:F){
  Y.rel <- rbind(Y.rel, Y.rel.list[[i]])
} 
Y.rel <- Y.rel[-1,]

Y.irrel <- c(0)
for (i in 1:F){
  Y.irrel <- rbind(Y.irrel, Y.irrel.list[[i]])
} 
Y.irrel <- Y.irrel[-1,]

## Checking to see what is the correlation like

cov.mat.2 <- cov(Y.rel,Y.irrel)

boxplot(cov.mat.2)

Y.full <- cbind(Y.rel,Y.irrel)

levelplot(cov(Y.full[,1:100]))

## Checking the kind of lower dimensional space
## True Labels for the points
c.true <- c(0)
for ( i in 1:F){
  c.true <- rbind(as.matrix(c.true) , as.matrix(c(rep(i, as.integer(N * p.dist[i])))))  
}
c.true <- as.factor(c.true[-1,])



pc <- prcomp(Y.rel)
pc.pred <- predict(pc,newdata = Y.rel)
plot(pc.pred[,1], pc.pred[,2], pch = 19,col = c.true)


############################################################################
########### CLEARLY WE NEED TO COMBINE THE TWO TO HAVE CLUSTERING AS WELL AS CORRELATION
########### ONLY RELEVANT FEATURES SHOULD BE USED TO DEFINE THE OVERLAP #####################


## Selcting the beta co-efficients

## The Co-efficients have to be obtained from uniform distribution between [-3,3]

beta.list <- list(0)
for ( i in 1:F){
  beta.list[[i]] <- as.vector(rbind(runif(half, min = -3, max = -0.1), runif(half, min = 0.1, max = 3)))
}

#
