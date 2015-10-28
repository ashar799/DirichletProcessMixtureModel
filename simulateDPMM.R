simulateDPMM = function(){
  

  ## Actual Number of Components and dimension  which are relevant
  rel.D  <- as.integer(D* (1-prob.noise.feature))
  ## Actual Number of Irrelevant Componenets
  irrel.D <-  D - rel.D
  
  
  
  
  ## The data
  A <- MixSim(BarOmega = prob.overlap ,K = F, p = rel.D, int =c(-1.0,1.0), lim = 1e09)
  data.mu = array(data = NA, dim =c(F,D))
  data.S = array(data = NA, dim =c(F,D,D))
 
## Data coming from a hypothetical population diagonal  precision matrix
 for( i in 1:F){
  data.mu[i,1:rel.D] <- A$Mu[i,1:rel.D]
  data.S[i,1:rel.D,1:rel.D] <- diag(diag(A$S[1:rel.D,1:rel.D,i]))
}
## Data coming from a hypothetical population NON diagoanl matrix
# for( i in 1:F){
#   data.mu[i,1:rel.D] <- A$Mu[i,1:rel.D]
#   data.S[i,1:rel.D,1:rel.D] <- A$S[1:rel.D,1:rel.D,i]
# }



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
  mean <- runif(irrel.D,-1.5,1.5)
  Y.irrel.list[[i]] <- mvrnorm(n = as.integer(N * p.dist[i]), mu = mean, Sigma = diag(x =1, nrow = irrel.D, ncol = irrel.D))
}

### Combining the data with relevant and irrelevant columns
data.old <- list(0) 
for (i in 1:F){
  data.old[[i]] <-  cbind(Y.rel.list[[i]], Y.irrel.list[[i]]) 
}

############################################### MAKING Y from the clusters data #####################3
Y.old <- c(0)
for (i in 1:F){
  Y.old <- rbind(Y.old, data.old[[i]])
} 
Y.old <- Y.old[-1,]

######################################################################################### Making the irrelevant features independent from the dependent features #############
X <- Y.old

rel.X <- as.matrix(X[,1:rel.D])

obj.qr <- qr(X)

rk <- obj.qr$rank

alpha <- qr.Q(obj.qr)[,1:rel.D]

gamma <- qr.Q(obj.qr)[,(1+rel.D):rk]

matT <- matrix(runif(n = rel.D*(rk -rel.D), min = -0.005, max= 0.005), nrow = rel.D, ncol = (rk -rel.D))

matP <- t(matT) %*% matT

max.eig <- eigen(matP)$values[1]

max.corr <- sqrt(max.eig)/sqrt(1 + max.eig)

linear.space <- gamma + alpha %*% matT

irrel.X <- matrix(NA, nrow = N, ncol = irrel.D)

for ( i in 1: irrel.D){
  
  matTemp <- matrix(runif(n = (rk -rel.D), min = -1.5, max= 1.5),  nrow = (rk-rel.D), ncol =1)
  irrel.X[,i] <- as.vector(linear.space %*% matTemp)
  
}

## Checking if the covariance is indeed small

cov.mat <- cov(rel.X,irrel.X)

boxplot(cov.mat)

## Building the full data matrix

X.full <- cbind(rel.X, irrel.X)


levelplot(cov(X.full))

Y <- X.full

#########################################################################################
##########################################################################################
##### Now WE DEAL WITH CLUSTERED DATA AND GENERATE NON RELEVANT FEATURES INDEPENENTLY ###


## Selcting the beta co-efficients

## The Co-efficients have to be obtained from uniform distribution between [-3,3]
beta.list <- list(0)
for ( i in 1:F){
  beta.list[[i]] <- runif(rel.D,-3,3)
}

## Let's See if the Clusters are separate WITH AND WITHOUT THE RELEVANT FEATURES
Y.rel <- c(0)
for (i in 1:F){
  Y.rel <- rbind(Y.rel, Y.rel.list[[i]])
} 
Y.rel <- Y.rel[-1,]


## True Labels for the points
c.true <<- c(0)
for ( i in 1:F){
  c.true <- rbind(as.matrix(c.true) , as.matrix(c(rep(i, as.integer(N * p.dist[i])))))  
}
c.true <- as.factor(c.true[-1,])






## All features
pc <- prcomp(Y)
pc.pred <- predict(pc,newdata = Y)
plot(pc.pred[,1], pc.pred[,2], pch = 19, col =c.true, main = 'Main Data')

## Only with the relevant features
pc <- prcomp(Y.rel)
pc.pred <- predict(pc,newdata = Y.rel)
plot(pc.pred[,1], pc.pred[,2], pch = 19, col =c.true, main = 'Time Relevant Data')


## The pure time is generated
time.pur.list <- list(0)
for ( i in 1:F){
  time.pur.list[[i]] <- t(beta.list[[i]]) %*% t(Y.rel.sc.list[[i]])
}


time.noise.list <- list(0)
for ( i in 1:F){
  time.noise.list[[i]] <- rnorm(length(which(c.true == i)), mean = 3*i, sd = 0.1*i )
}

time.list <- list(0)
for ( i in 1:F){
  time.list[[i]] <- time.pur.list[[i]] + time.noise.list[[i]]
}

#######################################MAKING TIME from cluster data ########################################################
## Real time without Censoring
time.real <- c(0)
for (i in 1:F){
  time.real <- cbind(time.real, time.list[[i]])
} 
time.real <- time.real[,-1]
time.real <- as.vector(unlist(time.real))



####################################### Adding CENSORING INFORMATION  ################################################################
## Adding the censoring information to the TIME

# censoring <- rbinom(n = NROW(Y), size =1, prob = 1- prob.censoring)
# right.censoring.time <- min(time.real)  
# 
# time <- time.real
# 
# index.time <- which(censoring==0)
# for ( q in 1:length(index.time)){
#   time[index.time[q]] <- right.censoring.time
#   
# }

time <- c(0)
censoring <- c(rep(1,N))
time <- time.real

## Boxplots for Vizualization of the time Data without censoring
boxplot(time.list)


### A Quick ManWhittney U / Kruskal test  test to check if the time's of the two cluster are significantly different
## wilcox.test(as.vector(time.list[[1]]), as.vector(time.list[[2]]), alternative = "two.sided")

kruskal.test(time, c.true)


### Visualization of the different survival curves
surv.ob <- Surv(time,censoring)
survfit <- survfit(surv.ob ~ c.true)
plot(survfit,)





assign("rel.D", rel.D, envir = .GlobalEnv)
assign("irrel.D", irrel.D, envir = .GlobalEnv)
assign("Y", Y, envir = .GlobalEnv)
assign("beta.list", beta.list, envir = .GlobalEnv)
assign("c.true", c.true, envir = .GlobalEnv)
assign("c.true", c.true, envir = .GlobalEnv)
assign("time.real", time.real, envir = .GlobalEnv)
assign("time", time, envir = .GlobalEnv)
assign("censoring", censoring, envir = .GlobalEnv)




}

