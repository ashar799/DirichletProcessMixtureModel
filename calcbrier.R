
nmrse <- c(0)
betahat.brier.list <- list(0) 
iter = 100
iter.thin =1
print("GIBB'S SAMPLING")
count = 1


for (o in 1:iter) {
  ################## PARAMETERS OF THE DP Mixture Model ######################################################
  ## Updating the parameters based on the observations 
  param <- posteriorGMMparametrs(c,Y,mu,S, alpha,K, epsilon, W, beta, ro,N,D )
  mu <- param$mean
  S <- param$precision
  paramtime <- posteriortimeparameters(c, That, lambda2,tau2,sigma2,beta0, betahat, Y, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.data)
  beta0 <- paramtime$beta0
  betahat <- paramtime$betahat
  sigma2 <- paramtime$sigma2
  lambda2 <- paramtime$lambda2
  tau2 <- paramtime$tau2
  
  ########################## THE HYPERPARAMETERS OF THE GMM #################################  
  source('posteriorhyper.R')  
  #  Updating the hyper paramters
  hypercognate <- posteriorhyper (c, Y, mu, S, epsilon, W, beta, ro )
  epsilon <- hypercognate$epsilon
  tmpW <- hypercognate$W
  W <- matrix(as.matrix(tmpW),nrow = D, ncol =D)
  ro <- hypercognate$ro
  
  ################# INDICATOR VARIABLE ##################################################################
  ## Updating the indicator variables and the parameters
  source('posteriorchineseAFT.R')
  cognate <- posteriorchineseAFT(c,Y,mu,S,alpha,That, beta0, betahat, sigma2, lambda2, tau2, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.dat)
  c <- cognate$indicator
  mu <- cognate$mean
  S <- cognate$precision
  beta0 <- cognate$beta0
  betahat <- cognate$betahat
  sigma2 <- cognate$sigma2
  lambda2 <- cognate$lambda2
  tau2 <- cognate$tau2
  
  ########################### The Concentration Parameter #################################################################
  
  
  source('posterioralpha.R') 
  # Updating the concentration parameter
  alpha <- posterioralpha(c, N, alpha, shape.alpha, rate.alpha)
  
  
  
  
  
  ######################## The Censored Times ###########################################################
  source('updatetime.R')
  # Updating the Time Variable
  ti <- NA
  ti <- updatetime(c, Y, Time,That, beta0, betahat, sigma2)
  That <- ti$time
  
  
  ################## PARAMETERS OF THE DP Mixture Model ######################################################
  ## Updating the parameters based on the observations 
  
  if(o%% iter.thin == 0 ){
   betahat.brier.list[[count]] <- betahat  
    count <- count +1
  }
  
  
  
  
  print(o/iter) 
  #   print(loglike[o])
  #   print(cindex)
} 

count =  count -1

list.brier.betahat <- list(0)

for ( i in 1:count){
  list.brier.betahat[[i]] <- (betahat.brier.list[[i]][1:2,] != 0) +0
}


betahat1.brier.final <- matrix(NA, nrow = count, ncol = D)
for ( i in 1:count){
  betahat1.brier.final[i,] <- list.brier.betahat[[i]][1,]
}
betahat2.brier.final <- matrix(NA, nrow = count, ncol = D)
for ( i in 1:count){
  betahat2.brier.final[i,] <- list.brier.betahat[[i]][2,]
}


relev <- rep(1, 10)
irrel <- rep(0, 10)

final.brier1 <- matrix(0, nrow = 10, ncol =D)
final.brier2 <- matrix(0, nrow = 10, ncol =D)
ten.briers <- list(0)
 
for ( i in 1:10){
  start = 10*(i-1) + 1
  end = 10*i 
  seq <- seq(start,end)
  final.brier1[i,] <- as.vector(apply(betahat1.brier.final[seq,],2,mean))
  final.brier2[i,] <- as.vector(apply(betahat2.brier.final[seq,],2,mean))
   }

score.brier.rel.1 <- matrix(NA, nrow = 1, ncol =rel.D)
score.brier.irrel.1 <- matrix(NA, nrow = 1, ncol =irrel.D)
score.brier.rel.2 <- matrix(NA, nrow = 1, ncol =rel.D)
score.brier.irrel.2 <- matrix(NA, nrow = 1, ncol =irrel.D)

for ( i in 1: rel.D){
   score.brier.rel.1[1,i] <- brier(obs = relev, pred = final.brier1[,i] )$bs
 }



for ( i in 1: irrel.D){
  temp <- i + rel.D
  score.brier.irrel.1[1,i] <- brier(obs = irrel,pred = final.brier1[,temp] )$bs
}

for ( i in 1: rel.D){
  score.brier.rel.2[1,i] <- brier(obs = relev,pred = final.brier2[,i] )$bs
}
for ( i in 1: irrel.D){
  temp <- i + rel.D
  score.brier.irrel.2[1,i] <- brier(obs = irrel, pred = final.brier2[,temp] )$bs
}

score.brier.final1 <-  cbind(score.brier.rel.1,score.brier.irrel.1)
score.brier.final2 <-  cbind(score.brier.rel.2,score.brier.irrel.2)


###### Plot the brier scores #####################################3
brier.join <- rbind(score.brier.final1, score.brier.final2)
rownames(brier.join) = c("cluster_1","cluster_2")
colnames(brier.join) =  c(rep("relevant",rel.D),rep("irrelevant",irrel.D))


heatmapdata <- as.data.frame(brier.join)
pdf("/home/bit/ashar/ExpressionSets/Simulations/Brier.pdf")
heatmap.2(t(as.matrix(heatmapdata)),dendrogram="none", col =cm.colors(10), margins=c(6,10), main = "Brier Scores \n over \n 10 repetitions ", cexCol = 0.85, cexRow = 0.7, trace = "none",Rowv =FALSE)
dev.off()





### Final bethats
final.betahat1 <- apply(betahat1.final,2,mean)
final.betahat2 <- apply(betahat2.final,2,mean)


### Probability of betahat of genes
final.betahat <- rbind(final.betahat1, final.betahat2)
rownames(final.betahat) = c("cluster_1","cluster_2")

