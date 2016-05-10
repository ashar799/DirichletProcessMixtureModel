############ This file takes in training and testing data for Glioblastoma #################
#################################################################################################
rm(list =ls())
################################################
######## Load the Data from the folder #########
load("/home/bit/ashar/ownCloud/DPMM_RESULTS/ONE_VIEW/Verhark/Final/DataVerhaak.RData")
Y.train.prelim <- relev$Y.train.prelim
Y.test.prelim  <- relev$Y.test.prelim
signature.vk <- relev$signature
signature.dpmm <- relev$signature.dpmm
pheno.train <- relev$pheno.train
pheno.test <- relev$pheno.test


Y <- t(Y.train.prelim[as.character(signature.dpmm),])
time <- pheno.train[,3]
censoring <- pheno.train[,2]
c.true <- pheno.train[,4]
levels(c.true) <- c(1,2,3,4)

Y.new <- t(Y.test.prelim[as.character(signature.dpmm),])
time.new <- pheno.test[,3]
censoring.new <- pheno.test[,2]
c.true.new <- pheno.test[,4]
levels(c.true.new) <- c(1,2,3,4)


##### checking Heatmaps ##################
hmcols<-colorRampPalette(c("green","black","red"))(512)
heatmap.2(x = t(Y), scale = "row",Rowv =TRUE ,Colv = FALSE, dendrogram = 'row', ColSideColors =c("blue","red","green","purple")[as.numeric(c.true)], labRow = colnames(Y), labCol = NA, main = ' \n Training Set \n ', col =hmcols, cexRow =0.40, trace ="none", key.title = "Color Key")






N <- nrow(Y)
D <- ncol(Y)
N.new <- nrow(Y.new)

############################# PARAMETERS for GIBB's SAMPLING ####
iter = 40
iter.burnin = 100
iter.thin  = 2
k = 4
Nps = 20

######################### Initialize the Parameters ################
source('initializeDPMM.R')
initializeDPMM()


########### Train the Model #########################################
source('burninDPMM.R')
burninDPMM()

source('gibbsDPMM.R')
gibbsDPMM()

########## Analyze the fit ##########################################
### Good feature selection from heatmap plus cindex plus randindex
source('SIManalyzeDPMM.R')
SIManalyzeDPMM()


############ Generating some Plots ##########################
pc <- prcomp(Y)
pc.pred <- predict(pc,newdata = Y)
p1 <- ggplot(as.data.frame(pc.pred), aes(x=pc.pred[,1], y= pc.pred[,2], colour= as.factor(c.final))) + ggtitle(" DPMM Clustering") + geom_point(shape=19) + labs(y = "PC1", x = "PC2", colour = "Classes") 

##### checking Heatmaps ##################
hmcols<-colorRampPalette(c("green","black","red"))(512)
heatmap.2(x = t(Y), scale = "row",Rowv =TRUE ,Colv = TRUE , dendrogram = 'row', ColSideColors =c("blue","red","green","purple")[as.numeric(c.final)], labRow = colnames(Y), labCol = NA, main = ' \n Training Set \n ', col =hmcols, cexRow =0.40, trace ="none", key.title = "Color Key")

####Fitting Survival Curves
surv.ob <- Surv(exp(time),censoring)
surv.fit <- survfit(surv.ob ~ c.final)
logrank <- survdiff(surv.ob ~ c.final)
p5 <- ggsurv(surv.fit, main = " DPMM \n Kaplan Meier Estimators") + ggplot2::guides(linetype = FALSE) + ggplot2::scale_colour_discrete(name = 'Classes',breaks = c(1,2),labels = c('1', '2'))
adjustedRandIndex(c.true, c.final)


#### Getting recovered C-Index ######
source('predictchineseAFTtime.R')
predictchineseAFTtime(Y)
### Check of the Predicted C-index 
recovered.cindex <- survConcordance(Surv(exp(time),censoring) ~ exp(-post.time.avg))[1]


################### Predicting C-Index ############3
source('predictchineseAFTtime.R')
predictchineseAFTtime(Y.new)
### Check of the Predicted C-index 
predicted.cindex <- survConcordance(Surv(exp(time.new),censoring.new) ~ exp(-post.time.avg))[1]

############################################################
surv.ob.new <- Surv(exp(time.new),censoring.new)
######## Predict on New Data Set  BASED ON JUST THE MOLECULAR DATA #####################################
source('predictCLUSTER.R')
predictCLUSTER(Y.new)
## Check the predicted Rand Index 
lr <- c(0)
for (j in 1:Nps){
  lr[j] <-  1 - pchisq(unlist(survdiff(surv.ob.new ~ c.new.list[[j]]))$chisq,df = length(table(c.new.list[[j]])) - 1 )
}
c.final.new <- c.new.list[[which.min(lr)]]



##### Generating some plots #####################################
pc <- prcomp(Y.new)
pc.pred <- predict(pc,newdata = Y.new)
p1 <- ggplot(as.data.frame(pc.pred), aes(x=pc.pred[,1], y= pc.pred[,2], colour= as.factor(c.final.new))) + ggtitle(" DPMM Clustering") + geom_point(shape=19) + labs(y = "PC1", x = "PC2", colour = "Classes") 
surv.fit <- survfit(surv.ob.new ~ c.final.new)
logrank <- survdiff(surv.ob.new ~ c.final.new)
p5 <- ggsurv(surv.fit, main = " DPMM \n Kaplan Meier Estimators") + ggplot2::guides(linetype = FALSE) + ggplot2::scale_colour_discrete(name = 'Classes',breaks = c(1,2,3,4),labels = c('1', '2','3',4))


############################################################################################################################################################
############################################################################################################################################################
######## ANLYSIS USING DPMM signature with K-Means ##################################################################
########### ANALYSIS USING DPMM SIGNATURE WITH K-NN + PCOX ##############################################################
#############################################################################################################################################################
#############################################################################################################################################################
label.train <- kmeans(Y, centers =4, nstart =10)$cluster
label.test <- knn(train = Y, test = Y.new, cl = label.train, k = 4)
adjustedRandIndex(label.train,c.true)
adjustedRandIndex(label.test, c.true.new)

surv.ob <- Surv(exp(time),censoring)
surv.ob.new <- Surv(exp(time.new),censoring.new)

survdiff(surv.ob ~ label.train)
survdiff(surv.ob.new ~ label.test)

linear.kkpcox.recovery <- c(0)
linear.kkpcox.prediction <- c(0)
for ( q in 1:F){
  ind <- which(label.train == q)
  ind.new <- which(label.test == q)
  reg.cox <- cv.glmnet(x = Y[ind,], y = Surv(exp(time[ind]),censoring[ind]), family = "cox", maxit = 10000000)
  linear.kkpcox.recovery[ind] <- predict(object =reg.cox, newx = Y[ind,])
  linear.kkpcox.prediction[ind.new] <- predict(object =reg.cox, newx = Y.new[ind.new,])
}
as.numeric(survConcordance(Surv(exp(time), censoring) ~ linear.kkpcox.recovery)[1])
as.numeric(survConcordance(Surv(exp(time.new), censoring.new) ~ linear.kkpcox.prediction)[1])








############################################################################################################################################################
############################################################################################################################################################
######## ANLYSIS USING VERHAAK SIGNATURE ##################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################

#### LogRank ####
logrank.vv <- survdiff(surv.ob ~ c.true)


####### Recovering C-Indexes ##############
### Verhaak Signature
linear.vv.prediction <- c(0)
Y.vk.train <- t(Y.train.prelim[signature.vk,])
Y.vk.test <- t(Y.test.prelim[signature.vk,])


for ( q in 1:F){
  ind <- which(c.true == q)
  reg.cox <- cv.glmnet(x = Y.vk.train[ind,], y = Surv(exp(time[ind]),censoring[ind]), family = "cox", maxit = 10000000)
  linear.vv.prediction[ind] <- predict(object =reg.cox, newx = Y.vk.train[ind,], s= "lambda.1se")
}
predCIndex.vv.pcox <<- as.numeric(survConcordance(Surv(exp(time), censoring) ~ linear.vv.prediction)[1])


####### COMPARISON METHOD Their Labels for training + KNN (Penalized Cox) ###############################
label.train <- c.true
label.test <- knn(train = Y.vk.train, test = Y.vk.test, cl = label.train, k = F)
adjustedRandIndex(label.test, c.true.new)

surv.ob <- Surv(exp(time.new),censoring.new)
surv.fit <- survfit(surv.ob ~ label.test)
logrank <- survdiff(surv.ob ~ label.test)

linear.kkpcox.prediction <- c(0)
for ( q in 1:F){
  ind <- which(label.train == q)
  ind.new <- which(label.test == q)
  reg.cox <- cv.glmnet(x = Y.vk.train[ind,], y = Surv(exp(time[ind]),censoring[ind]), family = "cox", maxit = 10000000)
  linear.kkpcox.prediction[ind.new] <- predict(object =reg.cox, newx = Y.vk.test[ind.new,])
}
predCIndex.kkpcox <<- as.numeric(survConcordance(Surv(exp(time.new), censoring.new) ~ linear.kkpcox.prediction)[1])


####################################################################################################################################################
####################################################################################################################################################
###################### VISUALIZATION ###############################################################################################################
####################################################################################################################################################

##########################################################
########## Creating a Ranking for the Points ####################
###### Creating a Ratio for the Points ########################
rank <- matrix(0, nrow = N, ncol =Nps)

for (j in 1:Nps){
  
  Y.scaled <- matrix(0, nrow = N, ncol =D)
  for ( v in 1:4){
    clust <- which(c.list[[j]] == v)
    Y.scaled[clust,1:D] <- scale(Y[clust,1:D], center = TRUE, scale = TRUE)
  }
  
  
  for ( i in 1:N){
    rank[i,j] <- dMVN(as.vector(t(Y[i,1:D])), mean = mu.list[[j]][c.final[i],1:D], Q= S.list[[j]][c.final[i],1:D,1:D], log = TRUE) - dMVN(as.vector(t(Y[i,1:D])), mean = mu.list[[j]][1,1:D], Q= S.list[[j]][1,1:D,1:D], log = TRUE) +  dnorm(x = That[i], mean = beta0.list[[j]][c.final[i]] + betahat.list[[j]][c.final[i],1:D] %*% as.vector(t(Y.scaled[i,1:D])), sd = sqrt(sigma2.list[[j]][c.final[i]]), log =TRUE) -  dnorm(x = That[i], mean = beta0.list[[j]][1] + betahat.list[[j]][1,1:D] %*% as.vector(t(Y.scaled[i,1:D])), sd = sqrt(sigma2.list[[j]][1]), log =TRUE) 
  }
}

avg.rank <- apply(rank,1,mean)
order.zo <- range01(avg.rank)

order.train <- sort(order.zo,index.return = TRUE, decreasing = TRUE)
Y.order <- Y[order.train$ix,]
c.final.order <- c.final[order.train$ix]


######## Reordering Again ################
order.2 <- c(which(c.final.order==1),which(c.final.order==3),which(c.final.order==4),which(c.final.order==2))
Y.order.2 <- Y.order[order.2,]
c.final.order.2 <- c.final.order[order.2]




surv.ob <- Surv(exp(time),censoring)
Classes <- c.final
Classes <- as.factor(Classes)
surv.fit <- survfit(surv.ob ~ Classes)
p5 <- ggsurv(surv.fit, main = " Kaplan Meier Estimators \n pvalue 5e -05", surv.col = c("green","red","blue","orange"), cens.col ="Blue")  


############ Generating some Plots ##########################
pc <- prcomp(Y)
pc.pred <- predict(pc,newdata = Y)
p1 <- ggplot(as.data.frame(pc.pred), aes(x=pc.pred[,1], y= pc.pred[,2], color = as.factor(c.final) )) + ggtitle(" DPMM Clustering \n 47 Gene Signature Verhaak \n 42% Variance Explained") + geom_point(shape=19) + labs(y = "PC1", x = "PC2") + scale_color_manual(breaks = c("1", "2", "3","4"),values=  c("green","red","blue","orange"))



pdf('DPMMSignature.pdf')
hmcols<-colorRampPalette(c("green","black","red"))(256)
heatmap.2(x = t(Y.order.2), Rowv =TRUE ,Colv = FALSE, dendrogram = c("row"), scale ="row",ColSideColors =c("green","red","blue","orange")[c.final.order.2], labRow = colnames(Y.order.2), labCol = NA, main = ' \n Training Set \n 47 Gene DPMM signature', col =hmcols, cexRow =0.40, trace ="none", key.title = "Color Key")
legend("topright", legend = c("Best","Good Moderate","Bad moderate","Worst"),fill = c("Green","Blue","Orange","Red"), cex = 0.4)
p5
p1
dev.off()




############## Similar Plots for Testing Data Set ######################
########################################################################


rank <- matrix(0, nrow = N, ncol =Nps)

for (j in 1:Nps){
  
  Y.scaled <- matrix(0, nrow = N, ncol =D)
  for ( v in 1:4){
    clust <- which(c.new.list.save[[j]] == v)
    Y.scaled[clust,1:D] <- scale(Y.new[clust,1:D], center = TRUE, scale = TRUE)
  }
  
  
  for ( i in 1:N){
    rank[i,j] <- dMVN(as.vector(t(Y.new[i,1:D])), mean = mu.list[[j]][c.final.new[i],1:D], Q= S.list[[j]][c.final.new[i],1:D,1:D], log = TRUE) - dMVN(as.vector(t(Y.new[i,1:D])), mean = mu.list[[j]][1,1:D], Q= S.list[[j]][1,1:D,1:D], log = TRUE) +  dnorm(x = time.new[i], mean = beta0.list[[j]][c.final.new[i]] + betahat.list[[j]][c.final.new[i],1:D] %*% as.vector(t(Y.scaled[i,1:D])), sd = sqrt(sigma2.list[[j]][c.final.new[i]]), log =TRUE) -  dnorm(x = time.new[i], mean = beta0.list[[j]][1] + betahat.list[[j]][1,1:D] %*% as.vector(t(Y.scaled[i,1:D])), sd = sqrt(sigma2.list[[j]][1]), log =TRUE) 
  }
}

avg.rank <- apply(rank,1,mean)
order.zo <- range01(avg.rank)

order.train <- sort(order.zo,index.return = TRUE, decreasing = TRUE)
Y.order.new <- Y.new[order.train$ix,]
c.final.order.new <- c.final.new[order.train$ix]


######## Reordering Again ################
order.2.new <- c(which(c.final.order.new==1),which(c.final.order.new==3),which(c.final.order.new==4),which(c.final.order.new==2))
Y.order.2.new <- Y.order.new[order.2.new,]
c.final.order.2.new <- c.final.order.new[order.2.new]


surv.ob <- Surv(exp(time.new),censoring.new)
Classes <- c.final.new
Classes <- as.factor(Classes)
surv.fit <- survfit(surv.ob ~ Classes)
p5 <- ggsurv(surv.fit, main = " DPMM \n Kaplan Meier Estimators \n pvalue 3e -02", surv.col = c("green","red","blue","orange"), cens.col ="Blue")  


############ Generating some Plots ##########################
pc <- prcomp(Y.new)
pc.pred <- predict(pc,newdata = Y.new)
p1 <- ggplot(as.data.frame(pc.pred), aes(x=pc.pred[,1], y= pc.pred[,2], color = as.factor(c.final.new) )) + ggtitle(" DPMM Clustering\n Testing Data \n 47 Gene Signature Verhaak \n 44% Variance Explained") + geom_point(shape=19) + labs(y = "PC1", x = "PC2") + scale_color_manual(breaks = c("1", "2", "3","4"),values=  c("green","red","blue","orange"))



pdf('TestingDPMMSignature.pdf')
hmcols<-colorRampPalette(c("green","black","red"))(256)
heatmap.2(x = t(Y.order.2.new), Rowv =TRUE ,Colv = FALSE, dendrogram = c("row"), scale ="row",ColSideColors =c("green","red","blue","orange")[c.final.order.2.new], labRow = colnames(Y.order.2.new), labCol = NA, main = ' \n Testing Set \n 47 Gene DPMM signature', col =hmcols, cexRow =0.40, trace ="none", key.title = "Color Key")
legend("topright", legend = c("Best","Good Moderate","Bad moderate","Worst"),fill = c("Green","Blue","Orange","Red"), cex = 0.4)
p5
p1
dev.off()

#########################################################################################################################################################
########################################################### Feature Selection #############################################################################
###########################################################################################################################################################

pdf("Heatmap_FeatureSelection_VerhaakDataSet.pdf")
hmcols<-colorRampPalette(c("white","black"))(128)
heatmap.2(t(as.matrix(heatmapdata[c(3,4),])),dendrogram="none", col =hmcols, margins=c(6,10), main = "Posterior prob.for Selection \n Verhaak Data Set ", cexCol = 0.85, cexRow = 0.7, Rowv = FALSE, labCol = c("Best","Worst","Moderate_Good","Moderate_Bad"))
dev.off()



####################################################################################################################################################
####################################################################################################################################################
########################### MAKING FUNCTIONAL /BIOLOGICAL SENSE OF THE DPMM SIGNATURE ##############################################################
####################################################################################################################################################
####################################################################################################################################################
