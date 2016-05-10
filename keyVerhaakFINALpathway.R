############ This file takes in training and testing data for Glioblastoma #################
#################################################################################################
rm(list =ls())
################################################
######## Load the Data from the folder #########
load("/home/bit/ashar/ownCloud/DPMM_RESULTS/ONE_VIEW/Verhark/30pathwayVerhaark.RData")
load("/home/bit/ashar/ownCloud/DPMM_RESULTS/ONE_VIEW/Verhark/phenoVerhaark.RData")


############################################################################################################################################################
##### Using the 840 Gene Signature For training and Testing Set and Filtering out Genes which have low survival information #################################
#############################################################################################################################################################
N.total <- nrow(data.combined)
index.train <- seq(1, N.total, by =2)
index.test <- seq(2, N.total, by =2)

Y.train <- data.combined[index.train,]
Y.test <- data.combined[index.test,]
pheno.train <- pheno[index.train,]
pheno.test <- pheno[index.test,]




Y <- Y.train
time <- pheno.train[,3]
censoring <- pheno.train[,2]
c.true <- pheno.train[,4]
levels(c.true) <- c(1,2,3,4)

Y.new <- Y.test
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
heatmap.2(x = t(Y), scale = "row",Rowv =TRUE ,Colv = TRUE, dendrogram = 'row', ColSideColors =c("blue","red","green","purple")[as.numeric(c.final)], labRow = colnames(Y), labCol = NA, main = ' \n Training Set \n ', col =hmcols, cexRow =0.40, trace ="none", key.title = "Color Key")

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
  lr[j] <-  1 - pchisq(unlist(survdiff(surv.ob.new ~ c.new.list.save[[j]]))$chisq,df =3)
}
c.final.new <- c.new.list.save[[which.min(lr)]]



##### Generating some plots #####################################
pc <- prcomp(Y.new)
pc.pred <- predict(pc,newdata = Y.new)
p1 <- ggplot(as.data.frame(pc.pred), aes(x=pc.pred[,1], y= pc.pred[,2], colour= as.factor(c.final.new))) + ggtitle(" DPMM Clustering") + geom_point(shape=19) + labs(y = "PC1", x = "PC2", colour = "Classes") 
surv.fit <- survfit(surv.ob.new ~ c.final.new)
logrank <- survdiff(surv.ob.new ~ c.final.new)
p5 <- ggsurv(surv.fit, main = " DPMM \n Kaplan Meier Estimators") + ggplot2::guides(linetype = FALSE) + ggplot2::scale_colour_discrete(name = 'Classes',breaks = c(1,2),labels = c('1', '2'))


############################################################################################################################################################
############################################################################################################################################################
######## ANLYSIS USING DPMM signature with K-Means ##################################################################
########### ANALYSIS USING DPMM SIGNATURE WITH K-NN + PCOX ##############################################################
#############################################################################################################################################################
#############################################################################################################################################################
label.train <- kmeans(Y, centers =2, nstart =10)$cluster
label.test <- knn(train = Y, test = Y.new, cl = label.train, k = 2)
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
######## ANLYSIS USING VAN'T VEER DATA SET ##################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################

#### LogRank ####
logrank.vv <- survdiff(surv.ob ~ c.true)


####### Recovering C-Indexes ##############
### Van't Veer Signature
linear.vv.prediction <- c(0)
for ( q in 1:F){
  ind <- which(c.true == q)
  reg.cox <- cv.glmnet(x = Y.vijver.train[ind,], y = Surv(exp(time[ind]),censoring[ind]), family = "cox", maxit = 10000000)
  linear.vv.prediction[ind] <- predict(object =reg.cox, newx = Y.vijver.train[ind,], s= "lambda.min")
}
predCIndex.vv.pcox <<- as.numeric(survConcordance(Surv(exp(time), censoring) ~ linear.vv.prediction)[1])


####### COMPARISON METHOD Their Labels for training + KNN (Penalized Cox) ###############################
label.train <- c.true
label.test <- knn(train = Y, test = Y.new, cl = label.train, k = 2)
adjustedRandIndex(label.test, c.true.new)

surv.ob <- Surv(exp(time.new),censoring.new)
surv.fit <- survfit(surv.ob ~ label.test)
logrank <- survdiff(surv.ob ~ label.test)

linear.kkpcox.prediction <- c(0)
for ( q in 1:F){
  ind <- which(label.train == q)
  ind.new <- which(label.test == q)
  reg.cox <- cv.glmnet(x = Y.vijver.train[ind,], y = Surv(exp(time[ind]),censoring[ind]), family = "cox", maxit = 10000000)
  linear.kkpcox.prediction[ind.new] <- predict(object =reg.cox, newx = Y.vijver.test[ind.new,])
}
predCIndex.kkpcox <<- as.numeric(survConcordance(Surv(exp(time.new), censoring.new) ~ linear.kkpcox.prediction)[1])



###### Using HC on Van't Veer data set ############
d <- dist(as.matrix(Y.vijver.train)) 
hc <- hclust(d) 
memb <- cutree(hc, k = F)
survdiff(Surv(exp(time),censoring) ~ memb)
adjustedRandIndex(c.true, memb)

### HClust Signature
linear.hc.prediction <- c(0)
for ( q in 1:F){
  ind <- which(memb == q)
  reg.cox2 <- cv.glmnet(x = Y.vijver.train[ind,], y = Surv(exp(time[ind]),censoring[ind]), family = "cox", maxit =100000000)
  linear.hc.prediction[ind] <- predict(object =reg.cox2, newx = Y.vijver.train[ind,])
}
predCIndex.hc.pcox <<- as.numeric(survConcordance(Surv(exp(time), censoring) ~ linear.hc.prediction)[1])


label.train <- as.factor(memb)
label.vijver.test <-  knn(train = Y.vijver.train, test = Y.vijver.test, cl = label.train, k = 5)
surv.ob.new <- Surv(exp(time.new),censoring.new)
surv.fit <- survfit(surv.ob.new ~ label.vijver.test)
logrank <- survdiff(surv.ob.new ~ label.vijver.test)

linear.kkpcox.prediction <- c(0)
for ( q in 1:F){
  ind <- which(label.train == q)
  ind.new <- which(label.vijver.test == q)
  reg.cox <- cv.glmnet(x = Y.vijver.train[ind,], y = Surv(exp(time[ind]),censoring[ind]), family = "cox", maxit =1000000000)
  linear.kkpcox.prediction[ind.new] <- predict(object =reg.cox, newx = Y.vijver.test[ind.new,], s = reg.cox$lambda[10])
}
predCIndex.kkpcox <<- as.numeric(survConcordance(Surv(exp(time.new), censoring.new) ~ linear.kkpcox.prediction)[1])



################# Calculating Hazard Ratio #################################

##### For the Classification Given by Van't Veer ########################

#### On Testing Data Set ####################

######### Hazard Ratio for the Van't Veer signature ######################################################3
hazard.ratio(x = c.true.new ,surv.time= exp(time.new), surv.event= censoring.new)

hazard.ratio(x = c.new ,surv.time= exp(time.new), surv.event= censoring.new)

hazard.ratio(x = label.vijver.test ,surv.time= exp(time.new), surv.event= censoring.new)


##### Visualization of the Results ###############################
#### USE Top 35 Molecular Features OBTAINED FROM THE MODEL ############
sum.diff <- c(0)
for ( i in 1:Nps){
  sum.diff <- sum.diff + abs(mu.list[[i]][1,] - mu.list[[i]][2,])  
}
sum.diff.sc <- (1/Nps) * (sum.diff)/(diag(W))

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
cluste <- range01(sum.diff.sc)

D.good <- (sum.diff.sc >  11)
Y.good <- Y[,D.good]
Y.new.good <- Y.new[,D.good]


##########################################################
########## Creating a Ranking for the Points ####################
###### Creating a Ratio for the Points ########################
rank <- matrix(0, nrow = N, ncol =Nps)

for (j in 1:Nps){
  
  Y.scaled <- matrix(0, nrow = N, ncol =D)
  for ( v in 1:2){
    clust <- which(c.list[[j]] == v)
    Y.scaled[clust,1:D] <- scale(Y[clust,1:D], center = TRUE, scale = TRUE)
  }
  
  
  for ( i in 1:N){
    rank[i,j] <- dMVN(as.vector(t(Y[i,1:D])), mean = mu.list[[j]][1,1:D], Q= S.list[[j]][1,1:D,1:D], log = TRUE) - dMVN(as.vector(t(Y[i,1:D])), mean = mu.list[[j]][2,1:D], Q= S.list[[j]][2,1:D,1:D], log = TRUE) +  dnorm(x = That[i], mean = beta0.list[[j]][1] + betahat.list[[j]][1,1:D] %*% as.vector(t(Y.scaled[i,1:D])), sd = sqrt(sigma2.list[[j]][1]), log =TRUE) -  dnorm(x = That[i], mean = beta0.list[[j]][2] + betahat.list[[j]][2,1:D] %*% as.vector(t(Y.scaled[i,1:D])), sd = sqrt(sigma2.list[[j]][2]), log =TRUE) 
  }
}

avg.rank <- apply(rank,1,mean)
order.zo <- range01(avg.rank)

order.train <- sort(order.zo,index.return = TRUE, decreasing = TRUE)
Y.order <- Y[order.train$ix,]
c.final.order <- c.final[order.train$ix]


surv.ob <- Surv(exp(time),censoring)
Classes <- c.final
Classes <- as.factor(Classes)
surv.fit <- survfit(surv.ob ~ Classes)
p5 <- ggsurv(surv.fit, main = " DPMM \n Kaplan Meier Estimators", surv.col = c("Red","Black"), cens.col ="Blue")  


pdf('DPMMSignature.pdf')
hmcols<-colorRampPalette(c("green","black","red"))(256)
heatmap.2(x = t(Y.order), Rowv =TRUE ,Colv = FALSE, dendrogram = c("row"), scale ="row",ColSideColors =c("red","black")[c.final.order], labRow = colnames(Y.vijver.train), labCol = NA, main = ' \n Training Set \n 58 Gene DPMM signature', col =hmcols, cexRow =0.40, trace ="none", key.title = "Color Key")
legend("topright", legend = c("Bad Prognosis","Good Prognosis"),fill = c("Red","Black"), cex = 0.4)
p5
dev.off()



############################ Let's See if we can git similar plots on the Testing Data #############################################
#### We just take That Model out of the MCMC samples which has the maximum likelihhod

rank.test <- c(0)

index.max <- which.max(model.weight)   

for ( i in 1:N.new){
  rank.test[i] <- dMVN(as.vector(t(Y.new[i,1:D])), mean = mu.list[[index.max]][1,1:D], Q= S.list[[index.max]][1,1:D,1:D], log = TRUE) - dMVN(as.vector(t(Y.new[i,1:D])), mean = mu.list[[index.max]][2,1:D], Q= S.list[[index.max]][2,1:D,1:D], log = TRUE) 
}



order.zo <- range01(rank.test)

order.train <- sort(order.zo,index.return = TRUE, decreasing = TRUE)
Y.order.new <- Y.new[order.train$ix,]
c.final.order.new <- c.new.list.save[[index.max]][order.train$ix]


surv.ob <- Surv(exp(time.new),censoring.new)
Classes <- c.new.list.save[[1]]
Classes <- as.factor(Classes)
levels(Classes) <- c("Bad","Good")
surv.fit <- survfit(surv.ob ~ Classes)
p5 <- ggsurv(surv.fit, main = " DPMM \n Kaplan Meier Estimators", surv.col = c("Red","Black"), cens.col ="Blue")  


pdf('TEST-Data.pdf')
hmcols<-colorRampPalette(c("green","black","red"))(256)
heatmap.2(x = t(Y.order.new), Rowv =TRUE ,Colv = FALSE, dendrogram = c("row"), scale ="row",ColSideColors =c("red","black")[c.final.order.new], labRow = colnames(Y.vijver.test), labCol = NA, main = ' \n Testing Set \n 58 Gene DPMM signature', col =hmcols, cexRow =0.40, trace ="none", key.title = "Color Key")
legend("topright", legend = c("Bad Prognosis","Good Prognosis"),fill = c("Red","Black"), cex = 0.4)
p5
dev.off()





### Make Heatmap Plots
hmcols<-colorRampPalette(c("green","black","red"))(256)
heatmap.2(x = t(Y.good), Rowv =TRUE ,Colv = TRUE, dendrogram = c("row"), scale ="row",ColSideColors =c("red","black")[as.factor(c.final)], labRow = colnames(Y.good), labCol = NA, main = 'VIJVER (2002) \n Training Set \n 12 Genes High Rank Genes', col =hmcols, cexRow =0.40, trace ="none", key.title = "Color Key")
legend("topright", legend = c("Good Prognosis","Bad Prognosis"),fill = c("Red","Black"), cex = 0.4)

heatmap.2(x = t(Y.new.good), Rowv =TRUE ,Colv = TRUE, dendrogram = c("row"), scale ="row",ColSideColors =c("red","black")[as.factor(c.new)], labRow = colnames(Y.new.good), labCol = NA, main = 'VIJVER (2002) \n Testing Set  \n 12 Genes High Rank Genes', col =hmcols, cexRow =0.40, trace ="none", key.title = "Color Key")
legend("topright", legend = c("Good Prognosis","Bad Prognosis"),fill = c("Red","Black"), cex = 0.4)


####### Plotting the Feature Importance

heatmapdata <- as.matrix(cbind(cluste,t(final.betahat)))
rownames(heatmapdata) <-  colnames(Y)
colnames(heatmapdata) <- c("Overall_Y","Cluster1_time","Cluster2_time")
hmcols<-colorRampPalette(c("white","black"))(128)
heatmap.2(heatmapdata ,dendrogram="none", margins=c(6,10),col = hmcols, main = "Feature Selection ", cexCol = 0.85, cexRow = 0.7, Rowv = FALSE, Colv = FALSE) 


###########################################
pdf('VIJVER.pdf')
hmcols<-colorRampPalette(c("green","black","red"))(256)
heatmap.2(x = t(Y.good), Rowv =TRUE ,Colv = TRUE, dendrogram = c("both"), scale ="row",ColSideColors =c("red","black")[as.factor(c.final)], labRow = colnames(Y.good), labCol = NA, main = 'VIJVER (2002) \n Training Set (n =175) \n 25 Genes High Rank Genes', col =hmcols, cexRow =0.40, trace ="none", key.title = "Color Key")
legend("topright", legend = c("Good Prognosis","Bad Prognosis"),fill = c("Red","Black"), cex = 0.4)

heatmap.2(x = t(Y.new.good), Rowv =TRUE ,Colv = TRUE, dendrogram = c("both"), scale ="row",ColSideColors =c("red","black")[as.factor(c.new)], labRow = colnames(Y.new.good), labCol = NA, main = 'VIJVER (2002) \n Testing Set (n =144) \n 25 Genes High Rank Genes', col =hmcols, cexRow =0.40, trace ="none", key.title = "Color Key")
legend("topright", legend = c("Good Prognosis","Bad Prognosis"),fill = c("Red","Black"), cex = 0.4)

######## Let us compare the survival curves #####################################
surv.ob <- Surv(exp(time.new),censoring.new)
survfit.dpmm <- survdiff(surv.ob ~ c.new)
survfit.kk <- survdiff(surv.ob ~ label.test)
surv.fit <- survfit(surv.ob ~ c.new)
p5 <- ggsurv(surv.fit, main = " Test Data \n Stratified with Prediction Labels \n p-value 0.06", surv.col = c("black","red"))  + ggplot2::guides(linetype = FALSE) + ggplot2::scale_colour_discrete(name = 'Classes',breaks = c(1,2),labels = c('Good Prognosis', 'Bad Prognosis'))
p5
########### Feature Selection ##############################################
heatmapdata <- as.matrix(cbind(cluste,t(final.betahat)))
rownames(heatmapdata) <-  colnames(Y)
colnames(heatmapdata) <- c("Overall_Y","Cluster1_time","Cluster2_time")
hmcols<-colorRampPalette(c("white","black"))(128)
heatmap.2(heatmapdata ,dendrogram="none", margins=c(6,10),col = hmcols, main = "Feature Selection ", cexCol = 0.85, cexRow = 0.7, Rowv = FALSE, Colv = FALSE) 

dev.off()


##################################################################################
################### Visualizing the Clustering ###################################
##################################################################################

######### TRAINING DATA ############################################################
########## VIJVER TRAINING DATA  ########################################################
pdf('VantVeer.pdf')
hmcols<-colorRampPalette(c("green","black","red"))(256)
heatmap.2(x = t(Y.vijver.train), Rowv =TRUE ,Colv = TRUE, dendrogram = c("both"), scale ="row",ColSideColors =c("red","black")[as.factor(c.true)], labRow = colnames(Y.vijver.train), labCol = NA, main = 'VIJVER (2002) Training Set \n with 70 Vant Veer (2001) Signature', col =hmcols, cexRow =0.40, trace ="none", key.title = "Color Key")
legend("topright", legend = c("Good Prognosis","Bad Prognosis"),fill = c("Red","Black"), cex = 0.4)
dev.off()




####################################################################################################################################################
####################################################################################################################################################
########################### MAKING FUNCTIONAL /BIOLOGICAL SENSE OF THE DPMM SIGNATURE ##############################################################
####################################################################################################################################################
####################################################################################################################################################

anno.sig <- fData(vanDeVijver)[(rownames(fData(vanDeVijver)) %in% colnames(Y) ),]
anno.genes <- anno.sig[,2]
anno.probes <- rownames(anno.sig)[!is.na(anno.sig[,2])]

### Look ONLY at probes which have annotation ##############
heatmap.annot <- heatmapdata[rownames(heatmapdata) %in% anno.probes,]

###### Plotting Heatmaps for Annotated Genes
genes.heatmap <- (anno.sig[match(rownames(heatmap.annot),rownames(anno.sig)),2])

### Plotting the actual heatmap
hmcols<-colorRampPalette(c("white","black"))(256)
pdf('FeatureSelection_BreastCancer.pdf')
heatmap.2(x = heatmap.annot, Rowv =TRUE ,Colv = FALSE, dendrogram = c("row"), scale ="none", labRow = genes.heatmap, labCol = c("Overall","Bad_Prognosis","Good_Prognosis"), main = ' \n Feature Selection \n 42 Annotated DPMM signature', col =hmcols, cexRow =0.40, cexCol =0.50,trace ="none", key.title = "Color Key")
dev.off()
