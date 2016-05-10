############ This file takes in training and testing data For the VIP data study #################
#################################################################################################
rm(list =ls())

setwd("~/Dropbox/Code/DPmixturemodel/DPplusAFT")
################################################
######## Load the Data from the folder #########
load("/home/bit/ashar/ExpressionSets/ONE_VIEW/VIPdataset/DPMM-FactorAnalysis/DataVIPOWNsignaturePLUSMeanSamples.RData")

###########################################################################
##### This file Runs the Actual VIP data set ##############################
###########################################################################
Y <- relev$Y.train
time <- relev$time.train
censoring<- rep(1, dim(Y)[1])

N <- nrow(Y)
D <- ncol(Y)


############################# PARAMETERS for GIBB's SAMPLING ####
iter = 100
iter.burnin = 100
iter.thin  = 2
k = 2
Nps = 50

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
p1 <- ggplot(as.data.frame(pc.pred), aes(x=pc.pred[,1], y= pc.pred[,2], colour= as.factor(c.final))) + ggtitle(" DPMM Clustering \n VIP patients (45) \n 60 Gene DPMM signature") + geom_point(shape=19) + labs(y = "PC1", x = "PC2", colour = "Classes") 


####Fitting Survival Curves
surv.ob <- Surv(exp(time),censoring)
surv.fit <- survfit(surv.ob ~ c.final)
logrank <- survdiff(surv.ob ~ c.final)
p5 <- ggsurv(surv.fit, main = " Kaplan Meier Estimators \n PFS for Peripheral Cells") + ggplot2::guides(linetype = FALSE) + ggplot2::scale_colour_discrete(name = 'Classes',breaks = c(1,2),labels = c('Mean: 569d', 'Mean 309d'))


####### Calculating the rank
rank <- matrix(0, nrow = N, ncol =Nps)
Nps =25
for (j in 1:Nps){
  for ( i in 1:N){
    rank[i,j] <- dMVN(as.vector(t(Y[i,1:D])), mean = mu.list[[j]][1,1:D], Q= S.list[[j]][1,1:D,1:D], log = TRUE) - dMVN(as.vector(t(Y[i,1:D])), mean = mu.list[[j]][2,1:D], Q= S.list[[j]][2,1:D,1:D], log = TRUE) 
  }
}

avg.rank <- apply(rank,1,mean)
order.zo <- range01(avg.rank)

order.train <- sort(order.zo,index.return = TRUE, decreasing = TRUE)
Y.order <- Y[order.train$ix,]
c.final.order <- c.final[order.train$ix]
rownames(Y.order) <- rownames(Y)[order.train$ix]

c.ultimate <- as.factor(c.final.order)
levels(c.ultimate) <- c(1,2)


pdf('DPMM-mean.pdf')
hmcols<-colorRampPalette(c("green","black","red"))(256)
heatmap.2(x = t(Y.order), Rowv =TRUE ,Colv = FALSE, dendrogram = c("row"), scale ="row",ColSideColors =c("red","black")[c.ultimate], labRow = colnames(Y), labCol = rownames(Y.order), main = 'DPMM clustering on VIP samples \n Based on P Cells (averaged) (45 P) \n 60 Gene DPMM signature', col =hmcols, cexRow =0.40, trace ="none", key.title = "Color Key")
legend("topright", legend = c("Cluster 1","Cluster 2"),fill = c("Red","Black"), cex = 0.4)
p1
dev.off()

pdf("VIPsurvival-curves.pdf")
plot(surv.fit, main = 'PF Survival Curves for VIP Patients \n p = 0.06 \n Clustering Based on P Samples', xlab= 'Mean PFS Cluster 1: 309d \n Mean PFS Cluster2: 569d')
dev.off()


