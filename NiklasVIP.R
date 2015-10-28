#### This file tries to find cluster membership of the Relapse patients by using parameters learnt from 
#### my model
rm(list =ls())

library("Biobase")
library('survival')

## Load the parameters of the model learnt from the data
load( file = '/home/bit/ashar/ExpressionSets/VIPdataset/VIPDP.RData')
library(xlsx)
table1 <- read.xlsx(file ='/home/bit/ashar/ExpressionSets/VIPdataset/Niklas/15-07-10 VIP Sample Collection Data_Clip Project.xlsx', sheetIndex =1)


## Just take the periphery samples and the recurrent samples
index.1 <- (table1[,5] == '91R' | table1[,5] == '127R' | table1[,5] == '129R' | table1[,5] == '132R' | table1[,5] == '143R') & (table1[,7] != 'center') 

table.rel <- table1[index.1,]

## Removing some NA samples

rel.tab <- table.rel[1:16,]

### Apparently some of the samples are already in the original data set
index.table <- which( rel.tab[,4] %in% pheno.sig[,3] == TRUE)
index.sig <-  which(pheno.sig[,3] %in% rel.tab[,4] == TRUE)

sub.table <- rel.tab[index.table,c(4,5,32)]
sub.phen <- pheno.sig[index.sig,c(3,31)]

## The conclusion is that 129 and 143 patients have the same .CELL files and have PFS1 and PFS2.
#### Need a clarification of this fact


###### PREPARING THE NEW DATA SET

Y.new <- eset.ec.norm[featureNames(eset.ec.norm) %in% colnames(Y),sampleNames(eset.ec.norm) %in% as.character(as.matrix(rel.tab[,4]))] 

new.table <- rel.tab[match(sampleNames(Y.new),rel.tab[,4]),]
sampleNames(Y.new) == new.table[,4]
Ynew <- t(as.matrix(Y.new))
### Preparing the new time points PFS2
time.new <- log(as.numeric(as.matrix(new.table[,32])) * 30)




#### We have 16 new points ####### 
### We would like to know what are the probabilities to belong to clusters
#### Just some vizualization
pc <- prcomp((Y))
pc.pred <- predict(pc,newdata = (Y))
plot(pc.pred[,1], pc.pred[,2], pch = 19, col =c)




Y.grand <- rbind(Y,Ynew)
That.new <- exp(time.new)
t.grand <- c(30*That,That.new)
c.grand <- c(c,as.vector(as.matrix(new.table[,5])))
ct.grand <- c(c, rep(3,16))
censoring.new <- rep(1,16)
censoring.grand <- c(censoring,censoring.new)


## Drawing survival curves
surv.days <- Surv((t.grand),censoring.grand)
fit <- survfit(surv.days ~ ct.grand)
pdf("PFS2survivalcurves.pdf")
plot(fit, col = c("black","magenta","green"), xlab = "PFS and PFS2 (green) (days) \n Median PFS black 259d  Median PFS Red 198d")
title("Survival Curves")
dev.off()

names <- c(as.character(pheno.sig[,4]),as.character(as.matrix(new.table[,5])))
pdf("PatientswithNewLabels.PCA")
pc <- prcomp(Y.grand)
pc.pred <- predict(pc,newdata = Y.grand)
pc.pred[,2] <- -pc.pred[,2]
plot(pc.pred[,1], pc.pred[,2], pch = 19, col =as.factor(c.grand))
text(pc.pred[,1], pc.pred[,2],labels=names, cex= 0.3, pos=2)
dev.off()


source('predictchineseAFT.R')
poste <- predictchineseAFT(Y.new, time.new, N.new, c.list ,Y, mu.list, S, alpha, beta0.list, betahat.list, sigma2.list, lambda2.list, tau2.list, K, D, N)
post <- cbind(as.vector(as.matrix(new.table[,5])), posteriorprob)  
write.csv(post, file = 'posterior.csv')



### Plotting the clusters also with the Patient Names
pdf("PatinetsLabelPCA.pdf")
pc <- prcomp(Y)
pc.pred <- predict(pc,newdata = Y)
plot(pc.pred[,1], pc.pred[,2], pch = 19, col = c("black","magenta")[c])
text(pc.pred[,1],pc.pred[,2],labels= pheno.sig[,4], cex= 0.2, pos=2)
dev.off()

