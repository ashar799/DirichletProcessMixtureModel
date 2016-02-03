######## This program visualizes the code that was used to generate plots for the ABC conference ########
time.dataframe <- list()
time.dataframe$time <-  time
time.dataframe$class <-  c.true
t <- as.data.frame(time.dataframe)


p1 <- ggplot(as.data.frame(Y), aes(x=Y[,1], y= Y[,2], colour= c.true)) + ggtitle("Relevant Molecular Features") + geom_point(shape=19) + labs(y = "Feature 1", x = "Feature 2", colour = "Classes") 
p2 <- ggplot(as.data.frame(Y), aes(x=Y[,3], y= Y[,4], colour= c.true)) + ggtitle("Irrelevant Molecular Noisy Features") + geom_point(shape=19) + labs(y = "Feature 3", x = "Feature 4", colour = "Classes") 
p3 <- ggplot(data = t) + geom_boxplot(aes(y = t$time, x= t$class, fill = (t$class))) + ggtitle("Simulated Survival Times") + labs(y = "Survival Times", x = "Class") + theme(legend.title = element_blank())

############### Taking PCA plots
pc <- prcomp(Y)
pc.pred <- predict(pc,newdata = Y)
p4 <- ggplot(as.data.frame(pc.pred), aes(x=pc.pred[,1], y= pc.pred[,2], colour= as.factor(c.final))) + ggtitle(" DPMM Clustering") + geom_point(shape=19) + labs(y = "PC1", x = "PC2", colour = "Classes") 

####Fitting Survival Curves
surv.ob <- Surv(time,censoring)
surv.fit <- survfit(surv.ob ~ c.final)
logrank <- survdiff(surv.ob ~ c.final)
p5 <- ggsurv(surv.fit, main = " DPMM \n Kaplan Meier Estimators") + ggplot2::guides(linetype = FALSE) + ggplot2::scale_colour_discrete(name = 'Classes',breaks = c(1,2),labels = c('1', '2'))

##### Feature Selection
df <- expand.grid(x = 0:1, y = 0:3)
df$Prob <- melt(heatmapdata)[,2]
df$variable <- c("Cluster1","Cluster2")
df$Name <- c("Feature1","Feature1","Feature2","Feature2","Feature3","Feature3","Feature4","Feature4")
p6 <- ggplot(df, aes(variable, Name)) + geom_tile(aes(fill = Prob)) + scale_fill_gradient(low = "white",high = "black") + ggtitle(" DPMM \n Posterior Probabilities \n Feature Selection") + labs(y = "Features", x = "Clusters")
  
###### Plot the results
source('multiplot.R')
multiplot(p1, p2, p3, cols=3)


###################################################################################################
########## SIMULATIONS AND COMPARISON WITH THE GROUND TRUTH #######################################
###################################################################################################


## Analyze the 6 simulations
rm(list =ls())
setwd('/home/bit/ashar/ExpressionSets/Simulations/D10K30C')

### Rand-Indices
rand.km <- rep(0,6)
rand.flx <- rep(0,6)
rand.DPGMM.train <- rep(0,6)
rand.DPGMM.pred <- rep(0,6)


#### C- Indices
cindex.flx.cox.sim <- rep(0,6)
cindex.km.cox.sim <- rep(0,6)
cindex.km.pcox.sim <- rep(0,6)
cindex.km.aft.sim <- rep(0,6)
cindex.km.paft.sim <- rep(0,6)
cindex.DPGMM.train <- rep(0,6)
cindex.DPGMM.predicted <- rep(0,6)

### P-value for Closeness
pvalue.pred <- rep(0,6)

load(file = 'D10K30C.RData')


rand.km[1] <- gr.km.rand.final
rand.flx[1] <- gr.flx.rand.final
rand.DPGMM.train[1] <- mean(final.rand)

cindex.flx.cox.sim[1] <- cindex.flx.cox.final 
cindex.km.cox.sim[1] <- cindex.km.cox.final
cindex.km.pcox.sim[1] <- cindex.km.pcox.final
cindex.km.aft.sim[1] <- cindex.km.aft.final
cindex.km.paft.sim[1] <- cindex.km.paft.final
cindex.DPGMM.train[1] <- mean(cindex.final)


pvalue.pred[1] <-  predicted.close.pvalue

for (i in 1:5){
  filename <- paste("D10K30Crep",i,".RData",sep = "")
  load(filename)
  rand.km[i+1] <- gr.km.rand.final
  rand.flx[i+1] <- gr.flx.rand.final
  rand.DPGMM.train[i+1] <- mean(final.rand)
  
  cindex.flx.cox.sim[i+1] <- cindex.flx.cox.final 
  cindex.km.cox.sim[i+1] <- cindex.km.cox.final
  cindex.km.pcox.sim[i+1] <- cindex.km.pcox.final
  cindex.km.aft.sim[i+1] <- cindex.km.aft.final
  cindex.km.paft.sim[i+1] <- cindex.km.paft.final
  cindex.DPGMM.train[i+1] <- mean(cindex.final)

}


### Plot some nice Box Plots

cindex.comp <- cbind(cindex.flx.cox.sim, cindex.km.cox.sim,cindex.km.pcox.sim,cindex.km.aft.sim,cindex.km.paft.sim,cindex.DPGMM.train)
colnames(cindex.comp) <- c("FLx+Cox","KM+Cox","KM+PenCox","KM-AFT","KM-PenAFT","DPMM")

randindex.comp <- cbind(rand.km,rand.flx,rand.DPGMM.train )
colnames(randindex.comp) <- c("KM","Flexmix","DPMM")


cind <- melt(cindex.comp)
p1 <- ggplot(data = as.data.frame(cind)) + geom_boxplot(aes(y = cind$value, x= as.factor(cind$X2), fill = (cind$X2))) + ggtitle("Simulation Results C-Index") + labs(y = "C-Index", x = "Models") + theme(legend.title = element_blank())

rind <- melt(randindex.comp)
p2 <- ggplot(data = as.data.frame(rind)) + geom_boxplot(aes(y = rind$value, x= as.factor(rind$X2), fill = (rind$X2))) + ggtitle("Simulation Adj. Rand-Index") + labs(y = "Adjusted Rand-Index", x = "Models") + theme(legend.title = element_blank())
source('multiplot.R')
multiplot(p1, p2)

########### Visualizing Verhaak Data Set
############### Taking PCA plots
pc <- prcomp(Y)
pc.pred <- predict(pc,newdata = Y)
p1<- ggplot(as.data.frame(pc.pred), aes(x=pc.pred[,1], y= pc.pred[,2], colour= as.factor(c.final))) + ggtitle(" DPMM Clustering \n on Verhaak Data Set \n Pathway Features") + geom_point(shape=19) + labs(y = "PC1", x = "PC2", colour = "Classes") 

####Fitting Survival Curves
surv.ob <- Surv(exp(time),censoring)
surv.fit <- survfit(surv.ob ~ c.final)
logrank <- survdiff(surv.ob ~ c.final)

####### Calculate the Mean Survival times
p2<- ggsurv(surv.fit, main = " DPMM \n Kaplan Meier Estimators \n LogRank p-value 0.12") + ggplot2::guides(linetype = FALSE) + ggplot2::scale_colour_discrete(name = 'Median Survival Time',breaks = c(1,2,3,4),labels = c('383d','372d','399d','357d')) + ggplot2::xlab('Time in Days')
pdf('Verhaak.pdf')
multiplot(p1, p2, cols =2)
dev.off()
######## Vusualizing the VIP data #############
pc <- prcomp(Y)
pc.pred <- predict(pc,newdata = Y)
p1<- ggplot(as.data.frame(pc.pred), aes(x=pc.pred[,1], y= pc.pred[,2], colour= as.factor(c.final))) + ggtitle(" DPMM Clustering \n on Periphey VIP Data Set \n DE Gene Features") + geom_point(shape=19) + labs(y = "PC1", x = "PC2", colour = "Classes") 

####Fitting Survival Curves
surv.ob <- Surv(30*(time),censoring)
surv.fit <- survfit(surv.ob ~ c.final)
logrank <- survdiff(surv.ob ~ c.final)
p2<- ggsurv(surv.fit, main = " DPMM \n Kaplan Meier Estimators \n LogRank p-value 0.06") + ggplot2::guides(linetype = FALSE) + ggplot2::scale_colour_discrete(name = 'Median PFS Time',breaks = c(1,2),labels = c('259d','198d')) + ggplot2::xlab('Time in Days')
pdf('VIP.pdf')
multiplot(p1, p2, cols =2)
dev.off()
