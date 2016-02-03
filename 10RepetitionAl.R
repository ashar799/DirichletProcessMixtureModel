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
rand.DPGMM.pred[1] <- as.numeric(predicted.cindex)

cindex.flx.cox.sim[1] <- cindex.flx.cox.final 
cindex.km.cox.sim[1] <- cindex.km.cox.final
cindex.km.pcox.sim[1] <- cindex.km.pcox.final
cindex.km.aft.sim[1] <- cindex.km.aft.final
cindex.km.paft.sim[1] <- cindex.km.paft.final
cindex.DPGMM.train[1] <- mean(cindex.final)
cindex.DPGMM.predicted[1] <- as.numeric(predicted.cindex)

pvalue.pred[1] <-  predicted.close.pvalue

for (i in 1:5){
  filename <- paste("D10K30Crep",i,".RData",sep = "")
  load(filename)
  rand.km[i+1] <- gr.km.rand.final
  rand.flx[i+1] <- gr.flx.rand.final
  rand.DPGMM.train[i+1] <- mean(final.rand)
  rand.DPGMM.pred[i+1] <- as.numeric(predicted.cindex)
  
  cindex.flx.cox.sim[i+1] <- cindex.flx.cox.final 
  cindex.km.cox.sim[i+1] <- cindex.km.cox.final
  cindex.km.pcox.sim[i+1] <- cindex.km.pcox.final
  cindex.km.aft.sim[i+1] <- cindex.km.aft.final
  cindex.km.paft.sim[i+1] <- cindex.km.paft.final
  cindex.DPGMM.train[i+1] <- mean(cindex.final)
  cindex.DPGMM.predicted[i+1] <- as.numeric(predicted.cindex)
  
  pvalue.pred[i+1] <-  predicted.close.pvalue
}


### Plot some nice Box Plots

cindex.comp <- cbind(cindex.flx.cox.sim,cindex.km.cox.sim,cindex.km.pcox.sim,cindex.km.aft.sim,cindex.km.paft.sim,cindex.DPGMM.train,cindex.DPGMM.predicted)
colnames(cindex.comp) <- c("FL-C","KM-C","KM-PC","KM-A","KM-PA","tr-DM","te-DM")
pdf("C-Index.pdf")
par(mfrow=c(2,1))
boxplot(cindex.comp[,1:3], main ="C-Index over 6 Runs",cex.lab=1.5)
boxplot(cindex.comp[,4:7], main ="C-Index over 6 Runs",cex.lab=1.5)
dev.off()

randindex.comp <- cbind(rand.km,rand.flx,rand.DPGMM.train,rand.DPGMM.pred )
colnames(randindex.comp) <- c("KM","Flex","tr-DP","te-DP")
pdf("Rand-Index.pdf")
boxplot(randindex.comp, main ="Rand Index over 6 Runs",cex.lab=1.5)
dev.off()

