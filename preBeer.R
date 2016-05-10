##### This file reads the Lung Cancer Gene Expression as in Beer et al. 2002
##### We use the library(lungExpression) to read the data set (michigan)
rm(list =ls())
library(org.Hs.eg.db)
library(lungExpression)
data(michigan)
data(harvard)

ind.michigan <- (pData(michigan)$stage == 1)
ind.harvard <- (pData(harvard)$stage == 1) & !is.na(pData(harvard)$dead)

Y.pre <- t(exprs(michigan))[ind.michigan,]
Y.new.pre <- t(exprs(harvard))[ind.harvard,]

time <- pData(michigan)[ind.michigan,9]
censoring <- pData(michigan)[ind.michigan,10]
censoring <- as.numeric(as.matrix(censoring))

time.new <- pData(harvard)[ind.harvard,6]
censoring.new <- pData(harvard)[ind.harvard,4]
########## Use Empirical Bayes Cross Platform Normalization ########################
dat.norm.connor <-  eb(t(as.data.frame(Y.pre)), t(as.data.frame(Y.new.pre)))

#################### CREATING MATRICES ####################################################
##########################################################################################

### Creating ready data matrices
Y.ready <- t(dat.norm.connor$x)
Y.new.ready <- t(dat.norm.connor$y)

### The Class Labels from Beer et al.
c.beer <- pData(michigan)[ind.michigan,3]



######### Using Univariate Cox's Model for Ranking of the Survival Power of the Genes ####

surv.obj <- Surv(time,censoring)
coeff.sig <- c(0)
pvalue.sig <- c(0)

for ( i in 1:ncol(Y.ready)){
  q <- unlist(summary(coxph(surv.obj ~ Y.ready[,i], data = as.data.frame(Y.ready))))
  pvalue.sig[i] <- q$logtest.pvalue  
}

####### Using ANOVA to Rank the Clustering Interpretability of the Genes ############
pvalue.anova <- c(0)
for ( i in 1:ncol(Y.ready)){
  q <- unlist(summary(aov(c.beer ~ Y.ready[,i], data=as.data.frame(Y.ready)) ))
  pvalue.anova[i] <- as.numeric(q[9]) 
}

sum(((pvalue.anova < 0.000001) & (pvalue.sig < 0.05)) + 0)

signature <- (((pvalue.anova < 0.000001) & (pvalue.sig < 0.05)))


###########################################################
Y <- Y.ready[,signature]
Y.new <- Y.new.ready[,signature]

################################################################
### Principal Componenets 
pc <- prcomp(Y)
pc.pred <- predict(pc,newdata = Y)
p1 <- ggplot(as.data.frame(pc.pred), aes(x=pc.pred[,1], y= pc.pred[,2], colour= as.factor(c.beer))) + ggtitle(" Beer el.al 2002") + geom_point(shape=19) + labs(y = "PC1", x = "PC2", colour = "Classes") 

pc <- prcomp(Y.new)
pc.pred <- predict(pc, newdata = Y.new)
p2 <- ggplot(as.data.frame(pc.pred), aes(x=pc.pred[,1], y= pc.pred[,2])) + ggtitle(" Bhattacharya el.al 2000") + geom_point(shape=19) + labs(y = "PC1", x = "PC2") 


########## Saving the preprocessed data ##########################
relev <- list('Y.train' =Y, 'time.train' = time,'censoring.train'= censoring, 'Y.test'= Y.new,  'time.test' =time.new, 'censoring.test'= censoring.new)
save(relev, file = 'DataLungAdeno.RData')
