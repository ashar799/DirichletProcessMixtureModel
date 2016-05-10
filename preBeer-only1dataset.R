##### THIS JUST USES ONE DATA SET
##### This file reads the Lung Cancer Gene Expression as in Beer et al. 2002
##### We ONLY use the library(lungExpression) to read the data set (michigan)
rm(list =ls())
library(org.Hs.eg.db)
library(lungExpression)
data(michigan)
ind.michigan <- (pData(michigan)$stage == 1)

Y.pre <- t(exprs(michigan))[ind.michigan,]

time.pre <- pData(michigan)[ind.michigan,9]
censoring.pre <-as.numeric(as.matrix( pData(michigan)[ind.michigan,10]))


######### Randomly dividing the Data in Training and Testing data
Y.ready <- Y.pre[1:35,]
Y.new.ready <- Y.pre[36:67,]

time <-  time.pre[1:35]
time.new <- time.pre[36:67]

censoring <- censoring.pre[1:35]
censoring.new <- censoring.pre[36:67]

### The Class Labels from Beer et al.
c.beer <- pData(michigan)[ind.michigan,3]
c.true <- c.beer[1:35]
c.true.new <- c.beer[36:67]

#################### CREATING MATRICES ####################################################
##########################################################################################

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
  q <- unlist(summary(aov(c.true ~ Y.ready[,i], data=as.data.frame(Y.ready)) ))
  pvalue.anova[i] <- as.numeric(q[9]) 
}

sum(((pvalue.anova < 0.0001) & (pvalue.sig < 0.05)) + 0)

signature <- (((pvalue.anova < 0.0001) & (pvalue.sig < 0.05)))


###########################################################
Y <- Y.ready[,signature]
Y.new <- Y.new.ready[,signature]

################################################################
### Principal Componenets 
pc <- prcomp(Y)
pc.pred <- predict(pc,newdata = Y)
p1 <- ggplot(as.data.frame(pc.pred), aes(x=pc.pred[,1], y= pc.pred[,2], colour= as.factor(c.true))) + ggtitle(" Beer el.al 2002 \n Training") + geom_point(shape=19) + labs(y = "PC1", x = "PC2", colour = "Classes") 

pc <- prcomp(Y.new)
pc.pred <- predict(pc, newdata = Y.new)
p2 <- ggplot(as.data.frame(pc.pred), aes(x=pc.pred[,1], y= pc.pred[,2], colour= as.factor(c.true.new))) + ggtitle(" Beer el.al 2002 \n Testing") + geom_point(shape=19) + labs(y = "PC1", x = "PC2") 


########## Saving the preprocessed data ##########################
relev <- list('Y.train' =Y, 'time.train' = time,'censoring.train'= censoring, 'Y.test'= Y.new,  'time.test' =time.new, 'censoring.test'= censoring.new,'c.true'= c.true, 'c.true.new' =c.true.new)
save(relev, file = 'DataLungAdeno-only1Datset.RData')
