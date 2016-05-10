########## This file prepares the TRaining and Test Data for the VIJVER Case ###############
########## Our own Feature Selection is Performed ##########################################
########## Our Feature Selection was done on Metastis or NOT ###############################
##### Use just the Training Data ###########################################################
rm(list =ls())
library('breastCancerNKI')
library('seventyGeneData')
library('gdata')
library('cancerdata')

data(nki)
data(vantVeer)
data(vanDeVijver)
data(VIJVER1)



######################## 70 GENE SIGNATURE #################################################
###  Geting Features Which Actually are Part of the 70 Gene Signature 
gns231 <- read.xls(xls ='/home/bit/ashar/ExpressionSets/ONE_VIEW/breastcancer/415530a-s9.xls', skip=0, header=TRUE, stringsAsFactors=FALSE)
###Remove special characters in the colums header, which are due to white spaces present in the Excel files
colnames(gns231) <- gsub("\\.\\.", "", colnames(gns231))
###Remove GO annotation
gns231 <- gns231[, -grep("sp_xref_keyword_list", colnames(gns231))]
###Reorder the genes in decreasing order by absolute correlation
gns231 <- gns231[order(abs(gns231$correlation), decreasing=TRUE),]
###Select the feature identifiers corresponding to the top 231 and 70 genes
gns231$genes231 <- TRUE
gns231$genes70 <- gns231$accession %in% gns231$accession[1:70]
######## Signature Probes ##########################################
signature.probes <- gns231$accession[gns231$genes70]
signature.names <- gns231$gene.name[gns231$genes70]


#################### Getting the Appropriate phenoData #################################
########################################################################################

pheno.vijver <- read.csv('/home/bit/ashar/ExpressionSets/ONE_VIEW/breastcancer/VIJVERpheno.csv')
train.index <- pheno.vijver$Label_Traing_and_Validation == 'Training'
test.index <- pheno.vijver$Label_Traing_and_Validation == 'Validation'



########## Preparing the Data Matrix ######################################
Y.full <- t(exprs(vanDeVijver[featureNames(vanDeVijver) %in% featureNames(VIJVER1),]))

time.pre <- pheno.vijver$Follow_up_time_or_metastasis 
censoring.pre <- pheno.vijver$event_metastasis

c.vijver <- pheno.vijver$X70_genes
levels(c.vijver)[1:2] <- c(1:2)


########### Defining the data for the VIJVER case ###########################
Y.vijver.train <- Y.full[train.index, colnames(Y.full) %in% signature.probes]
Y.vijver.test <-  Y.full[test.index, colnames(Y.full) %in% signature.probes]


############### Defining Our data without prefiltering #########################
Y.ashar.pretrain <- Y.full[train.index, ]
Y.ashar.pretest <-  Y.full[test.index,]


c.true <- c.vijver[train.index]
c.true.new <- c.vijver[test.index]

time <- time.pre[train.index]
censoring <- censoring.pre[train.index]

time.new <- time.pre[test.index]
censoring.new <- censoring.pre[test.index]

######## Prefiltering of the Genes ############################### ###########################
######### Using Univariate Cox's Model for Ranking of the Survival Power of the Genes ########
######## Using T-test (or K way Anova) for Ranking of Clustering Power of the Genes ###########
# surv.obj <- Surv(time.pre,censoring.pre)
# coeff.sig <- c(0)
# 
# pvalue.sig <- c(0)
# pvalue.anova <- c(0)
# 
# for ( i in 1:ncol(Y.full)){
#   q1 <- unlist(summary(coxph(surv.obj ~ Y.full[,i], data = as.data.frame(Y.full))))
#   pvalue.sig[i] <- q1$logtest.pvalue 
#   q2 <- unlist(summary(aov(as.numeric(c.vijver) ~ as.numeric(Y.full[,i]), data=as.data.frame(Y.full)) ))
#   pvalue.anova[i] <- as.numeric(q2[9]) 
#   
# }

######## Prefiltering of the Genes ############################### ###########################
######### Using Univariate Cox's Model for Ranking of the Survival Power of the Genes ########
######## Using T-test (or K way Anova) for Ranking of Clustering Power of the Genes ###########
######## The Class Labels ARE WHETHER METASTATIS OCCURED OR NOT ###############################
surv.obj <- Surv(time,censoring)
coeff.sig <- c(0)

pvalue.sig <- c(0)
pvalue.anova <- c(0)

for ( i in 1:ncol(Y.ashar.pretrain)){
  q1 <- unlist(summary(coxph(surv.obj ~ Y.ashar.pretrain[,i], data = as.data.frame(Y.ashar.pretrain))))
  pvalue.sig[i] <- q1$logtest.pvalue 
  q2 <- unlist(summary(aov(as.numeric(censoring) ~ as.numeric(Y.ashar.pretrain[,i]), data=as.data.frame(Y.ashar.pretrain)) ))
  pvalue.anova[i] <- as.numeric(q2[9]) 
  
}

###### Adjusting p-values for Multiple Test Correction
pvalue.sig.adj <- p.adjust(pvalue.sig, method = "fdr")
pvalue.anova.adj <-  p.adjust(pvalue.anova, method = "fdr")



sum(((pvalue.anova.adj < 0.04) & (pvalue.sig.adj < 0.05)) + 0)

signature <-  ((pvalue.anova.adj < 0.04) & (pvalue.sig.adj < 0.05))

probes.signature <- colnames(Y.ashar.pretrain)[signature]



########## Creating training and Test Data ##########################
Y <- Y.ashar.pretrain[,colnames(Y.ashar.pretrain) %in% probes.signature ]
Y.new <- Y.ashar.pretest[,colnames(Y.ashar.pretest) %in% probes.signature ] 

################################################################################################
########### Apparently Just 4 Genes Are common between the Vant Veer Signature and Ours #########
#################################################################################################
pc <- prcomp(Y)
pc.pred <- predict(pc,newdata = Y)
p<- ggplot(as.data.frame(pc.pred), aes(x=pc.pred[,1], y= pc.pred[,2], colour= as.factor(c.true))) + ggtitle(" Van't Veer Clustering") + geom_point(shape=19) + labs(y = "PC1", x = "PC2", colour = "Classes") 

pc <- prcomp(Y.new)
pc.pred <- predict(pc,newdata = Y.new)
p<- ggplot(as.data.frame(pc.pred), aes(x=pc.pred[,1], y= pc.pred[,2], colour= as.factor(c.true.new))) + ggtitle(" Van't Veer Clustering") + geom_point(shape=19) + labs(y = "PC1", x = "PC2", colour = "Classes") 

surv.ob <- Surv(time,censoring)
surv.fit <- survfit(surv.ob ~ c.true)
logrank <- survdiff(surv.ob ~ c.true)
p2 <- ggsurv(surv.fit, main = " DPMM \n Kaplan Meier Estimators") + ggplot2::guides(linetype = FALSE) + ggplot2::scale_colour_discrete(name = 'Classes',breaks = c(1,2),labels = c('1', '2'))


relev <- list('Y.train' =Y, 'time.train' = time,'censoring.train'= censoring, 'Y.test'= Y.new,  'time.test' =time.new, 'censoring.test'= censoring.new, 'c.true' = c.true, 'c.true.new'= c.true.new, 'Y.vijver.train'= Y.vijver.train, 'Y.vijver.test'= Y.vijver.test)
save(relev, file = 'DataVijverOWNsignatureBEST.RData')
