### The File takes 218 samples ###############################################################################
#### This file does Prefiltering also based on those genes which discriminate PvsZ and ARE releated to survival
##### This script CLUSTERS Periphery Cells ON THE PvsZ signature and PFS survival  ###########################
library('affy')
library('xlsx')
library('limma')
library('survival')
rm(list =ls())

## load data ##
load("/home/abidata/Dinis/VIP/Main.Data/ExpressionConsole.normalised.RData")


## The phenoData
tab <- read.xlsx(file  = '/home/abidata/Dinis/VIP/Sample.Info/15-07-14 VIP Sample Collection PvsZ DEG.xlsx', sheetIndex =1)
list.patients <- tab[,2]


### Only those patients who have a NON NA is PFS
pheno  <- pData(eset.ec.norm) 
pheno.ready <- pheno[!is.na(pheno$PFS),]
list.patients.final <- list.patients[(list.patients %in% pheno.ready[,3])]


## Include ONLY those features that have correspoding annotation ~ 27148 out of 70000
## Include ONLY those Samples which HAVE PFS and PvsZ Annotation
exprs <- eset.ec.norm[featureNames(eset.ec.norm) %in% rownames(anno.GENENAME),sampleNames(eset.ec.norm) %in% list.patients.final]



#############################################################################################################################
######## SIGNATURE WAS IDENTIFIED WITH THOSE FEATURES WHICH DISTINGUISHED PvsZ Clustering AND HAD SURvival Information ######
#############################################################################################################################
## LIMMA PREFILTERING BASED ON PvsZ #################################################
pheno <-pData(exprs)
pheno$Case <- as.factor(as.matrix(pheno$Case))
pheno$Topo <- as.factor(as.matrix(pheno$Topo))
pheno$Class <- as.factor(as.matrix(pheno$Class))


mm <- model.matrix( ~ 0 + Topo + Case  + Class , pheno)

cc <- colnames(mm)

cc <-  gsub("_","", cc)

cc <- gsub("-","", cc)

colnames(mm) <- cc

cont.matrix <- makeContrasts(ZvsP = Topocenter - Topoperiphery, levels= mm)

fit <- lmFit(exprs, design = mm)
fit2 <- contrasts.fit(fit, cont.matrix)
fit3 <- eBayes(fit2)

DEG_table <- topTable(fit3, adjust="BH", coef='ZvsP', number = Inf)

list <-  (DEG_table)[DEG_table$adj.P.Val < 0.05,1]

########### Fitting Univariate Cox's Models ######################################################

Xsig.ready <- as.data.frame(t(exprs(exprs[featureNames(exprs) %in% list,])))



######## Getting survival times (PFS) and status #####################
Time <- (as.numeric(as.matrix(pheno$PFS)))
status <- rep(1, dim(pheno)[1]) # 1 because the relapse always occurs
surv.obj <- Surv(Time, status)

##### FITTING A UNIVARIATE COX REGRESSION MODEL ################################
pvalue.sig <- c(0)
pvalue.adj <- c(0)

for ( i in 1:ncol(Xsig.ready)){
  q <- unlist(summary(coxph(surv.obj ~ Xsig.ready[,i], data = Xsig.ready)))
  pvalue.sig[i] <- q$logtest.pvalue  
}
pvalue.adj <- p.adjust(pvalue.sig, method = "fdr")
list.final <- colnames(Xsig.ready)[(pvalue.adj < 0.2)]


##############################################################################################
###############################################################################################
###### OUR SIGNATURE ###########################################################################
signature <- list.final

############# Let's Compare Our Signature with the Master List to see if we have any hits ##############
master <-  read.xlsx(file = '/home/abidata/Dinis/VIP/Main.Data/Master_File_FINAL.xlsx', sheetIndex =1)

signature %in% master[,2]

##########################################
exprs.final <- exprs[featureNames(exprs) %in% signature,]

annot <- anno.GENENAME[anno.GENENAME[,1] %in% signature,]

featureNames(exprs.final) ==  annot[,1]

###################################################
#### Prepare the data set for clustering ##########

Y <- t(exprs(exprs.final[,pData(exprs.final)$Topo == 'periphery']))
colnames(Y) <-  annot[,3]

time <- as.numeric(as.matrix(pData(exprs.final[,pData(exprs.final)$Topo == 'periphery'])$PFS))
censoring <- rep(1, nrow(Y))


Y.info <- pData(exprs.final[,pData(exprs.final)$Topo == 'periphery'])

relev <- list('Y.train' =Y, 'time.train' = time,'censoring.train'= censoring, 'Y.info' = Y.info)

save(relev, file = 'DataVIPOWNsignature.RData')

