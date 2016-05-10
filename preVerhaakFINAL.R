#######################################################################################################################
####################################### This file prepares the Verhaak classes #########################################
########################################################################################################################

######## The Original Verhaak Data Set is split up into Training and Testing ###########################################
########################################################################################################################


rm(list =ls())
############ Load Training set ###################################
load('/home/bit/ashar/ExpressionSets/ONE_VIEW/Verhark/OriginalVerhaakData.RData')
train.prelim <- t(exprs.norm)
mode(train.prelim) <- "numeric"



###### Load Verhaark Gene Signature #################################
signature.prelim <- read.xlsx('/home/bit/ashar/ExpressionSets/ONE_VIEW/Verhark/ClaNC840_centroids.xls', sheetIndex =1)
signature <- signature.prelim[3:nrow(signature.prelim),1]
signature <- signature[signature %in% rownames(train.prelim)]



##### Load Pheno Data for the training data set #####################
load('/home/bit/ashar/ExpressionSets/ONE_VIEW/Verhark/phenoVerhaark.RData')


### checking the pattern in training data #####################
hmcols<-colorRampPalette(c("green","black","red"))(512)
heatmap.2(x = (train.prelim[as.character(signature),]), scale = "row",Rowv =FALSE ,Colv = FALSE,ColSideColors =c("blue","red","green","purple")[pheno[,4]], labRow = signature, labCol = NA, main = ' \n Training Set \n ', col =hmcols, cexRow =0.40, trace ="none", key.title = "Color Key")

## Checking Separability on PCA
pc <- prcomp(t(train.prelim[as.character(signature),]))
pc.pred <- predict(pc,newdata = t(train.prelim[as.character(signature),]))
plot(pc.pred[,1], pc.pred[,2], pch = 19,col = pheno[,4])


############################################################################################################################################################
##### Using the 840 Gene Signature For training and Testing Set and Filtering out Genes which have low survival information #################################
#############################################################################################################################################################
N.total <- ncol(train.prelim)
index.train <- seq(1, N.total, by =2)
index.test <- seq(2, N.total, by =2)

Y.train.prelim <- train.prelim[,index.train]
Y.test.prelim <- train.prelim[,index.test]
pheno.train <- pheno[index.train,]
pheno.test <- pheno[index.test,]


### checking the pattern in training data #####################
hmcols<-colorRampPalette(c("green","black","red"))(512)
heatmap.2(x = (Y.train.prelim[as.character(signature),]), scale = "row",Rowv =FALSE ,Colv = FALSE,ColSideColors =c("blue","red","green","purple")[pheno.train[,4]], labRow = signature, labCol = NA, main = ' \n Training Set \n ', col =hmcols, cexRow =0.40, trace ="none", key.title = "Color Key")

### checking the pattern in testing data #####################
hmcols<-colorRampPalette(c("green","black","red"))(512)
heatmap.2(x = (Y.test.prelim[as.character(signature),]), scale = "row",Rowv =FALSE ,Colv = FALSE,ColSideColors =c("blue","red","green","purple")[pheno.test[,4]], labRow = signature, labCol = NA, main = ' \n Testing Set \n ', col =hmcols, cexRow =0.40, trace ="none", key.title = "Color Key")


####################################################################################################################
######################### How to reduce FEATURES i.e FEATURE SELECTION #############################################
####################################################################################################################



Y_train <- t(Y.train.prelim[as.character(signature),])
Y_test  <- t(Y.test.prelim[as.character(signature), ])


######## Getting survival times and status #####################
Time <- as.numeric(exp(pheno.train[,3]))
status <- pheno.train[,2]
surv.obj <- Surv(Time, status)

##### FITTING A UNIVARIATE COX REGRESSION MODEL ################################
pvalue.sig <- c(0)
pvalue.adj <- c(0)

X_train <- as.data.frame(t(Y.train.prelim))
for ( i in 1:ncol(X_train)){
  q <- unlist(summary(coxph(surv.obj ~ X_train[,i], data = X_train)))
  pvalue.sig[i] <- q$logtest.pvalue  
}
pvalue.adj <- p.adjust(pvalue.sig, method = "fdr")
signature.dpmm<- colnames(X_train)[(pvalue.sig < 0.005)]



############################################################################
Y.train.final <- t(Y.train.prelim[as.character(signature.dpmm),])
Y.test.final <- t(Y.test.prelim[as.character(signature.dpmm),])

###########################################################################
relev <- list('Y.train.prelim' =Y.train.prelim, 'Y.test.prelim' = Y.test.prelim, 'signature' = signature, 'signature.dpmm' = signature.dpmm, 'pheno.train' =pheno.train, 'pheno.test'= pheno.test)

save(relev, file = 'DataVerhaak.RData')


####################################### As the Genes Do not contain information on survival times, we look for pathways ###################
library(globaltest)
library(org.Hs.eg.db)
library(GO.db)
library(KEGG.db)


exprs.norm <- ExpressionSet(Y.train.prelim)
xx <- as.list(org.Hs.egALIAS2EG)
# Remove pathway identifiers that do not map to any entrez gene id
xx <- xx[!is.na(xx)]

glob <- gtKEGG( Surv(exp(pheno.train[,3]), pheno.train[,2]), exprs.norm, annotation = 'org.Hs.eg.db', multtest = "BH",probe2entrez = xx)


pathway.names <- names(glob)[1:30]


yy <- as.list(KEGGPATHID2EXTID)

######### TO GET GENE NAMES of GENES WITHIN PATHWAYS ##############################3

# For the reverse map:
# Convert the object to a list
dd <- as.list(org.Hs.egPATH2EG)
# Remove pathway identifiers that do not map to any entrez gene id
dd <- dd[!is.na(dd)]
if(length(dd) > 0){
  # The entrez gene identifiers for the first two elements of XX
  dd[1:2]
  # Get the first one
  dd[[1]]
}

path2entrez <- dd

path2entrez.subset <- path2entrez[pathway.names]

library('biomaRt')
# listDatasets(ensembl)
ensembl = useMart('ensembl')
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

names.genes <- list(0)
for( i in 1:length(pathway.names)){
  names.genes[[i]] = getBM(attributes = c("entrezgene","hgnc_symbol"), filters = "entrezgene",values =path2entrez.subset[i], mart=ensembl)
}

### Data frame FROM the EXPRESSION #################

data.exprs <- as.data.frame(exprs.norm)
data.subsets <- list(0)
for ( i in 1:length(pathway.names) ){
  
  temp.names <- names.genes[[i]][,2][c(names.genes[[i]][,2])%in% colnames(data.exprs)]
  data.subsets[[i]] <- as.data.frame(data.exprs[, temp.names]) 
}

data.combined <- matrix(0, nrow = nrow(data.exprs), ncol =length(pathway.names) )
for ( i in 1:length(pathway.names)){
  #   data.combined[,i] <- as.vector(apply(data.subsets[[i]],1, mean))
  pc <- prcomp(data.subsets[[i]])
  data.combined[,i] <- predict(pc,newdata = data.subsets[[i]])[,1]
}

rownames(data.combined) <- rownames(data.exprs)
colnames(data.combined) <- pathway.names

############ A Small Vizualization ###################
pc <- prcomp(data.combined)
pc.pred <- predict(pc,newdata = data.combined)
plot(pc.pred[,1], pc.pred[,2], pch = 19, col = pheno.train[,4])


### checking the Heatmap #######################
### checking the pattern in training data #####################
hmcols<-colorRampPalette(c("green","black","red"))(512)
heatmap.2(x = t(data.combined), scale = "row",Rowv =TRUE ,Colv = FALSE, dendrogram = "row", ColSideColors =c("blue","red","green","purple")[c.true], labRow = colnames(data.combined), labCol = NA, main = ' \n Training Set \n ', col =hmcols, cexRow =0.40, trace ="none", key.title = "Color Key")



######### Saving the Top 30 pathways #################
save(data.combined,file = '30pathwayVerhaark.RData')



