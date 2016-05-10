########## This file prepares the Test Data for the VIJVER Case ###############
######### All Lymph negative patients were classified as Traning case
######## All lymph positive patients were classified as Testing Case

rm(list =ls())
library('breastCancerNKI')
library('gdata')
data(nki)

###  Geting Features Which Actually are Part of the 70 Gene Signature 
gns231 <- read.xls(xls ='/home/bit/ashar/ExpressionSets/breastcancer/415530a-s9.xls', skip=0, header=TRUE, stringsAsFactors=FALSE)
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

########## Creating Training and Testing Data Set ##########################################
####### The Model is Based on the Time to Distant Metastatis Free

ind.exl <- is.na(pData(nki)$e.dmfs)

nod.negative <- pData(nki)$node == 0

nod.positive <- pData(nki)$node == 1


train.index <- nod.negative & !ind.exl

test.index <- nod.positive & !ind.exl

########## Creating training and Test Data ##########################

data.breast <- nki[featureNames(nki) %in% signature.probes, ]


Y.train <- t(exprs(data.breast))[train.index,]
time.train <- pData(data.breast)$t.dmfs[train.index]
censoring.train <- pData(data.breast)$e.dmfs[train.index]


Y.test <- t(exprs(data.breast))[test.index,]
time.test <- pData(data.breast)$t.dmfs[test.index]
censoring.test <- pData(data.breast)$e.dmfs[test.index]


relev <- list('Y.train' =Y.train, 'time.train' = time.train,'censoring.train'= censoring.train, 'Y.test'= Y.test,  'time.test' =time.test, 'censoring.test'= censoring.test)
save(relev, file = 'DataVijver.RData')
