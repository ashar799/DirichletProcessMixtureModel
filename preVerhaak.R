#######################################################################################################################
####################################### This file prepares the Verhaak classes #########################################
########################################################################################################################
rm(list =ls())
############ Load Training set ###################################
load('/home/bit/ashar/ExpressionSets/ONE_VIEW/Verhark/OriginalVerhaakData.RData')
train.prelim <- t(exprs.norm)
mode(train.prelim) <- "numeric"
###### Load Testing Set ############################################
test.prelim <- read.csv('/home/bit/ashar/ExpressionSets/ONE_VIEW/Verhark/CancerCell_GBM_validation.data.txt', sep = "\t")
test.prelim <- test.prelim[24:NROW(test.prelim),]


train.pre <- train.prelim[rownames(train.prelim) %in% test.prelim[,1],]
test.pre <- test.prelim[test.prelim[,1] %in% rownames(train.prelim),]
test.man <- test.pre[match(rownames(train.pre), test.pre[,1]),]

sum((rownames(train.pre) %in% test.man[,1]) +0 )
sum((test.man[,1]  %in% rownames(train.pre) )  +0 )


rownames(test.man) <- test.man[,1]
test.man <- test.man[,2:NCOL(test.man)]


#### Just Checking
sum((rownames(test.man) == rownames(train.pre)) +0)



###### Load Verhaark Gene Signature #################################
signature.prelim <- read.xlsx('/home/bit/ashar/ExpressionSets/ONE_VIEW/Verhark/ClaNC840_centroids.xls', sheetIndex =1)
signature <- signature.prelim[3:nrow(signature.prelim),1]

signature <- signature[signature %in% rownames(test.man)]

##### Load Pheno Data for the training data set #####################
load('/home/bit/ashar/ExpressionSets/ONE_VIEW/Verhark/phenoVerhaark.RData')


### Checking for Any batch Effects ###########
data.combined <- t(cbind(train.pre[as.character(signature),],test.man[as.character(signature),]))
lab <- c(rep(1,NCOL(train.pre)),rep(2, NCOL(test.man)))
pc <- prcomp(data.combined)
pc.pred <- predict(pc,newdata = data.combined)
plot(pc.pred[,1], pc.pred[,2], pch = 19,col = lab)


### checking the pattern in training data #####################
hmcols<-colorRampPalette(c("green","black","red"))(512)
heatmap.2(x = (train.pre[as.character(signature),]), scale = "row",Rowv =FALSE ,Colv = FALSE,ColSideColors =c("blue","red","green","purple")[pheno[,4]], labRow = signature, labCol = NA, main = ' \n Training Set \n ', col =hmcols, cexRow =0.40, trace ="none", key.title = "Color Key")


index.final <- c(which(cl == 3), which(cl == 1), which(cl ==2))
heatmap <- as.matrix((test.man[as.character(signature),]))
heatmap_test <- heatmap[index.final,]
test.labels <- cl[index.final]

heatmap.2(x = heatmap_test , scale = "row", Rowv =FALSE , Colv = FALSE, ColSideColors =c("blue","red","green")[test.labels], labRow = signature, labCol = NA, main = ' \n Testing Set \n ', col =hmcols, cexRow =0.40, trace ="none", key.title = "Color Key")



### Clearly as there are batch effects ######################
#### We should apply cross platform normalization ##########
## Make sure your R working directory is set to the folder containing the 


##################################################################################################################
##### Using the 840 Gene Signature to predict classes for the validation data set #################################
###################################################################################################################
Y_train <- t(train.pre[as.character(signature),])
Y_test  <- t(test.man[as.character(signature),])


#### Use Cross Platform Normalization ###########################
eset.train <- ExpressionSet(as.matrix(t(Y_train)))
eset.test <- ExpressionSet(as.matrix(t(Y_test)))

## As both columns matrices should contain the same amount of genes in the rows for binding using cbind
# cross-platform integration
list.eset <- list(eset.train,eset.test)
## Using Empirical Bayes to perform cross platform Normalization
merged.eset <- merge(list.eset, method ="COMBAT")


### Creating ready data matrices
Y.ready <- t(exprs((merged.eset[,sampleNames(eset.train)])))
Y.test.ready <- t(exprs((merged.eset[,sampleNames(eset.test)])))

## Just Checking ###########
data.combined <- rbind(Y.ready,Y.test.ready)
pc <- prcomp(data.combined)
pc.pred <- predict(pc,newdata = data.combined)
plot(pc.pred[,1], pc.pred[,2], pch = 19,col = lab)



## Create vector for class membership.  For these example data, we just 
## need 10 ones followed by 10 twos, then 10 threes, then 10 fours.
id <- pheno[,4]
levels(id) <- c(1,2,3,4)

knn.model <- knn(train = Y.ready, test = Y.test.ready, id, k = 4, l = 0, prob = FALSE, use.all = TRUE)


####### Read the Validation data classification annotation from other file ###########################################
index.final <- c(which(knn.model == 4), which(knn.model == 3), which(knn.model == 1), which(knn.model ==2))
heatmap_test <- Y.test.ready[index.final,]
test.labels <- knn.model[index.final]


### checking the pattern in training data #####################
hmcols<-colorRampPalette(c("green","black","red"))(256)
heatmap.2(x = t(heatmap_test), scale = "row",Rowv = FALSE ,Colv = FALSE, labRow = signature, labCol = NA, main = ' \n Testing Set \n ', col =hmcols, ,ColSideColors =c("blue","red","green","purple")[test.labels], cexRow =0.40, trace ="none", key.title = "Color Key")


