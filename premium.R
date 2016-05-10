######### This R-file is used to run the PReMiuM package in R ###########################
######### It only needs the directory where the output will be stored
premium = function(dir){

  
setwd(dir)
data.premium <- list(NULL)
data.premium$inputData <- as.data.frame(cbind(Y,time))
data.premium$covNames <- paste0(rep("Variable",D),as.character(c(1:D)))
names(data.premium$inputData) <- c(data.premium$covNames,"outcome")
data.premium$xModel <- "Normal"
data.premium$yModel <- "Normal"
data.premium$ncovariates <- D
data.premium$outputData <- as.data.frame(Y.new)
names(data.premium$outputData) <- paste0(rep("Variable",D),as.character(c(1:D)))


####### Run the Profile Regression

runInfoObj <-profRegr(yModel=data.premium$yModel, xModel= data.premium$xModel, nSweeps= 200, nBurn= 100, data= data.premium$inputData, output= "outcome", nFilter= 5 ,covNames= data.premium$covNames,nClusInit= k,reportBurnIn=FALSE, predict = data.premium$outputData)
dissimObj<- calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)


######## Rand Indices ##############################
rand.recovery.premium <<- adjustedRandIndex(c.true,clusObj$clustering)
rand.predicted.premium <<-  adjustedRandIndex(c.true.new,clusObj$clusteringPred)

########## Making predictions #####################
riskProfileObj <- calcAvgRiskAndProfile(clusObj)

output_predictions <- calcPredictions(riskProfileObj)

cindex.predicted.premium <<-  as.numeric(survConcordance(Surv(exp(time.new),censoring.new) ~ exp(-output_predictions$predictedY))[1]) 

}



