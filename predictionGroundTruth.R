########## This file compares prediction on test cases with different methods
####### So we use k-means + SVM to predict labels on new data points ###########

predictionGroundTruth = function(){
  

############ Predicting New Class Labels using SVM #################################  
# gr.km <- kmeans(Y, F, nstart =10)
# label.train <- gr.km$cluster
# svms <- sapply(2^(-10:14), function(cost) cross(ksvm(Y, factor(label.train), C=cost, kernel="vanilladot", kpar=list(), cross=5)))
# mysvm <- ksvm(Y, factor(label.train), C=2^(-10:14)[which.min(svms)], kernel="vanilladot", kpar=list(), cross=10) # accuracy ~97%
# pred.svm <-  predict(mysvm, Y.new)
# predRandIndex.svm <<- adjustedRandIndex(c.true.new, pred.svm)


############ Predicting the New Class Labels using kNN #############################
gr.km <- kmeans(Y, F, nstart =10)
label.train <- gr.km$cluster
knear <- knn(train = Y, test = Y.new, cl = label.train, k = F)
predRandIndex.knear <<- adjustedRandIndex(c.true.new, knear)

######## Predicting New C-Indices based on a Penalized Cox or AFT model#################### 
######## Penalized Cox PH ###########################################

linear.pred.cox <- c(0)
### see if we can use glmnet
reg.pcox <- cv.glmnet(x = Y, y = Surv(exp(time), censoring), family = "cox")
linear.pred.cox <- predict(object =reg.pcox, newx = Y.new, s= "lambda.min")
smod <-  Surv(exp(time.new), censoring.new)
predCIndex.cox <<- as.numeric(survConcordance(smod ~ linear.pred.cox)[1])

##### Penalized AFT Model #############################################
linear.pred.paft <- c(0)
### see if we can use glmnet
reg.paft <- cv.glmnet(x = Y, y = time, family = "gaussian")
linear.pred.paft <- predict(object = reg.paft, newx = Y.new, s= "lambda.min")
smod <-  Surv(exp(time.new), censoring.new)
predCIndex.aft <<- as.numeric(survConcordance(smod ~ exp(-linear.pred.paft))[1])

  
###### K-Means + KNN 

gr.km <- kmeans(Y, F, nstart =10)
label.train <- gr.km$cluster
label.test <- knn(train = Y, test = Y.new, cl = label.train, k = F)


#### penAFT ###################################################################

linear.kkpaft.prediction <- c(0)
for ( q in 1:F){
  ind <- which(label.train == q)
  ind.new <- which(label.test == q)
  reg.aft <- cv.glmnet(x = Y[ind,], y = time[ind], family = "gaussian")
  linear.kkpaft.prediction[ind.new] <- predict(object =reg.aft, newx = Y.new[ind.new,], s= "lambda.min")
}
predCIndex.kkpaft <<- as.numeric(survConcordance(smod ~ exp(-linear.kkpaft.prediction))[1])


###### penCox ###################################################################

linear.kkpcox.prediction <- c(0)
for ( q in 1:F){
  ind <- which(label.train == q)
  ind.new <- which(label.test == q)
  reg.aft <- cv.glmnet(x = Y[ind,], y = Surv(exp(time[ind]),censoring[ind]), family = "cox")
  linear.kkpcox.prediction[ind.new] <- predict(object =reg.aft, newx = Y.new[ind.new,], s= "lambda.min")
}
predCIndex.kkpcox <<- as.numeric(survConcordance(smod ~ linear.kkpcox.prediction)[1])


#### PenFLXMIX  to make predictions on Clusters ##################################################################
data <- data.frame(time, Y)
data.new <- data.frame(time.new, Y.new)
## The cross validation folds for choosing lambda
fo <- sample(rep(seq(10), length = nrow(data)))
gr.flx <- flexmix(time ~ ., data = data, k = F, cluster = gr.km$cluster, model = FLXMRglmnet(foldid = fo, adaptive= FALSE), control = list(iter.max = 500))
cl.flx <- clusters(gr.flx, newdata = data.new)
rand.flx.prediction <<-  adjustedRandIndex(c.true.new,as.factor(cl.flx))


  
  
  
}