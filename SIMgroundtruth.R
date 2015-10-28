############ Ground Truth on TRAINING DATA ###################################
##############################################################
###########


### K-means + CoxPH
### K-means + AFT

### K-means + Penalized CoxPH
### K-means + Penalized AFT

### FlexMix +  CoxPH
### FlexMix +  AFT


SIMgroundtruth = function(){
  
smod <-  Surv(exp(time), censoring)
  
#############################################
########### K-means #########################
#############################################
#############################################
gr.km <- kmeans(Y, F, nstart =10)
gr.km.rand <- adjustedRandIndex(c.true,as.factor(gr.km$cluster))






########## CoxPH #############################
fit.cox.km <- coxph(smod ~ Y[,1:D] + strata(as.factor(gr.km$cluster)), data = as.data.frame(Y))
## C-Index
cindex.km.cox <- survConcordance(smod ~ predict(fit.cox.km))[1]
## Brier Score
fit.coxph <- survfit(fit.cox.km, newdata = as.data.frame(Y[,1:D]))
brier.km.cox <- sbrier(Surv(fit.coxph$time,fit.coxph$n.event), fit.coxph$surv)[[1]]
                     


######## AFT ###################################
fit.aft.km <- survreg(smod ~ Y[,1:D] + strata(as.factor(gr.km$cluster)) , dist="lognormal")
cindex.km.aft <- rcorr.cens(predict(fit.aft.km), smod)[1]
# predict.km.aft <- predict(fit.aft.km,newdata = as.data.frame(Y),type ="quantile",p=seq(.01,.99,by=.01))
# predict.km.aft <- predict(fit.aft.km,newdata = as.data.frame(Y),type ="response")



### Brier Score
brier.km.aft <- c(0)
for ( q in 1:F){
  ind <- which((gr.km$cluster) == q)
  time.tmp <- time[ind]
  censoring.tmp <- censoring[ind]
  Y.tmp <- Y[ind,]
  rownames(Y.tmp) <- as.character(c(1:nrow(Y.tmp)))
  smod.tmp <-  Surv(exp(time.tmp), censoring.tmp) 
  f1 <- psm(smod.tmp ~ Y.tmp[,1:D]  , dist="lognormal")
  S1 <- Survival(f1)
  L = length(ind)
  mat.tmp <- matrix(NA, nrow = L, ncol = L)
  for (j in 1:L){
    mat.tmp[,j] <- S1(exp(time.tmp[j]),f1$linear.predictors)
  }
  brier.km.aft[q] <- sbrier(smod.tmp,mat.tmp,exp(time.tmp))[1]
  
}




brier.km.pcox <- c(0)
cindex.km.pcox <- c(0)
cindex.km.pcox2 <- c(0)
brier.km.paft <- c(0)
cindex.km.paft <- c(0)


######## Penalized Cox PH ###########################################
for ( q in 1:F){
ind <- which((gr.km$cluster) == q)
time.tmp <- time[ind]
censoring.tmp <- censoring[ind]
Y.tmp <- Y[ind,]
coxreg <- list(0)
coxreg$x <- Y.tmp
coxreg$time <- exp(time.tmp)
coxreg$status <- censoring.tmp
path <- coxpath(data = coxreg)
f.reg <- predict(object = path, data = coxreg, s =5,type =  "coxph", mode = "lambda")
fit.coxregph <-  survfit(f.reg, newdata = as.data.frame(Y.tmp[,1:D]))
brier.km.pcox[q] <- sbrier(Surv(fit.coxregph$time,fit.coxregph$n.event), fit.coxregph$surv)[[1]]
cindex.km.pcox[q]  <-  survConcordance(Surv(coxreg$time,coxreg$status) ~ predict(f.reg))[1]
 ### see if we can use glmnet
 reg.pcox <- cv.glmnet(x = Y.tmp, y = Surv(coxreg$time, coxreg$status), family = "cox")
linear.pred <- predict(object =reg.pcox, newx = Y.tmp, s= "lambda.min")
cindex.km.pcox2[q] <- survConcordance(Surv(coxreg$time,coxreg$status) ~ linear.pred)[1]
}


######## Penalized AFT ######################################################

for ( q in 1:F){
ind <- which((gr.km$cluster) == q)
L= length(ind)

time.tmp <- time[ind]
censoring.tmp <- censoring[ind]
Y.tmp <- Y[ind,]
  
reg <- cv.glmnet(x = Y.tmp, y = time.tmp, family = "gaussian")
linear.pred <- predict(object =reg, newx = Y.tmp, s= "lambda.min")
coeff.pred <- coef(object =reg, newx = Y.tmp, s= "lambda.min")
rel.coeff <- coeff.pred[2:(D+1)] 
ind.rel <- which(rel.coeff !=0)

### Calculate the error variance as it is required to calcuate the Survival Probabilities for AFT model
predicted.tmp <- predict(object = reg, newx = Y.tmp, s = "lambda.min") 
sum.err <-c(0)
cindex.km.paft[q] <- survConcordance(Surv(time.tmp,censoring.tmp) ~ exp(-predicted.tmp))[1]

sum.err <- c(0)
for ( j in 1:L){
  sum.err[j] <- (time.tmp[j] - predicted.tmp[j])^2 
}
variance = sqrt(median(sum.err))

## This is just to get a survival function1 - pnorm((t.trans - lp)/exp(parms))
aft_survival = function (times = NULL, lp = NULL, parms = variance) 
{
  t.trans <- logb(times)
  names(t.trans) <- format(times)
  1 - pnorm((t.trans - lp)/(parms))
}

mat.reg.tmp <- matrix(NA, nrow = L, ncol = L)
for (j in 1:L){
  mat.reg.tmp[,j] <- aft_survival(exp(time.tmp[j]),predicted.tmp)
}
brier.km.paft[q] <- sbrier(Surv(time.tmp,censoring.tmp),mat.reg.tmp,exp(time.tmp))[1]
}



##### Save some Ground truth statistics

gr.km.rand.final <<- gr.km.rand
cindex.km.cox.final <<- as.numeric(cindex.km.cox)
cindex.km.aft.final <<-  as.numeric(cindex.km.aft)

brier.km.cox.final <<- as.numeric(brier.km.cox)
brier.km.aft.final <<- mean(brier.km.aft)

cindex.km.pcox.final <<- mean(unlist(cindex.km.pcox))
cindex.km.paft.final <<- mean(unlist(cindex.km.paft))

brier.km.pcox.final <<- mean(unlist(brier.km.pcox))
brier.km.paft.final <<- mean(unlist(brier.km.paft))



#################################################################################
##############################################################################
############### FlexMix #######################################################
################################################################################

gr.flx <- flexmix(time ~ Y, k =F)
gr.flx.rand <- adjustedRandIndex(c.true,clusters(gr.flx))



########## CoxPH #############################
fit.cox.flx <- coxph(smod ~ Y[,1:D] + strata(as.factor(clusters(gr.flx))), data = as.data.frame(Y))
## C-Index
cindex.flx.cox <- survConcordance(smod ~ predict(fit.cox.flx))[1]
## Brier Score
fit.coxph <- survfit(fit.cox.flx, newdata = as.data.frame(Y[,1:D]))
brier.flx.cox <- sbrier(Surv(fit.coxph$time,fit.coxph$n.event), fit.coxph$surv)[[1]]


gr.flx.rand.final <<- gr.flx.rand
cindex.flx.cox.final <<- as.numeric(cindex.flx.cox)
brier.flx.cox.final <<-  brier.flx.cox






}

