##### This cript finds out about the c-index in the case of a AFT model 
### Let's generate Molecular Data 5 dimensional

X.aft <- mvrnorm(n = 100 , mu = c(0.8,-0.2,1.5,0.7,-1.2), Sigma = diag(rep(1,5)))
X.sc.aft <- scale(X.aft, center = TRUE, scale = TRUE)
beta0.aft <- 5
betahat.aft <- c(2,-2,3,-3,5)



time.pur<- as.vector(t(betahat.aft) %*% t(X.sc.aft))
sigma.aft <- 2
time.noise <- rnorm(100, mean = beta0.aft, sd = sigma.aft )

time.aft <-  as.vector(time.pur + time.noise)

time.act <- exp(time.aft)
status <- rep(1,100)
surv.aft <- Surv(time.act, status)
#### Calculate the hazazrd rate function

time.mu <- time.pur + beta0.aft

#### Fitting A SurvReg Model
fit <- survreg(surv.aft ~ X.sc.aft[,1:5], dist="lognormal")
invers.predict <- (predict(fit))^-1
## Getting C-indices
library(survival)
survConcordance(surv.aft ~ invers.predict)
library(Hmisc)
rcorr.cens(predict(fit), surv.aft)

### Fitting Co-Proportional Hazard's Model
fit.coxph <- coxph(surv.aft ~ X.sc.aft[,1:5] )
survConcordance(surv.aft ~ predict(fit.coxph))
rcorr.cens(exp(-predict(fit.coxph)), surv.aft)

### Fitting our own Model
rcorr.cens(time.mu, surv.aft)
