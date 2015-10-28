## Gibb's sampling 
burninDPMM = function(){
  
source('flexposteriorchineseAFT.R')
source('posteriorGMMparametrs.R')
source('posteriortimeparameters.R')
source('updatetime.R')
source('priordraw.R')
source('likelihood.R')

iter.burnin = iter.burnin
cognate <- NA
param <- NA
paramtime <- NA
loglike<- rep(0, iter)  
timeparam <- NA
time.predicted <- c(0)
cindex <- c(0)

init.likli <- loglikelihood(c,Y,mu,S,alpha,That, beta0, betahat, sigma2, lambda2, tau2, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.dat) 
rmse <- c(0)
randy <- c(0)
likli <- c(0)


o =1
#################### BURNIN PHASE ###################################################
o.iter = o
print("BURNIN...PHASE")

pb <- txtProgressBar(min = o.iter, max = iter.burnin , style = 3)
for (o in o.iter:iter.burnin) {
  
  
  ################## PARAMETERS OF THE DP Mixture Model ######################################################
  ## Updating the parameters based on the observations 
  param <- posteriorGMMparametrs(c,Y,mu,S, alpha,K, epsilon, W, beta, ro,N,D )
  mu <- param$mean
  S <- param$precision
  paramtime <- posteriortimeparameters(c, That, lambda2,tau2,sigma2,beta0, betahat, Y, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.data)
  beta0 <- paramtime$beta0
  betahat <- paramtime$betahat
  sigma2 <- paramtime$sigma2
  lambda2 <- paramtime$lambda2
  tau2 <- paramtime$tau2
  
  ########################## THE HYPERPARAMETERS OF THE GMM #################################  
  source('posteriorhyperPLUS.R')  
  #  Updating the hyper paramters
  hypercognate <- posteriorhyperPLUS (c, Y, mu, S, epsilon, W, beta, ro )
  epsilon <- hypercognate$epsilon
  tmpW <- hypercognate$W
  W <- matrix(as.matrix(tmpW),nrow = D, ncol =D)
  ro <- hypercognate$ro
  
  source('posteriorbeta.R')
  if( o%%10 == 0){
    res <- try(posteriorbeta(c, beta, D, S, W))
    if (class(res) == "try-error"){
      beta = beta
    } else{
      beta <- posteriorbeta(c, beta, D, S, W)
      
    }
  } 
  
  ################# INDICATOR VARIABLE ##################################################################
  ## Updating the indicator variables and the parameters
  source('posteriorchineseAFT.R')
  cognate <- posteriorchineseAFT(c,Y,mu,S,alpha,That, beta0, betahat, sigma2, lambda2, tau2, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.dat)
  c <- cognate$indicator
  mu <- cognate$mean
  S <- cognate$precision
  beta0 <- cognate$beta0
  betahat <- cognate$betahat
  sigma2 <- cognate$sigma2
  lambda2 <- cognate$lambda2
  tau2 <- cognate$tau2
  
  ########################### The Concentration Parameter #################################################################
  
  
  source('posterioralpha.R') 
  # Updating the concentration parameter
  alpha <- posterioralpha(c, N, alpha, shape.alpha, rate.alpha)
  
  
  
  
  
  ######################## The Censored Times ###########################################################
  source('updatetime.R')
  # Updating the Time Variable
  ti <- NA
  ti <- updatetime(c, Y, Time,That, beta0, betahat, sigma2)
  That <- ti$time
  
  
  ##################### Print SOME Statistics #####################################################
  randy[o] <- adjustedRandIndex(c.true,as.factor(c))
#   print(randy[o])
  rmse[o] <- calcrmse(time.real,That)$rmse
#   print(rmse[o])
  likli[o] <- loglikelihood(c,Y,mu,S,alpha,That, beta0, betahat, sigma2, lambda2, tau2, K, epsilon, W, beta, ro,D, r, si, Time,N, sig2.dat)
#   print(likli[o])
#   print(o/iter.burnin)
  
  
  ##### Print the status bar
  Sys.sleep(0.1)
  setTxtProgressBar(pb, o)
} 

assign("alpha", alpha, envir = .GlobalEnv)
assign("ro", ro, envir = .GlobalEnv)
assign("That", That, envir = .GlobalEnv)
assign("c", c, envir = .GlobalEnv)
assign("epsilon", epsilon, envir = .GlobalEnv)
assign("W", W, envir = .GlobalEnv)
assign("mu", mu, envir = .GlobalEnv)
assign("S", S, envir = .GlobalEnv)
assign("beta0", beta0, envir = .GlobalEnv)
assign("betahat", betahat, envir = .GlobalEnv)
assign("sigma2", sigma2, envir = .GlobalEnv)
assign("lambda2", lambda2, envir = .GlobalEnv)
assign("tau2", tau2, envir = .GlobalEnv)
assign("randy.burnin", randy, envir = .GlobalEnv)
assign("rmse.burnin", rmse, envir = .GlobalEnv)
assign("likli.burnin", likli, envir = .GlobalEnv)





}

