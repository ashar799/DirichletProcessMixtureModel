## Let's See if the C-index with the correct co-efficinets are good
c.initial <-0

cumalative.survival <- c(0)
for ( i in 1:N){  
  cumalative.survival[i] <- 1 - pnorm(time.real[i], mean = time.cluster$Mu[c.true[i]] + t(beta.list[[c.true[i]]]) %*% Ytemp[i,1:rel.D], sd = sqrt(time.cluster$S[c.true[i]]))
}

risk.survival <- 1 - cumalative.survival

for ( i in 1:F){
  clust <- which(c.true == activeclass[i])
  c.initial[i] <- concordance.index(x = risk.survival[clust],  surv.time = time.real[clust], surv.event= c(rep(1, as.integer(N * p.dist[i]))))$c.index 
}

c.initial2 <- c(0)
for ( i in 1:F){
  clust <- which(c.true == activeclass[i])
  c.initial2[i] <- survConcordance(Surv(time.real[clust], censoring[clust]) ~ risk.survival[clust]) 
}

c.initial.kp <- c(0)

for (j in 1:length(activeclass)) {
  clust <- which(c.true ==activeclass[j])
  ob.surv <- Surv(as.vector(time.predicted[clust]), censoring[clust])
  km <- survfit(ob.surv~1)
  survest <- stepfun(km$time, c(1, km$surv))
  predicted.survival <- survest(time.predicted[clust])
  cis[j] <- concordance.index(x = predicted.survival, surv.time = time.real[clust], surv.event= censoring[clust])$c.index 
}
