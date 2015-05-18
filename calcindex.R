calcindex = function(c, Time, time.predicted ) {
  numclust <- table(factor(c, levels = 1:K))
  activeclust <- which(numclust!=0)

  cis <- c(0)
  
  for (j in 1:length(activeclust)) {
    clust <- which(c==activeclust[j])
    ob.surv <- Surv(as.vector(time.predicted[clust]), Time[clust,2])
    km <- survfit(ob.surv~1)
    survest <- stepfun(km$time, c(1, km$surv))
    predicted.survival <- survest(time.predicted[clust])
    cis[j] <- concordance.index(x = predicted.survival, surv.time = Time[clust,1], surv.event= Time[clust,2])$c.index 
    }

  
  list('cindex' = cis)
}
