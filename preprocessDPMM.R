preprocessDPMM = function(){

########### Mixture of Factor Analyzers ######################################################   
# mcfa.fit<- mcfa(Y.dat, g=3, q=2, itmax=250, nkmeans=5, nrandom=5, tol=1.e-3)
# adjustedRandIndex(macfa.fit$clust.c.true)
#   
######### sparse NMF + Kmeans ####################################################################
### Making the Matrix poistive
# Y.temp <- Y.dat - matrix(data = min(Y.dat), nrow = N, ncol =D)
# res <- nmf(t(Y.temp), 2, method = 'snmf/l',nrun = 20) 
# low.var <- which.min(apply(coef(res),1,var))  
# good.basis <- c(1:2)[- low.var]
# basis.good <-  basis(res)[,good.basis]
# ind.rel <<-  which(basis.good > 0)
# km <- kmeans(Y.temp[,ind.rel],centers = F,nstart =5)
# adjustedRandIndex(km$cluster,c.true)
# Y <<- Y.dat[,ind.rel]
# D <<- ncol(Y)

############## Sparse COX LASSO model This needs to hand tweaked   ############################





####### Sparse clustering ############################################
# per <- KMeansSparseCluster.permute(Y.dat,K=F,wbounds=seq(1.1,1.5,len=100),nperms=2)
# km.perm <- KMeansSparseCluster(x = Y.dat, K= F, wbounds = per$bestw, nstart = 20, silent = FALSE, maxiter=6, centers=NULL)
# weight <- unlist(km.perm)[1:D]
# clustering <-unlist(km.perm)[(D+1):(D+N)] 
#   
# Y <<- Y.dat[,which(weight!=0)]
# Y.new <<- Y.new.dat[,which(weight!=0)]
# D <<- ncol(Y)  
#   


}
