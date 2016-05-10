############ Ground Truth on TRAINING DATA ###################################
########### Some other methods which take care of noise ######################


########## Mixture of Factor Analyzers ##################
########## Sparse Hierarchical Clustering #################
### Sparse k-means #######################################


iSIMgroundtruth = function(){
  
  ############## Testing Mixture of Factor Analyzers ##########################3
  
  mcfa.fit<- mcfa(Y, g= k, q=2, itmax=250, nkmeans=5, nrandom=5, tol=1.e-3)
  
  ########## Seeing if the PCA plot with the corresponding features with releevant features makes sense
  randindexMCFA <<- adjustedRandIndex(mcfa.fit$clust, c.true)
  
  
  
  #############################################
  ########### sparse K-means #########################
  #############################################
  #############################################
  km.perm <- KMeansSparseCluster.permute(x = Y, K= k ,wbounds=c(1.5,2:6),nperms=5)
  km.out <- KMeansSparseCluster(x = Y, K=k,wbounds=km.perm$bestw)
  randindexSKM <<- adjustedRandIndex(km.out[[1]]$Cs, c.true)
  
  ###################################################
  ########### sparse hierarchical clustering #########################
  #############################################
  #############################################
  perm.out <- HierarchicalSparseCluster.permute(x = Y, wbounds=c(1.5,2:6), nperms=5)
  # Perform sparse hierarchical clustering
  sparsehc <- HierarchicalSparseCluster(dists=perm.out$dists, wbound=perm.out$bestw, method="complete")
  randindexSHC <<- adjustedRandIndex(cutree(sparsehc$hc, k = k), c.true)
  
 
  
  
  
}
