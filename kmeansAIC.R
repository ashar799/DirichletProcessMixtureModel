kmeansAIC = function(fit){
m.aic = ncol(fit$centers)
n.aic = length(fit$cluster)
k.aic = nrow(fit$centers)
d.aic = fit$tot.withinss
return(d.aic + 2*m.aic*k.aic)
}

