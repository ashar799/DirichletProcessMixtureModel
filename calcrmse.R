calcrmse = function(time.real,time.predicted ) {
  error.time <- c(0)
  error.time <- time.real - time.predicted
  rmse.time <- sqrt((1/N)*sum(error.time^2))
  nmse <- rmse.time/(max(time.real)- min(time.real))
  
  list('rmse' = nmse)
}
