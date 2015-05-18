posterioralpha = function(c, N, alpha, shape.alpha, rate.alpha) {
  
  ## Calculate the number of Active clusters
  ## Calculate the number of Active clusters
  
  g <- table(factor(c, levels = 1:K))
  active <- which(g!=0)
  Kplus <- length(active)
  
  
  f = function(x, N = N, Kp = Kplus){
    (Kp - 1.5)* log (x) + lgamma(x) - lgamma(N+x) - (0.5/(x))
    
  }
  
  fprima = function(x, N = N, Kp = Kplus){
    (Kp - 1.5)* (1/x) + digamma(x) - digamma(N+x) + (0.5/ (x)^2)
    
  }
  
  ars(1, f, fprima, x = alpha , m =1, N = N, Kp = Kplus ) 
  

}

