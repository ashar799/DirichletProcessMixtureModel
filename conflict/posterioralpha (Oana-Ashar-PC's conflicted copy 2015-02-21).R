posterioralpha = function(c, N, alpha) {
  
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

## Just to check if the the function is log concave
# y = c(rep(0,5000))
# i = 0.001
# for ( j in 1:5000) {
#  y[j] <- f(i, N , Kplus)
#   i = i +0.001
# }
# 
# xval <- seq(0.001, 5, by = 0.001,)
# plot(xval, y)

## Test
# f1<-function(x,shape,scale=1){(shape-1)*log(x)-x/scale}
# f1prima<-function(x,shape,scale=1) {(shape-1)/x-1/scale}
# mysample1<-ars(20,f1,f1prima,x=4.5,m=1,lb=TRUE,xlb=0,shape=2,scale=0.5)
