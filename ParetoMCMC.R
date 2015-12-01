# Pareto dist
library(MASS)
n<-1000
normMu<-c(0,0)
coVar<-matrix(c(1,0,0,1),2,2)
lamda<-2.38^2/length(normMu)
theta<-matrix(rep(0,2*n),2)
for(i in 1:n){
  Xtemp<-mvrnorm(n=1,normMu,lamda*coVar)
  
  
}

RGDP<-function(Xstar,X){
  
}