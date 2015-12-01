# Pareto dist
library(MASS)
n<-1000
normMu<-c(0,0)
coVar<-matrix(c(1,0,0,1),2,2)
lamda<-2.38^2/length(normMu)
theta<-matrix(rep(0,2*n),2)
for(i in 1:n){
  Xtemp<-mvrnorm(n=1,normMu,lamda*coVar)
  if(u<RGDP(data,Xtemp,theta[,(i-1)])){
    theta[,i]<-Xtemp
  }else{
    theta[,i]<-theta[,(i-1)]
  }
  normMu<-XXX
  covar<-XXX
  
}


lnGPD->function(data,si,xi){
  return(-length(data)*log(si)-(1+1/xi)*sum(log(1+xi*data/si)))
}

RGDP<-function(data,Xtemp,X){
lnGPNtemp<--length(data)*Xtemp[2]-(1+1/Xtemp[1])*sum(log(1+Xtemp[1]*data/exp(Xtemp[2])))
lnGPN<--length(data)*X-(1+1/X[1])*sum(log(1+X[1]*data/exp(X[2])))

return(exp( ))  
}