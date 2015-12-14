# 
# thetaOpt<-optim(c(1,1), lnGDP,data=dataVec)$par

# Maximum liklihood estimator of GPD. returns a vector c(xi,sigma)
mleGPD<-function(data,u,start=c(1,1)){
  if(is.matrix(data)){data<-as.vector(data)}
  y=data-u
  y=y[y>0]
  X<-optim(start,nlGPD,data=y,hessian=TRUE)
  return(list(par=X$par, covar=solve(X$hessian)))
}
# negative likelihood function for GPD (negative for optim)
nlGPD<-function(data,X){
  if(X[1]<(-X[2]/max(data))||X[2]<=0){return(-Inf)}
  return(length(data)*log(X[2])+(1+1/X[1])*sum(log(1+X[1]*data/X[2])))
}