# Pareto dist
# MCMC prior xi (v=100), prior phi (v=10000) [ref onenote book page 174]
dataVec=as.vector(data)
library(MASS)
n<-10000
normMu<-c(-0.014,0.437)#normMu<-rep(.Machine$double.eps*100,2)
coVar<-matrix(c(1.032e-6,-3.552e-6,-3.552e-6,1.01e-4),2)#coVar<-matrix(c(1,0,0,1),2,2)
lamda<-2.38^2/length(normMu)
theta<-matrix(rep(0,2*n),2)
theta[,1]<-normMu
gamma<-1
j=0
for(i in 2:n){
  Xtemp<-mvrnorm(n=1,theta[,(i-1)],lamda*coVar)
  u=log(runif(1))
  R=lnRGDP(dataVec,Xtemp,theta[,(i-1)])
  if(u<R){
    #cat(u,"<",R,"\n")
    theta[,i]<-Xtemp
    j=j+1
  }else{
    theta[,i]<-theta[,(i-1)]
  }
  gamma<-1/i
  #lamda<-lamda*exp(gamma*(r-a)) #where a is the raget acceptance rate
  coVar<-coVar+gamma*((theta[,i]-normMu)%*%t(theta[,i]-normMu)-coVar)
  normMu<-normMu+gamma*(theta[,i]-normMu)
}



# l(xi,si)=n*log(si)-(1+1/xi)*sum(log(1+xi*data/si)))


lnRGDP<-function(data,Xtemp,X){
if(Xtemp[1]< (-exp(Xtemp[2])/max(data))){return(-Inf)}
lnGPNtemp<--length(data)*Xtemp[2]-(1+1/Xtemp[1])*sum(log(1+Xtemp[1]*data/exp(Xtemp[2])))
lnpxitemp<-dnorm(Xtemp[1], mean = 0, sd = 100, log = TRUE)
lnpphitemp<-dnorm(Xtemp[2], mean = 0, sd = 10000, log = TRUE)
lnGPN<--length(data)*X[2]-(1+1/X[1])*sum(log(1+X[1]*data/exp(X[2])))
lnpxi<-dnorm(X[1], mean = 0, sd = 100, log = TRUE)
lnpphi<-dnorm(X[2], mean = 0, sd = 10000, log = TRUE)
return(lnGPNtemp+lnpxitemp+lnpphitemp-lnGPN-lnpxi-lnpphi)  
}



