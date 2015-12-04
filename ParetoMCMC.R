# Pareto dist
# MCMC prior xi (v=100), prior phi (v=10000) [ref onenote book page 174]
dataVec=as.vector(t(data))
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
  gamma<-1/sqrt(i)
  #lamda<-lamda*exp(gamma*(r-a)) #where a is the raget acceptance rate
  coVar<-coVar+gamma*((theta[,i]-normMu)%*%t(theta[,i]-normMu)-coVar)
  normMu<-normMu+gamma*(theta[,i]-normMu)
}




mcmcGDP<-function(data,start=NULL,mu=NULL,var=NULL,n=1000,lamda=NULL,gamma=NULL){
  # Fix starting values
  if(is.null(start)){
    start<-matrix(rep(NA,10),2)
    start[,1]<-rep(.Machine$double.eps*100,2)
    start[,2]<-c(1,log(1))
    start[,3]<-c(1,log(5))
    start[,4]<-c(-1,log(1))
    start[,5]<-c(-1,log(5))
  }
  if(all(start==0)){
    start<-rep(.Machine$double.eps*100,2)
  }
  
  # Fix mu and var for MCMC
  if(is.null(mu)){mu=c(0,0)}
  if(is.null(var)){var<-matrix(c(1,0,0,1),2,2)}
  
  # Choose starting value
  if(length(start)>2){
    startBest<-rep(.Machine$double.eps*100,2)
    lnBest<-lnGDP(data,startBest)
    temp<-NA
    for(i in 1:length(start[1,])){
      temp<-mcmcGDP(data,start[,i],mu=mu,var=var,n=200,gamma=0.5)
      if(temp$lnMax>lnBest){
        lnBest<-temp$lnMax
        startBest<-temp$mle
      }
    }
    start<-startBest
  }
  
  #MCMC code
  #
  #
  #####################################################################
  
  
  MCMC<-list()
  class(MCMC)<-'MCMC'
  return(MCMC)
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

lnGDP<-function(data,X){
  if(is.matrix(X)){
    temp<-rep(NA,length(X[1,]))
    for(i in 1:length(X[1,])){
      temp[i]<-sum(log(1+X[1,i]*data/exp(X[2,i])))
    }
    return(-length(data)*X[2,]-(1+1/X[1,])*temp)
  }else if(is.numeric(X)){
    return(-length(data)*X[2]-(1+1/X[1])*sum(log(1+X[1]*data/exp(X[2]))))
  }
}












# c++ :D :D !!!
# http://adv-r.had.co.nz/Rcpp.html
library(Rcpp)

cppFunction('double sumC(NumericVector x) {
  int n = x.size();
            double total = 0;
            for(int i = 0; i < n; ++i) {
            total += x[i];
            }
            return total;
            }')


sumR <- function(x) {
  total <- 0
  for (i in seq_along(x)) {
    total <- total + x[i]
  }
  total
}

x=runif(100000000)
tm1 <- system.time(
  {
    sum(x)
  })

tm2 <- system.time(
  {
    sumC(x)
  })
