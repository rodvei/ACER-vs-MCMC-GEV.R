#library(Rcpp)
#sourceCpp("Cpp/MCMC_GPD_CPP.cpp")

# data is matrix/vector of data (X), u aka thershold, start is starting values ([type,number])
# type is xi and sigma (=exp(phi)), 
# mu and var is mean and variance for the random walk starting distribution (g(x*|x)), 
# gamma is the adaptive paramter (exp(-x/t)), 
# a is the optimal accaptance rate.
MCMCGPD<-function(data,u,start=NULL,mu=NULL,var=NULL,n=1000,gamma=NULL,a=NULL){
  require(MASS)
  
  if(is.matrix(data)){data=as.vector(t(data))}
  y<-data-u
  y=y[y>0]
  dataLength<-length(data)
  
  # Fix starting values
  if(is.null(start)){
    start<-matrix(rep(NA,10),2)
    start[,1]<-rep(.Machine$double.eps*100,2)
    start[,2]<-c(1.5,0.3)
    start[,3]<-c(0.5,3)
    start[,4]<-c(-0.1,0.3)
    start[,5]<-c(-0.5,3)
  }else if(all(start==0)){
    start<-rep(.Machine$double.eps*100,2)
  }else{
    if(length(start)==2){
      start<-c(start[1],log(start[2]))
    }else{
      start<-rbind(start[1,],log(start[2,]))
    }
  }
  
  
  # Fix mu and var for MCMC
  if(is.null(mu)){mu=c(0,0)}
  if(is.null(var)){var<-matrix(c(1,0,0,1),2,2)}
  
  # Choose starting value
  if(length(start)>2){
    startBest<-c(0,0)
    lnBest<--Inf
    temp<-NA
    for(i in 1:length(start[1,])){
      if(start[1,i]< (-start[2,i]/max(data))){
        start[1,i]<- 0.9*(-start[2,i]/max(data))
      }
      temp<-MCMCGPD(data,u,start[,i],mu=mu,var=var,n=500,gamma=0.1,a=a)
      lnTemp<-lnGPD(y,c(temp$theta[1,500],temp$theta[2,500]))
      if(lnTemp>lnBest){
        lnBest<-lnTemp
        startBest<-c(temp$theta[1,500],log(temp$theta[2,500]))
      }
    }
    start<-startBest
  }
  
  ##############################MCMC###################################
  lamda<-2.38^2/2
  theta<-matrix(rep(NA,2*n),2)
  theta[,1]<-start
  uni=log(runif(n-1))
  if(is.null(gamma)){
    gamma<-0.5*exp(-(1:n)*log(10)/(0.1*n))
  }else{
    gamma<-rep(gamma,n)
  }
  for(i in 2:n){
    Xtemp<-mvrnorm(n=1,theta[,(i-1)],lamda*var)
    R=lnRGPD(y,Xtemp,theta[,(i-1)])
    if(uni[i-1]<R){
      theta[,i]<-Xtemp
    }else{
      theta[,i]<-theta[,(i-1)]
    }
    if(!is.null(a)){
      if(R>=0){R=1
      }else{R=exp(R)}
      lamda<-lamda*exp(gamma[i]*(R-a)) #better way? 2 difference for?
    }
    var<-var+gamma[i]*((theta[,i]-mu)%*%t(theta[,i]-mu)-var)
    mu<-mu+gamma[i]*(theta[,i]-mu)
  }
  #####################################################################
  theta
  MCMC<-list(theta=rbind(theta[1,],exp(theta[2,])), u=u, var=var, mu=mu, burnin=NA, aRate=NA, MLE=NA, MLEest=NA)
  class(MCMC)<-'MCMC'
  return(MCMC)
}

# l(xi,si)=n*phi-(1+1/xi)*sum(log(1+xi*data/exp(phi)))) # phi=log(si)
lnRGPD<-function(data,Xtemp,X){
  dataLength<-length(data)
  if(Xtemp[1]< (-exp(Xtemp[2])/max(data))){return(-Inf)}
  lnGPNtemp<--dataLength*Xtemp[2]-(1+1/Xtemp[1])*sum(log(1+Xtemp[1]*data/exp(Xtemp[2])))
  lnpxitemp<-dnorm(Xtemp[1], mean = 0, sd = 100, log = TRUE)
  lnpphitemp<-dnorm(Xtemp[2], mean = 0, sd = 10000, log = TRUE)
  lnGPN<--dataLength*X[2]-(1+1/X[1])*sum(log(1+X[1]*data/exp(X[2])))
  lnpxi<-dnorm(X[1], mean = 0, sd = 100, log = TRUE)
  lnpphi<-dnorm(X[2], mean = 0, sd = 10000, log = TRUE)
  return(lnGPNtemp+lnpxitemp+lnpphitemp-lnGPN-lnpxi-lnpphi)  
}

# X is parameters (aka theta) matrix of thetas or single theta(c(xi,sigma)), data is a vector.
lnGPD<-function(data,X){
  if(X[1]<(-X[2]/max(data))){return(-Inf)}
  if(is.matrix(X)){
    temp<-rep(NA,length(X[1,]))
    for(i in 1:length(X[1,])){
      temp[i]<-sum(log(1+X[1,i]*data/X[2,i]))
    }
    return( -length(data)*log(X[2,])-(1+1/X[1,])*temp)
  }else if(is.numeric(X)){
    return(-length(data)*log(X[2])-(1+1/X[1])*sum(log(1+X[1]*data/X[2])))
  }
}

