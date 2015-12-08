# Pareto dist
# MCMC prior xi (v=100), prior phi (v=10000) [ref onenote book page 174]

################
# y=X-u (u>=1) #
################
dataVec=as.vector(t(data))-1
library(MASS)
n<-10000
#normMu<-rep(.Machine$double.eps*100,2)
normMu<-rep(1,2)
coVar<-matrix(c(1,0,0,1),2,2)
#normMu<-c(-0.014,0.437)
#coVar<-matrix(c(1.032e-6,-3.552e-6,-3.552e-6,1.01e-4),2)
lamda<-2.38^2/length(normMu)
theta<-matrix(rep(0,2*n),2)
theta[,1]<-normMu
tau<--0.1*n/log(0.2)
j<-0
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
  gamma<-0.5*exp(-i/tau)
  #lamda<-lamda*exp(gamma*(r-a)) #where a is the raget acceptance rate
  coVar<-coVar+gamma*((theta[,i]-normMu)%*%t(theta[,i]-normMu)-coVar)
  normMu<-normMu+gamma*(theta[,i]-normMu)
}
X<-list(data=dataVec, estimators=theta, var=coVar, mu=normMu, burnin=NA, aRate=NA, MLE=NA )
class(X)<-'MCMC'
plot(density(X$estimators[1,1000:n]))

prob<-lnGDP(dataVec,theta)
which(prob==max(prob))
X$estimators[,7060]


# data is matrix/vector of data (X), threshold aka u, start is starting values, mu and var is mean and
# variance for the random walk starting distribution (g(x*|x)), gamma is the adaptive paramter (exp(-x/t)), 
# a is the optimal accaptance rate.
MCMCGDP<-function(data,threshold,start=NULL,mu=NULL,var=NULL,n=1000,gamma=NULL,a=NULL){
  require(MASS)
  data<-data-threshold
  dataLength<-length(data)
  if(is.matrix(data)){data=as.vector(t(data))}
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
  if(length(start[1,])>1){
    startBest<-c(0,0)
    lnBest<--Inf
    temp<-NA
    for(i in 1:length(start[1,])){
      temp<-mcmcGDP(data,0,start[,i],mu=mu,var=var,n=200,gamma=0.5,a=a)
      if(temp$lnMax>lnBest){
        lnBest<-temp$lnMax
        startBest<-temp$mle
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
    gamma<-exp(-(1:n)/(-0.1*n/log(0.05)))
  }else{
    gamma<-rep(gamma,n)
  }
  for(i in 2:n){
    Xtemp<-mvrnorm(n=1,theta[,(i-1)],lamda*coVar)
    R=lnRGDP(data,dataLength,Xtemp,theta[,(i-1)])
    if(uni[i-1]<R){
      theta[,i]<-Xtemp
    }else{
      theta[,i]<-theta[,(i-1)]
    }
    if(!is.null(a)){
      lamda<-lamda*exp(gamma[i]*(R-a)) #better way? 2 difference for?
    }
    var<-var+gamma[i]*((theta[,i]-mu)%*%t(theta[,i]-mu)-var)
    mu<-mu+gamma[i]*(theta[,i]-mu)
  }
  #####################################################################
  
  MCMC<-list(estimators=theta, threshold=threshold, var=var, mu=mu, burnin=NA, aRate=NA, MLE=NA, MLEest=NA)
  class(MCMC)<-'MCMC'
  return(MCMC)
}



# l(xi,si)=n*log(si)-(1+1/xi)*sum(log(1+xi*data/si)))
lnRGDP<-function(data,dataLength,Xtemp,X){
  if(Xtemp[1]< (-exp(Xtemp[2])/max(data))){return(-Inf)}
  lnGPNtemp<--dataLength*Xtemp[2]-(1+1/Xtemp[1])*sum(log(1+Xtemp[1]*data/exp(Xtemp[2])))
  lnpxitemp<-dnorm(Xtemp[1], mean = 0, sd = 100, log = TRUE)
  lnpphitemp<-dnorm(Xtemp[2], mean = 0, sd = 10000, log = TRUE)
  lnGPN<--dataLength*X[2]-(1+1/X[1])*sum(log(1+X[1]*data/exp(X[2])))
  lnpxi<-dnorm(X[1], mean = 0, sd = 100, log = TRUE)
  lnpphi<-dnorm(X[2], mean = 0, sd = 10000, log = TRUE)
  return(lnGPNtemp+lnpxitemp+lnpphitemp-lnGPN-lnpxi-lnpphi)  
}

lnGDP<-function(data,X){
  if(X[1]<(-exp(X[2])/max(data))){return(-Inf)}
  if(is.matrix(X)){
    temp<-rep(NA,length(X[1,]))
    for(i in 1:length(X[1,])){
      temp[i]<-sum(log(1+X[1,i]*data/exp(X[2,i])))
    }
    return( -length(data)*X[2,]-(1+1/X[1,])*temp)
  }else if(is.numeric(X)){
    return( -(-length(data)*X[2]-(1+1/X[1])*sum(log(1+X[1]*data/exp(X[2]))))) #minus since optim finds minimum
  }
}
thetaOpt<-optim(c(1,1), lnGDP,data=dataVec)$par
lines(1:200,(1+thetaOpt[1]*(1:200-1)/exp(thetaOpt[2]))^(-1/thetaOpt[1]),col=2)
lines(1:200,(1+theta[1,200]*(1:200-1)/exp(theta[2,200]))^(-1/theta[1,200]),col=3)


lines.mcmc<-function(data,X,xlim){
  # which estmate=max prob
  if(X[1]<0){
    maxvalue<- -exp(X[2])/X[1]-.Machine$double.eps
    if(is.null(xlim)){
      xlim<-c(min(data),maxvalue)
    }else if(xlim[2]>maxvalue){
      xlim[2]<-maxvalue
    }
  }
  if(is.null(xlim)){
    xlim<-c(min(data),3*(max(data)-min(data))+min(data))
  }
  x<-seq(xlim[1],xlim[2],length.out=500)
  print(length(x))
  y=(1+X[1]*x/exp(X[2]))^(-1/X[1])
  print(length(y))
  lines(x,y,type='l',col='blue')
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
