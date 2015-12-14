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
