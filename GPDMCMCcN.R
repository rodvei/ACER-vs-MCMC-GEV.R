test<-function(data){
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
  lnGPNtemp=0
  lnpxitemp=0
  lnpphitemp=0
  lnGPN=0
  lnpxi=0
  lnpphi=0
  
  for(i in 2:n){
    Xtemp<-mvrnorm(n=1,theta[,(i-1)],lamda*coVar)
    u=log(runif(1))
    
    
    
    if(Xtemp[1]< (-exp(Xtemp[2])/max(data))){R=-Inf
    }else{
      lnGPNtemp<--length(data)*Xtemp[2]-(1+1/Xtemp[1])*sum(log(1+Xtemp[1]*data/exp(Xtemp[2])))
      lnpxitemp<-dnorm(Xtemp[1], mean = 0, sd = 100, log = TRUE)
      lnpphitemp<-dnorm(Xtemp[2], mean = 0, sd = 10000, log = TRUE)
      lnGPN<--length(data)*theta[2,(i-1)]-(1+1/theta[1,(i-1)])*sum(log(1+theta[1,(i-1)]*data/exp(theta[2,(i-1)])))
      lnpxi<-dnorm(theta[1,(i-1)], mean = 0, sd = 100, log = TRUE)
      lnpphi<-dnorm(theta[2,(i-1)], mean = 0, sd = 10000, log = TRUE)
      R=lnGPNtemp+lnpxitemp+lnpphitemp-lnGPN-lnpxi-lnpphi
    }
    
    
    
    
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
}