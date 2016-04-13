# PARETO DISTRIBUTION(B)
library(ACER)
itL<-5

day<-365
year<-25

pYearly<-1/365
pDecade<-1/(10*365)
pCentury<-1/(100*365)
pMillennium<-1/(1000*365)
p=c(pYearly,pDecade,pCentury,pMillennium)

acerList<-list(BestVaR=rep(0,4),inCI=rep(0,4),nBestCIv=rep(0,4),proCI=matrix(NA,4,itL),nBestCIp=rep(0,4))
mcmcList<-acerList

realVaR<-NA
acerVaR<-NA
acerVarCI<-c(NA,NA)
mcmcVaR<-NA
mcmcVaRCI<-c(NA,NA)
best<-NA

for(i in 1:itL){
  print(i)
  u<-runif(day*year)
  beta<-runif(1,3,5)
  data<-matrix((1/(u))^(1/beta),year,day,byrow=TRUE)
  tresh<-sort(as.vector(data))[round(length(data)/3)]
  acerMod<-ACER(data,k=1,eta1=tresh,CI=0.95,stationary=TRUE,method='general', check.weibull = FALSE)
  mcmcMod<-mcmc.gpd(as.vector(data),u=tresh,n=15500)
  
  for(j in 1:4){
    # VaR RELATED
    realVaR<-1/p[j]^(1/beta)
    acerVaR<-reverseGeneral(p[j], acerMod$coef[1:5])
    acerVarCI<-c(reverseGeneral(p[j], acerMod$lowerCIcoef[1:5]),reverseGeneral(p[j], acerMod$upperCIcoef[1:5]))
    
    VaRdist<-reversePOTmcmc(p[j],mcmcMod)
    mcmcVaR<-mean(VaRdist)
    mcmcVaRCI<-distCI(VaRdist,alpha=0.95)

    if(abs(acerVaR-realVaR)<abs(mcmcVaR-realVaR)){
      acerList$BestVaR[j]<-acerList$BestVaR[j]+1
      best<-'a'
    }else if(abs(acerVaR-realVaR)>abs(mcmcVaR-realVaR)){
      mcmcList$BestVaR[j]<-mcmcList$BestVaR[j]+1
      best<-'m'
    }else{print('ERROR#1')}
    
    if((acerVarCI[1]<=realVaR)&&(realVaR<=acerVarCI[2])){
      acerList$inCI[j]<-acerList$inCI[j]+1
    }
    if((mcmcVaRCI[1]<=realVaR)&&(realVaR<=mcmcVaRCI[2])){
      mcmcList$inCI[j]<-mcmcVaRCI$inCI[j]+1
    }
    if(best=='a'){
      acerList$nBestCIv[j]<-acerList$nBestCIv[j]+(diff(acerVarCI)<diff(mcmcVaRCI))*1
    }else if(best=='m'){
      mcmcList$nBestCIv[j]<-mcmcList$nBestCIv[j]+(diff(acerVarCI)>diff(mcmcVaRCI))*1
    }else{print('ERROR#2')}
    
    # PRED FUTURE

  }

  
  
}






###################################################
############## Pareto functions####################
###################################################
# pareto density plot
reversePareto<-function(u,alpha=3){
  return(1/u^(1/alpha))
}
paretoProb<-function(z,prob=1e-6,alpha=3){
  return(1/prob*alpha*z^(-alpha-1)*(1-z^(-alpha))^(1/prob-1))
}
# how much of the furute prediction interval contains the real
futureCIcontain<-function(prob=1e-6,alpha=3,PI){
  return((1-PI[2]^(-alpha))^(1/prob)-(1-PI[1]^(-alpha))^(1/prob))
}


# From ACER package
reverseGeneral<-function(eps,est){
  return(est[2]+(1/(est[1]*est[5])*((eps/est[4])^(-est[5])-1))^(1/est[3]))
}
General<-function(eta,est){
  if(est[5]<0){eta[eta>(est[2]-1/(est[1]*est[5]))]=NA}
  return(est[4]*(1+est[1]*est[5]*(eta-est[2])^est[3])^(-1/est[5]))
}
# z_p pot MCMC
reversePOTmcmc<-function(prob,MCMC){
  return(MCMC$theta[2,]/MCMC$theta[1,]*((prob/MCMC$theta[3,])^(-MCMC$theta[1,])-1)+MCMC$u)
}
POTmcmc<-function(z,MCMC){
  return(MCMC$theta[3,]*(1+MCMC$theta[1,]*(z-MCMC$u)/MCMC$theta[2,])^(-1/MCMC$theta[1,]))
}

# narrowest CI for future predicted value
predCI<-function(mod,prob=1e-6,alpha=0.9,nACER=100,nmcmc=10000){
  minCI<-c(0,.Machine$integer.max)
  tempCI<-c(NA,NA)
  if(class(mod)=='ACER'){
    low<-(1-alpha)/2
    up<-(1+alpha)/2
    rseq<-seq(0,1-alpha,length.out=nACER+2)[2:(nACER+1)]
    for(i in 1:nACER){
      tempCI[1]<-reverseGeneral(eps=1-rseq[i]^prob,est=mod$coef[1:5])
      tempCI[2]<-reverseGeneral(eps=1-(rseq[i]+alpha)^prob,est=mod$coef[1:5])
      if(diff(minCI)>diff(tempCI)){
        minCI<-tempCI
      }
    }
  }else if(class(mod)=='mcmc'){
    dist<-predDist(mod=mod,prob=prob,nmcmc=nmcmc)
    minCI<-distCI(dist=dist,alpha=alpha)
  }
  return(minCI)
}

# Narrowest CI from a numerical distribution
distCI<-function(dist,alpha=0.9){
  minCI<-c(0,.Machine$integer.max)
  tempCI<-c(NA,NA)
  dist<-sort(dist)
  ndist<-length(dist)
  interval<-round(ndist*alpha)
  for(i in 1:(ndist-interval)){
    tempCI<-c(dist[i],dist[i+interval])
    if(diff(minCI)>diff(tempCI)){
      minCI<-tempCI
    }
  }
  return(minCI)
}

# Distribution of future predicted value
predDist<-function(mod,prob=1e-6,nACER=10000,nmcmc=10000){
  if(class(mod)=='ACER'){
    u<-seq(0,1,length.out=nACER+2)[2:(nACER+1)]#runif(nACER)
    return(reverseGeneral(eps=1-u^prob,est=mod$coef[1:5]))
  }else if(class(mod)=='mcmc'){
    n<-length(mod$theta[1,])
    nloop<-ceiling(nmcmc/n)
    dist<-rep(NA,nloop*n)
    for(i in 1:nloop){
      dist[(1+(i-1)*n):(i*n)]<-reversePOTmcmc(prob=1-runif(n)^prob,MCMC=mod)
    }
    return(dist)
  }
}

acerMcmcPlot<-function(acerMod,mcmcMod,acerCol,mcmcCol,ylim,xlim,n=200,xlab,ylab,alpha=0.95){
  eta<-seq(xlim[1],xlim[2],length.out=n)
  plot(eta,General(eta,acerMod$coef),type='l',col=acerCol,ylim=ylim,xlim=xlim,xlab=xlab,ylab=ylab,log='y')
  lines(eta,General(eta,acerMod$upperCIcoef),lty=2,col=acerCol)
  lines(eta,General(eta,acerMod$lowerCIcoef),lty=2,col=acerCol)
  
  meanCI<-matrix(NA,3,n)
  for(i in 1:n){
    dist<-POTmcmc(eta[i],mcmcMod)
    meanCI[2,i]<-mean(dist)
    meanCI[c(1,3),i]<-distCI(dist,alpha=alpha)
  }
  lines(eta,meanCI[2,],col=mcmcCol)
  lines(eta,meanCI[1,],lty=2,col=mcmcCol)
  lines(eta,meanCI[3,],lty=2,col=mcmcCol)
}