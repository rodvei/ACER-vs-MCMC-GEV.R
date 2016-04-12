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

###################################################
############## Pareto analysis ####################
###################################################
day<-365
year<-25
set.seed(5000)
u<-runif(day*year)
a<-3
data<-matrix((1/(u))^(1/a),year,day,byrow=TRUE) #10 realisations (one for each day over 10 years)
plot(as.vector(data),type='h',ylab=bquote(paste(x[i])),xlab='i')
#points((1:(365*25))[as.vector(data)>4],as.vector(data)[as.vector(data)>4],pch=16,cex=0.5)
#points(as.vector(data),pch=16,cex=0.5)

plot(1:200,1/(1:200)^(a),log='y',type='l',ylim=c(1e-7,1),ylab='P(X>z)', xlab='z')

pYearly<-1/365
pDecade<-1/(10*365)
pCentury<-1/(100*365)
pMillennium<-1/(1000*365)

library(ACER)
acerMod<-ACER(data,k=1,eta1=1.5,CI=0.95,stationary=TRUE,method='general', check.weibull = FALSE)
mcmcMod<-mcmc.gpd(as.vector(data),u=1.5,n=10000,start=c(1,1),a=0.35)
old.par <- par(mfrow=c(2, 1))
plot(mcmcMod$theta[1,1:1000],type='l',ylab=expression(xi^t),xlab='t',mgp=c(2,1,0))
plot(mcmcMod$theta[2,1:1000],type='l',ylab=expression(phi^t),xlab='t',mgp=c(2,1,0))
par(old.par)
aRate(mcmcMod$theta[2,350:10000])



mcmcMod<-mcmc.gpd(as.vector(data),u=1.5,n=100000)
mcmcMod$theta=mcmcMod$theta[1:3,500:100000]
effsampSize(mcmcMod$theta[1,])
effsampSize(mcmcMod$theta[2,])
# plot est Pr(X>z)
acerMcmcPlot(acerMod=acerMod,mcmcMod=mcmcMod,acerCol='blue',mcmcCol='red',ylim=c(1e-8,1e-3),xlim=c(10,150),n=200,ylab='Pr(Z>z)', xlab='z',alpha=0.95)
lines(seq(10,150,length.out=200),1/(seq(10,150,length.out=200))^(3))
legend("topright", inset=0,c('Pareto Dist','ACER','POT MCMC'), fill=c('black','blue','red'), horiz=FALSE)

old.par <- par(mfrow=c(1, 3))
plot(density(mcmcMod$theta[1,],bw=0.005),xlab=expression(xi),main='')
plot(density(mcmcMod$theta[2,],bw=0.005),xlab=expression(sigma),main='')
plot(density(mcmcMod$theta[3,],bw=0.005),xlab=expression(beta),main='')
par(old.par)



# plot year est Pr(M_n=z)  # '#FF003322' first 6 is color, last two is transparancy.
m <- matrix(c(1,2,3,4,5,5),nrow = 3,ncol = 2,byrow = TRUE)
layout(mat = m,heights = c(0.4,0.4,0.2))
#old.par<-par(mar = c(2,2,1,1))
xlim=c(1,30)
ylim=NULL
z=seq(xlim[1],xlim[2],length.out=200)
plot(z,paretoProb(z=z,prob=pYearly,alpha=a),type='l',col=1,ylim=ylim,ylab=bquote(paste('Pr(',M[year],' = z)')),xlab='z',mgp=c(2,1,0),cex.lab=1.8)
lines(density(predDist(acerMod,prob=pYearly,nACER=10000),bw=0.25,from=xlim[1],to=xlim[2]),col='blue')
abline(v=predCI(acerMod,prob=pYearly,alpha=0.90,nACER=100),lty=2,col='blue')
lines(density(predDist(mcmcMod,prob=pYearly,nmcmc=99501),bw=0.25,from=xlim[1],to=xlim[2]),col='red')
abline(v=predCI(mcmcMod,prob=pYearly,alpha=0.90,nmcmc=99501),lty=2,col='red')
#legend("topright", inset=0,c('Pareto Dist','ACER','POT MCMC'), fill=c('black','blue','red'), horiz=FALSE,cex = 0.75)

xlim=c(1,60)
ylim=NULL
z=seq(xlim[1],xlim[2],length.out=200)
plot(z,paretoProb(z=z,prob=pDecade,alpha=a),type='l',col=1,ylim=ylim,ylab=bquote(paste('Pr(',M[decade],' = z)')),xlab='z',mgp=c(2,1,0),cex.lab=1.8)
lines(density(predDist(acerMod,prob=pDecade,nACER=10000),bw=0.25,from=xlim[1],to=xlim[2]),col='blue')
abline(v=predCI(acerMod,prob=pDecade,alpha=0.90,nACER=100),lty=2,col='blue')
lines(density(predDist(mcmcMod,prob=pDecade,nmcmc=99501),bw=0.25,from=xlim[1],to=xlim[2]),col='red')
abline(v=predCI(mcmcMod,prob=pDecade,alpha=0.90,nmcmc=99501),lty=2,col='red')
#legend("topright", inset=0,c('Pareto Dist','ACER','POT MCMC'), fill=c('black','blue','red'), horiz=FALSE,cex = 0.75)
#par(old.par)

#old.par <- par(mfrow=c(2, 1))
xlim=c(1,150)
ylim=NULL
z=seq(xlim[1],xlim[2],length.out=200)
plot(z,paretoProb(z=z,prob=pCentury,alpha=a),type='l',col=1,ylim=ylim,ylab=bquote(paste('Pr(',M[century],' = z)')),xlab='z',mgp=c(2,1,0),cex.lab=1.8)
lines(density(predDist(acerMod,prob=pCentury,nACER=10000),bw=0.25,from=xlim[1],to=xlim[2]),col='blue')
abline(v=predCI(acerMod,prob=pCentury,alpha=0.90,nACER=100),lty=2,col='blue')
lines(density(predDist(mcmcMod,prob=pCentury,nmcmc=99501),bw=0.25,from=xlim[1],to=xlim[2]),col='red')
abline(v=predCI(mcmcMod,prob=pCentury,alpha=0.90,nmcmc=99501),lty=2,col='red')
#legend("topright", inset=0,c('Pareto Dist','ACER','POT MCMC'), fill=c('black','blue','red'), horiz=FALSE)


xlim=c(1,300)
ylim=NULL
z=seq(xlim[1],xlim[2],length.out=200)
plot(z,paretoProb(z=z,prob=pMillennium,alpha=a),type='l',col=1,ylim=ylim,ylab=bquote(paste('Pr(',M[millennium],' = z)')),xlab='z',mgp=c(2,1,0),cex.lab=1.8)
lines(density(predDist(acerMod,prob=pMillennium,nACER=10000),bw=0.25,from=xlim[1],to=xlim[2]),col='blue')
abline(v=predCI(acerMod,prob=pMillennium,alpha=0.90,nACER=100),lty=2,col='blue')
lines(density(predDist(mcmcMod,prob=pMillennium,nmcmc=99501),bw=0.25,from=xlim[1],to=xlim[2]),col='red')
abline(v=predCI(mcmcMod,prob=pMillennium,alpha=0.90,nmcmc=99501),lty=2,col='red')
#legend("topright", inset=0,c('Pareto Dist','ACER','POT MCMC'), fill=c('black','blue','red'), horiz=FALSE)
#par(old.par)
par(mar = c(2,2,1,1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x="top", inset=0,legend=c('Pareto Dist ','ACER ','POT MCMC '), fill=c('black','blue','red'),lwd=1,horiz = TRUE)
par(mar=c(5.1,4.1,4.1,2.1))

futureCIcontain(prob=pYearly,alpha=3,PI=predCI(mcmcMod,prob=pYearly,alpha=0.90,nmcmc=99501))
futureCIcontain(prob=pDecade,alpha=3,PI=predCI(mcmcMod,prob=pDecade,alpha=0.90,nmcmc=99501))
futureCIcontain(prob=pCentury,alpha=3,PI=predCI(mcmcMod,prob=pCentury,alpha=0.90,nmcmc=99501))
futureCIcontain(prob=pMillennium,alpha=3,PI=predCI(mcmcMod,prob=pMillennium,alpha=0.90,nmcmc=99501))


futureCIcontain(prob=pYearly,alpha=3,PI=predCI(acerMod,prob=pYearly,alpha=0.90,nACER=100))
futureCIcontain(prob=pDecade,alpha=3,PI=predCI(acerMod,prob=pDecade,alpha=0.90,nACER=100))
futureCIcontain(prob=pCentury,alpha=3,PI=predCI(acerMod,prob=pCentury,alpha=0.90,nACER=100))
futureCIcontain(prob=pMillennium,alpha=3,PI=predCI(acerMod,prob=pMillennium,alpha=0.90,nACER=100))

diff(predCI(mcmcMod,prob=pYearly,alpha=0.90,nmcmc=99501))
diff(predCI(mcmcMod,prob=pDecade,alpha=0.90,nmcmc=99501))
diff(predCI(mcmcMod,prob=pCentury,alpha=0.90,nmcmc=99501))
diff(predCI(mcmcMod,prob=pMillennium,alpha=0.90,nmcmc=99501))

print('ACER START')
diff(predCI(acerMod,prob=pYearly,alpha=0.90,nACER=100))
diff(predCI(acerMod,prob=pDecade,alpha=0.90,nACER=100))
diff(predCI(acerMod,prob=pCentury,alpha=0.90,nACER=100))
diff(predCI(acerMod,prob=pMillennium,alpha=0.90,nACER=100))



###################################################
############### SIM ANALYSE #######################
###################################################
coil<-read.table("plotcode/coil.txt",header=TRUE)
plot(coil[,2],type='l')
n<-length(coil[,2])
noMeanCoil<-coil[,2]-mean(coil[,2])
#gamma<-(diff(coil[,2])>0)*1
library(fGarch)
test0<-garchFit(formula=~arma(1,0)+aparch(1,1),cond.dist="sstd",data=noMeanCoil,leverage=TRUE,include.mean=FALSE)
test0@fit$ics
test0@fit$llh

testA2<-garchFit(formula=~arma(2,0)+aparch(1,1),cond.dist="sstd",data=noMeanCoil,leverage=TRUE,include.mean=FALSE)
testA2@fit$ics
testA2@fit$llh

testA3<-garchFit(formula=~arma(3,0)+aparch(1,1),cond.dist="sstd",data=noMeanCoil,leverage=TRUE,include.mean=FALSE)
testA3@fit$ics
testA3@fit$llh

testA4<-garchFit(formula=~arma(4,0)+aparch(1,1),cond.dist="sstd",data=noMeanCoil,leverage=TRUE,include.mean=FALSE)
testA4@fit$ics
testA4@fit$llh

2*(-testA4@fit$llh+testA3@fit$llh)

# AIC=>A3
# BIC=>A2
# SIC=>A4
# HQIC=>A$



garchSim(garchSpec(model=test@fit$params$params,cond.dist="sstd"),n=100)

plot(garchSim(garchSpec(model=test@fit$params$params,cond.dist="sstd"),n=10000)+mean(mean(coil[2:n,2])),type='l')

# residuals:
# (standardize=FALSE) difference between Yt and its conditional expectation (est \hat{a}_t)
# (standardize=TRUE) \hat{eps}_t=\hat{a}_t / \hat{\sigma}_t
plot(residuals(test, standardize = FALSE),type='l')








#EKSTRA??
acermod<-ACERm(data,k=1,stationary=TRUE)
acermod<-CI(acermod,level=alpha[i])
acermod<-updatecond(acermod,eta1=2)
acermod<-MSEoptimization(acermod,metgid='general')


# effective sample sice vs accaptance rate a 
eff=rep(NA,31)
for(i in 1:31){
  print(i)
  u<-runif(day*year)
  data<-(1/(u))^(1/3)
  temp1=effsampSize(mcmc.gpd(as.vector(data),u=2,n=10000,start=c(1,1),a=0.19+i/100,cpp=TRUE)$theta[1,500:10000])
  u<-runif(day*year)
  data<-(1/(u))^(1/3)
  temp2=effsampSize(mcmc.gpd(as.vector(data),u=2,n=10000,start=c(1,1),a=0.19+i/100,cpp=TRUE)$theta[1,500:10000])
  eff[i]=(temp1+temp1)/2
}
plot(20:50,eff,type='p')
