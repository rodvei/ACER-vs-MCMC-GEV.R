#library(ismev)
#?mrl.plot
#?gpd.fitrange
# from page 83 stuart
# \xi =constant
# \sigma=lin, \sigma*=\sigma-\xi*u=const




##################### MY OWN #####################

# plot mean exceedance vs u (thresholds) with confidence lines.
# For valid u's, need linearity.
# X=data set
# r=numbers of data below threshold, before defined as new cluster.
# k=plot ends when only k extreme is left
# iid Independent and Identical Distributed => r=0
# numPoint= number of points used in plot
uplot<-function(X,k=3,r=1,xlim=NULL, CI=0.95,numPoint=100,...){
  if(CI>1||CI<0){stop("CI should be 0 < CI < 1")}
  if(!(is.vector(xlim)&&length(xlim)==2)){
    Xsort<-sort(X)
    if(is.null(xlim)){
      xlim=c(Xsort[floor(length(Xsort)*0.75)],NA)
    }else{
      xlim<-c(xlim,NA)
    }
    if(r==0){
      xlim[2]<-Xsort[length(Xsort)-k+1]
    }else{
      xlim[2]<-umax(X,k=k,r=r)
    }
  }
  uStep<-seq(from=xlim[1],to=xlim[2], length=numPoint)
  ymean<-rep(NA,numPoint)
  tSD<-rep(NA,numPoint)
  #the variable is tSD=t*sd
  lowCI<-rep(NA,numPoint)
  pb <- txtProgressBar(min = 0, max = numPoint, style = 3)
  for(i in 1:numPoint){
    ytemp<-ufilt(X,uStep[i],r=r)-uStep[i]
    ymean[i]=mean(ytemp)
    tSD[i]=qt(p=(1-CI)/2, df=length(ytemp)-1, lower.tail = FALSE)*sd(ytemp)/sqrt(length(ytemp))
    setTxtProgressBar(pb, i)
  }
  close(pb)
  plot(uStep,ymean,type='l',xlim=xlim,ylim=c(min(ymean-tSD),max(ymean+tSD)),xlab = "u", ylab = "Mean Excess")
  lines(uStep, ymean+tSD, col='blue', lty=2)
  lines(uStep, ymean-tSD, col='blue', lty=2)
}

# parameter plot (xi vs u and sigma*=sigma-xi*u \sim const vs u)
# maxtime= maximum time calculating at each u (by defualt 5 seconds)
# numPoint= number of points used in plot
parplot<-function(X,k=3,r=1,xlim=NULL,CI=0.95,maxtime=5,burnin=500, numPoint=15,...){
  if(CI>1||CI<0){stop("CI should be 0 < CI < 1")}
  if(!(is.vector(xlim)&&length(xlim)==2)){
    Xsort<-sort(X)
    if(is.null(xlim)){
      xlim=c(Xsort[floor(length(Xsort)*0.75)],NA)
    }else{
      xlim<-c(xlim,NA)
    }
    if(r==0){
      xlim[2]<-Xsort[length(Xsort)-k+1]
    }else{
      xlim[2]<-umax(X,k=k,r=r)
    }
  }
  
  names<-c(expression(paste("shape (", xi,")")), expression(paste("Modified Scale (", sigma,"*)")))
  uStep<-seq(from=xlim[2],to=xlim[1], length=numPoint)
  para<-matrix(NA,2,numPoint) # para[1,]=xi while para[2,]=sigma*
  upCI<-matrix(NA,2,numPoint)
  lowCI<-matrix(NA,2,numPoint)
  
  
  tempmcmc<-NA
  tempstart<-rep(NA,2)
  tempdata<-matrix(NA,2,9000)
  burnintemp=burnin
  
  pb <- txtProgressBar(min = 0, max = numPoint, style = 3)
  for(j in 1:numPoint){
    if(j==1){
      tempmcmc<-mcmc.gpd(X,u=uStep[j],n=10000)
    }else{
      if((sum(X>uStep[j])>(maxtime*1200))&&(sum(X>uStep[j-1])>30)){
        timescale<-maxtime*1200/(sum(X>uStep[j]))
        tempmcmc<-mcmc.gpd(X=X,u=uStep[j],start=tempstart,n=floor(timescale*10000),mu=tempmcmc$mu,var=tempmcmc$var)
        burnintemp=floor(0.2*(timescale*10000))
      }else if(sum(X>uStep[j-1])>30){
        tempmcmc<-mcmc.gpd(X=X,u=uStep[j],start=tempstart,n=10000,mu=tempmcmc$mu,var=tempmcmc$var)
        burnintemp=burnin
      }else{
        tempmcmc<-mcmc.gpd(X=X,u=uStep[j],n=10000)
        burnintemp=burnin
      }
    }

    tempstart<-tempmcmc$theta[1:2,tempmcmc$n]
    tempsort<-rbind(sort(tempmcmc$theta[1,burnintemp:tempmcmc$n]),sort(tempmcmc$theta[2,burnintemp:tempmcmc$n]-tempmcmc$theta[1,burnintemp:tempmcmc$n]*uStep[j]))
    
    para[1:2,j]<-c(mean(tempsort[1,]),mean(tempsort[2,]))
    upCI[1:2,j]<-c(tempsort[1,floor((tempmcmc$n-burnintemp)*(1+CI)/2)],tempsort[2,floor((tempmcmc$n-burnintemp)*(1+CI)/2)])
    lowCI[1:2,j]<-c(tempsort[1,floor((tempmcmc$n-burnintemp)*(1-CI)/2)],tempsort[2,floor((tempmcmc$n-burnintemp)*(1-CI)/2)])
    setTxtProgressBar(pb, j)
  }
  close(pb)
  oldpar <- par(mfrow = c(2, 1))
  for(i in 1:2){
    plot(uStep,para[i,], xlab='u', ylab=names[i],ylim=c(min(lowCI[i,]),max(upCI[i,])),pch=16,type='b')
    arrows(uStep,lowCI[i,],uStep,upCI[i,],code=3,length=0.1,angle=90,col='blue')
  }
  par(oldpar)
}


# Filtering out threshold exceedance y, which is the max of each cluster exceeding the treshold.
# Cluster is define where r value is below the threshold (repeatedly)
# x=data set (vector)
# u=Threshold
# r=numbers of data below, before defined as new cluster 
# y=return x-u
ufilt<-function(X,u,r=1){
  if(r==0){
    y<-X[X>=u]
  }else{
    y<-rep(NA,sum(X>=u))
    max<-u*(1-2*.Machine$double.eps)
    rcount<-0
    clust<-FALSE
    j<-1
    for(i in 1:length(X)){
      if(X[i]>=u){
        rcount<-0
        clust<-TRUE
        if(X[i]>max){
          max<-X[i]
        }
      }else if(X[i]<u&&clust){
        rcount<-rcount+1
      }
      if(rcount>=r&&max>=u){
        y[j]<-max
        clust<-FALSE
        rcount<-0
        max<-u*(1-2*.Machine$double.eps)
        j<-j+1
      }
    }
    if(max>=u){
      y[j]<-max
    }
    y=y[!is.na(y)]
  }
  return(y)
}

# find u-max which gives only 
# k block exceedance
# r=numbers of data below threshold, before defined as new cluster.
# return u-max
umax<-function(X,k,r){
  Xmean<-mean(X)
  uLeng<-0
  Xsort<-sort(X)
  n <- length(Xsort)
  i=-2
  while(uLeng<k&&Xsort[n-k-i]>Xmean){
    i=i+1
    uLeng<-length(ufilt(X,u=Xsort[n-k-i],r))
  }
  if(uLeng>=k){
    return((Xsort[n-k-i]+Xsort[n-k-i-1])/2)
  }else{
    return(umax(X,k-1,r))
  }
}
