###
# Havy tail distribuiton
###########################
n<-365
u<-runif(n*10)
alpha<-3
data<-matrix((1/(u))^(1/alpha),10,n,byrow=TRUE) #10 realisations (one for each day over 10 years)

plot(1:200,(1:200)^(-alpha),log='y',type='l')


# Speed Test
test<-mcmc.gpd(X=data,u=1,n=1000,cpp=FALSE)
test<-mcmc.gpd(X=test,n=1000,cpp=FALSE)
plot(test$theta[1,],type='l')
testC<-mcmc.gpd(X=data,u=1,n=1000,cpp=TRUE)
testC<-mcmc.gpd(X=testC,n=1000,cpp=TRUE)
lines(testC$theta[1,],col=2)
effsampSize(test$theta[1,])
effsampSize(testC$theta[1,])

Rprof("out.out")
test<-mcmc.gpd(X=data,u=1,mu=mu,var=var,n=10000,gam=0.1,cpp=FALSE)
Rprof(NULL)
summaryRprof("out.out")

Rprof("out.out")
testC<-mcmc.gpd(X=data,u=1,mu=mu,var=var,n=10000,gam=0.1,cpp=TRUE)
Rprof(NULL)
summaryRprof("out.out")
# Speed Test



library(ACER)
mod1<-ACER(data,k=1)
plot(mod1,xlim=c(0,200))
lines(1:200,(1:200)^(-alpha),col='green')


#Rcpp
library(Rcpp)
sourceCpp("cpp/MCMC_GPD_CPP.cpp")
(muC<-c(1,1))
(varC<-matrix(c(1,0.2,0.2,1),2,2))
test=mcmcGpdC(y=1:10,nstart=1,n=10,start=c(1,2),muR=mu,varR=var, tau=0,a=0, gam=0.1,lamda=0)
test2=mcmcGpd.internal(y=1:10,nstart=1,n=10,start=c(1,2),mu=muC,var=varC, tau=0,a=0, gam=0.1,lamda=0)



#MCMC
# 1.Declustering (r)
# 2.max within each cluster
# ?(3.) find u (peak over threshold)
# 4. assume each kluster to be independent with pareto distribution
# 5.fitting cluster extreme through general pareto dist
# Need both u and r


# In this case u=1 (Kai Erik)
# Priors:
# pi=log(sigma)
# fmy,fpi,fxi ~ normal
# var of pi and mu = 1e04, var of xi = 100
# random walk:
# my*=my+Emy, pi*=.....
# (9.3) seqentually for each param
# GEV => 3 param, while POT GPD => 2 (+n~poisson(lambda))

#var of random walk for each with generator
n=1000
vmy<-10
vpi<-10
vxi<-10
rnorMy<-rnorm(n,0,vmy)
rnorPi<-rnorm(n,0,vpi)
rnorXi<-rnorm(n,0,vxi)

my<-rep(0,n)
pi<-rep(0,n)
xi<-rep(0,n)

# R=f(X)/f(x) (symetric)=LGPD*poisson(lamda)/LGPD*poisson(lamda)=exp(lnGPD+lnPoi-lnGPD+lnPoi)
#for GPD LGPD*poisson(lamda)=f(X)






lnGPD->function(data,si,xi){
  return(-length(data)*log(si)-(1+1/xi)*sum(log(1+xi*data/si)))
}

lnGEV->function(data,my,si,xi){
  temp<-1+xi*(data-my)/si
  return(-length(data)*log(si)-(1+1/xi)*sum(log(temp))-sum(temp^(-1/xi)))
}
