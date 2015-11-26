###
# Havy tail distribuiton
###########################
n<-3650
u<-runif(n)
alpha<-3
data<-matrix((1/(u))^(1/alpha),10,3650) #10 realisations (one for each day over 10 years)

plot(1:200,(1:200)^(-alpha),log='y',type='l')

library(ACER)
mod1<-ACER(data,k=2)
plot(mod1,xlim=c(0,200))
lines(1:200,(1:200)^(-alpha),col='green')





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

