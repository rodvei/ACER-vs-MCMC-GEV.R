###
# Havy tail distribuiton
###########################
n<-100000
u<-runif(n)
alpha<-2
data<-(1/(1-u))^(1/alpha)
plot(1:n,sort(data),log='y',type='l')
