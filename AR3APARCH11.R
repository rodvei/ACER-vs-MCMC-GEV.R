### GENERATED AR(3)-APARCH(1,1)
coil<-read.table("plotcode/coil.txt",header=TRUE)
plot(coil[,2],type='l')

testA3<-garchFit(formula=~arma(3,0)+aparch(1,1),cond.dist="sstd",data=coil[,2],leverage=TRUE)#,include.mean=FALSE)

model3=list(ar=c(-2.546e-02,-4.335e-02,-1.723e-02),beta=0.95,mu=4.126e-04,omega=9.769e-05,alpha=5.545e-02,gamma=2.333e-01,delta=1.127e+00,skew=9.617e-01,shape=7.342e+00)
garchSpecModel3<-garchSpec(model=model3,cond.dist = "sstd",presample = cbind(z=testA3@fit$series$z[1:10],h=testA3@fit$series$h[1:10],y=testA3@fit$series$x[1:10]))

plot(garchSim(garchSpecModel3,n=length(coil[,2]))$garch,type='l',xlab='t',ylab=expression(z[t]),mgp=c(2,1,0),cex.lab=1.4)
