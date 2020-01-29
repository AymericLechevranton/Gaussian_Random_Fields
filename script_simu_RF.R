#     Simulation avec geoR

library(geoR)
library(spam)
library(fields)
library(gstat)
library(Rcpp)
library(pracma)
library(RandomFieldsUtils)
library(RGeostats)
#demo(RGeostats.start)


    # Semivariogram of different covariances

x=seq(0,1, length.out=100)
a=0.6 #range
sigma=1 #sill

#Pure Nugget effect
f1=function(x,sigma){sigma*(1*(x>=0))} #we use greater or equal to have not continuous where it is discontinuous
plot(x,f1(x,sigma), lty=)
#Spherical
f2=function(x,sigma,a){sigma*(1.5*x/a-0.5*(x/a)^3)*(x<=a)+sigma*(x>a)}
#Gaussian
f3=function(x,sigma,a){sigma*(1-exp(-(x/a)^2))}
#Exponential
f4=function(x,sigma,a){sigma*(1-exp(-(x/a)))}

#Variogram théorique avec portées effectives
matplot(x,cbind(f1(x,sigma),f2(x,sigma,a),f3(x,sigma,a),f4(x,sigma,a)), type="l", 
        lty=2:5, col=2:5, xlab="h", ylab=expression(paste(gamma,"(h)")))
abline(v=0.6, col="gray", lwd=0.5)
legend("bottomright", legend=c("Pure nugget","Sperical","Gaussian","Exponential"), 
       lty=2:5, col=2:5)

#Variogram théorique avec portées ajustées
matplot(x,cbind(f1(x,sigma),f2(x,sigma,a),f3(x,sigma,a*0.5777614),f4(x,sigma,a*0.3338082)), type="l", 
        lty=2:5, col=2:5, xlab="h", ylab=expression(paste(gamma,"(h)")))
abline(v=0.6, col="gray", lwd=0.5)
legend("bottomright", legend=c("Pure nugget","Sperical","Gaussian","Exponential"), 
       lty=2:5, col=2:5)







#################### Cholesky decomposition ###################

# Covariance function: pure nugget (weack white noise)
set.seed(25)
n=50
tab=expand.grid(x=seq(0,200,,n+1),y=seq(0,200,,n+1))
d = as.matrix(dist(tab,diag=T,upper=T))
K=eye(dim(d)[1])
L = chol(K)
#Calculate Z=LY with Y a standard gaussian sample
Z = L %*% rnorm(((n+1)^2))
#Distribution
hist(Z, freq=F)
lines(x=seq(-3,3,,50), y=dnorm(seq(-3,3,,50)), col='red')
#Graphic
res <- list(x=seq(0,200,,n+1),y=seq(0,200,,n+1),z=matrix(as.numeric(Z),ncol=n+1, byrow=T))
image(res,main = 'Pure nugget')
#variogram
dbv = db.create(flag.grid=T,nx=c(n+1,n+1),dx=c(4,4),x0=c(0,0))
dbv = db.add(db=dbv,array(res$z),loctype="z")
var = vario.grid(dbv,dirincr=c(0,1),nlag=100)
plot(var,col="gray",lwd=1.5,ylim=c(0,1.5),xlab="Distance",ylab="Semi-variance",main="Pure nugget effect variogram")
abline(h=1,col='red',lty=2);legend("bottomright",lwd=c(1.5,1,1),lty=c(1,2,2),col=c("gray",1,"red"),legend=c("Experimental variogram","Average of the experimental variogram","Theoritical variogram"))


# Covariance function: spherical
set.seed(25)
n=50
#Range
a=20
#Create regular grid
tab=expand.grid(x=seq(0,200,,n+1),y=seq(0,200,,n+1))
d = as.matrix(dist(tab,diag=T,upper=T))
K=(1-(1.5*(d/a)-0.5*((d/a)^3)))*(d<a)
L = chol(K)
#Calculate Z=LY with Y a standard gaussian sample
Z = L %*% rnorm((n+1)^2)
#Distribution
hist(Z, freq=F)
lines(x=seq(-3,3,,50), y=dnorm(seq(-3,3,,50)), col='red')
#Graphic
res <- list(x=seq(0,200,,n+1),y=seq(0,200,,n+1),z=matrix(as.numeric(Z),ncol=n+1, byrow=T))
image(res,main = paste('Spherical of range',a))
#variogram
dbv = db.create(flag.grid=T,nx=c(n+1,n+1),dx=c(4,4),x0=c(0,0))
dbv = db.add(db=dbv,array(res$z),loctype="z")
var = vario.grid(dbv,dirincr=c(0,1),nlag=100)
plot(var,col="gray",lwd=1.5, ylim=c(0,1.5),xlab="Distance",ylab="Semi-variance",main="Spherical variogram")
x1=seq(0,a,,a+1);lines(x=x1,y=(3/2)*(x1/a)-0.5*((x1/a)^3),col='red',lty=2)
x2=seq(a+1,200,,200-a);lines(x=x2,y=rep_len(1,length(x2)),col='red',lty=2)
legend("bottomright",lwd=c(1.5,1,1),lty=c(1,2,2),col=c("gray",1,"red"),legend=c("Experimental variogram","Average of the experimental variogram","Theoritical variogram"))


# Covariance function: gaussian
set.seed(25)
n=50
#Range
#we couldn't go high enough due to the approximation of the cholesky decomposition!
a=sqrt(-(20^2)/log(0.05)) #In order to have the same range (that sperical)
#Create regular grid
tab=expand.grid(x=seq(0,200,,n+1),y=seq(0,200,,n+1))
d = as.matrix(dist(tab,diag=T,upper=T))
K=exp(-((d/a)^2))
L = chol(K)
#Calculate Z=LY with Y a standard gaussian sample
Z = L %*% rnorm(((n+1)^2))
#Distribution
hist(Z, freq=F)
lines(x=seq(-3,3,,50), y=dnorm(seq(-3,3,,50)), col='red')
#Graphic
res <- list(x=seq(0,200,,n+1),y=seq(0,200,,n+1),z=matrix(as.numeric(Z),ncol=n+1, byrow=T))
image(res,main = 'Gaussian of range 20')
#variogram
dbv = db.create(flag.grid=T,nx=c(n+1,n+1),dx=c(4,4),x0=c(0,0))
dbv = db.add(db=dbv,array(res$z),loctype="z")
var = vario.grid(dbv,dirincr=c(0,1),nlag=100)
plot(var,col="gray",lwd=1.5, ylim=c(0,1.5),xlab="Distance",ylab="Semi-variance",main="Gaussian variogram")
x=seq(0,200,,n+1);lines(x=x,y=1-exp(-((x/a)^2)),col='red',lty=2)
legend("bottomright",lwd=c(1.5,1,1),lty=c(1,2,2),col=c("gray",1,"red"),legend=c("Experimental variogram","Average of the experimental variogram","Theoritical variogram"))


# Covariance function: exponential
set.seed(25)
n=50
#Range
a=(-20/log(0.05))
#Create regular grid
tab=expand.grid(x=seq(0,200,,n+1),y=seq(0,200,,n+1))
#ptm<-proc.time()
d = as.matrix(dist(tab,diag=T,upper=T))
K=exp(-(d/a))
L = chol(K)
#Calculate Z=LY with Y a standard gaussian sample
Z = L %*% rnorm(((n+1)^2))
#proc.time() -ptm
#Distribution
hist(Z, freq=F,ylim=c(0,0.4),xlab="z")
lines(x=seq(-3,3,,50), y=dnorm(seq(-3,3,,50)), col='red')
#Graphic
res <- list(x=seq(0,200,,n+1),y=seq(0,200,,n+1),z=matrix(as.numeric(Z),ncol=n+1, byrow=T))
image(res,main='Exponential of range 20')
#variogram
dbv = db.create(flag.grid=T,nx=c(n+1,n+1),dx=c(4,4),x0=c(0,0))
dbv = db.add(db=dbv,array(res$z),loctype="z")
var = vario.grid(dbv,dirincr=c(0,1),nlag=100)
plot(var,col="gray",lwd=1.5, ylim=c(0,1.5),xlab="Distance",ylab="Semi-variance",main="Exponential variogram")
x=seq(0,200,,n+1);lines(x=x,y=1-exp(-(x/a)),col='red',lty=2)
legend("bottomright",lwd=c(1.5,1,1),lty=c(1,2,2),col=c("gray",1,"red"),legend=c("Experimental variogram","Average of the experimental variogram","Theoritical variogram"))



#   Pour connaitre le temps d'execution :
#ptm<-proc.time()
#proc.time() -ptm





################# Average variograms #####################

#Average of variograms of the exponential model of range 28
set.seed(23)
n=50
a=(-28/log(0.05))
tab=expand.grid(x=seq(0,200,,n+1),y=seq(0,200,,n+1))
d = as.matrix(dist(tab,diag=T,upper=T))
K=exp(-(d/a))
L = chol(K)
result=NA
for(i in 1:20){
  Z = L %*% rnorm(((n+1)^2))
  res = list(x=seq(0,200,,n+1),y=seq(0,200,,n+1),z=matrix(as.numeric(Z),ncol=n+1, byrow=T))
  geo = variog(as.geodata(cbind(tab,array(res$z)), coords.col = 1:2, data.col = 3))
  result = cbind(result,geo$v)
}
result=result[,-1]
matplot(x=geo$u,y=result, type = "l",lty=3,ylim=c(0,1.5),col=1,xlab="Distance",ylab="Semi-variance",main="Exponential variogram of range 28")
quant = apply(result,1,quantile,probs=c(0.05,0.95)) #confidence intervals
polygon(x=c(geo$u,geo$u[13:1]),y=c(quant[1,],quant[2,13:1]),col=rgb(0,0,0.8,0.3),lwd=0.3,border=NA)
lines(x=geo$u,y=rowMeans(result), type = "l",lty=2,lwd=1.5,col="blue")
x=seq(0,270,,n+1);lines(x=x,y=1-exp(-(x/a)),col='red',lwd=2)
abline(v=28, col="gray", lwd=0.5)
legend("bottomright", legend=c("Theoritical variogram","Average variogram","Experimental variograms","Confidence interval"), lty=c(1:3,1),lwd=c(1,1,1,10),col=c("red","blue",1,rgb(0,0,0.8,0.3)))

    #Average of variograms of the exponential model of range 80
set.seed(23)
n=50
a=(-80/log(0.05))
tab=expand.grid(x=seq(0,200,,n+1),y=seq(0,200,,n+1))
d = as.matrix(dist(tab,diag=T,upper=T))
K=exp(-(d/a))
L = chol(K)
result=NA
for(i in 1:20){
  Z = L %*% rnorm(((n+1)^2))
  res = list(x=seq(0,200,,n+1),y=seq(0,200,,n+1),z=matrix(as.numeric(Z),ncol=n+1, byrow=T))
  geo = variog(as.geodata(cbind(tab,array(res$z)), coords.col = 1:2, data.col = 3))
  result = cbind(result,geo$v)
}
result=result[,-1]
matplot(x=geo$u,y=result, type = "l",lty=3,ylim=c(0,1.5),col=1,xlab="Distance",ylab="Semi-variance",main="Exponential variogram of range 80")
quant = apply(result,1,quantile,probs=c(0.05,0.95)) #confidence intervals
polygon(x=c(geo$u,geo$u[13:1]),y=c(quant[1,],quant[2,13:1]),col=rgb(0,0,0.8,0.3),lwd=0.3,border=NA)
lines(x=geo$u,y=rowMeans(result), type = "l",lty=2,lwd=1.5,col="blue")
x=seq(0,270,,n+1);lines(x=x,y=1-exp(-(x/a)),col='red',lwd=2)
abline(v=80, col="gray", lwd=0.5)
legend("bottomright", legend=c("Theoritical variogram","Average variogram","Experimental variograms","Confidence interval"), lty=c(1:3,1),lwd=c(1,1,1,10),col=c("red","blue",1,rgb(0,0,0.8,0.3)))

    #Average of variograms of the exponential model of range 160
set.seed(23)
n=50
a=(-160/log(0.05))
tab=expand.grid(x=seq(0,200,,n+1),y=seq(0,200,,n+1))
d = as.matrix(dist(tab,diag=T,upper=T))
K=exp(-(d/a))
L = chol(K)
result=NA
for(i in 1:20){
  Z = L %*% rnorm(((n+1)^2))
  res = list(x=seq(0,200,,n+1),y=seq(0,200,,n+1),z=matrix(as.numeric(Z),ncol=n+1, byrow=T))
  geo = variog(as.geodata(cbind(tab,array(res$z)), coords.col = 1:2, data.col = 3))
  result = cbind(result,geo$v)
}
result=result[,-1]
matplot(x=geo$u,y=result, type = "l",lty=3,ylim=c(0,1.5),col=1,xlab="Distance",ylab="Semi-variance",main="Exponential variogram of range 160")
quant = apply(result,1,quantile,probs=c(0.05,0.95)) #confidence intervals
polygon(x=c(geo$u,geo$u[13:1]),y=c(quant[1,],quant[2,13:1]),col=rgb(0,0,0.8,0.3),lwd=0.3,border=NA)
lines(x=geo$u,y=rowMeans(result), type = "l",lty=2,lwd=1.5,col="blue")
x=seq(0,270,,n+1);lines(x=x,y=1-exp(-(x/a)),col='red',lwd=2)
abline(v=160, col="gray", lwd=0.5)
legend("bottomright", legend=c("Theoritical variogram","Average variogram","Experimental variograms","Confidence interval"), lty=c(1:3,1),lwd=c(1,1,1,10),col=c("red","blue",1,rgb(0,0,0.8,0.3)))








#################### Turning bands ###################

melem.name()

#Pure nugget effect:
n=50
mygrid <- db.create(flag.grid = T,nx=c(n+1,n+1),dx=c(4,4),x0=c(0,0))
mymodel <- model.create(vartype=1)
simpne <- simtub(dbout=mygrid,model=mymodel,nbtuba=1000)
plot(simpne, col=heat.colors(6))
#variogram
var = vario.grid(simpne,dirincr=c(0,1),nlag=100)
plot(var,col="gray",lwd=1.5,ylim=c(0,1.5),xlab="Distance",ylab="Semi-variance",main="Pure nugget effect variogram")
abline(h=1,col='red',lty=2);legend("bottomright",lwd=c(1.5,1,1),lty=c(1,2,2),col=c("gray",1,"red"),legend=c("Experimental variogram","Average of the experimental variogram","Theoritical variogram"))


#Spherical:
n=50
a=20
mygrid <- db.create(flag.grid = T,nx=c(n+1,n+1),dx=c(4,4),x0=c(0,0))
mymodel <- model.create(vartype=3,range=a)
simsphe <- simtub(dbout=mygrid,model=mymodel,nbtuba=1000)
plot(simsphe, col=heat.colors(6))
#variogram
var = vario.grid(simsphe,dirincr=c(0,1),nlag=100)
plot(var,col="gray",lwd=1.5, ylim=c(0,1.5),xlab="Distance",ylab="Semi-variance",main="Spherical variogram")
x1=seq(0,a,,a+1);lines(x=x1,y=(3/2)*(x1/a)-0.5*((x1/a)^3),col='red',lty=2)
x2=seq(a+1,200,,200-a);lines(x=x2,y=rep_len(1,length(x2)),col='red',lty=2)
legend("bottomright",lwd=c(1.5,1,1),lty=c(1,2,2),col=c("gray",1,"red"),legend=c("Experimental variogram","Average of the experimental variogram","Theoritical variogram"))


#Gaussian:
n=50
a=20
mygrid <- db.create(flag.grid = T,nx=c(n+1,n+1),dx=c(4,4),x0=c(0,0))
mymodel <- model.create(vartype=4,range=a)
simgaus <- simtub(dbout=mygrid,model=mymodel,nbtuba=100)
plot(simgaus, col=heat.colors(6))
#variogram
var = vario.grid(simgaus,dirincr=c(0,1),nlag=100)
plot(var,col="gray",lwd=1.5, ylim=c(0,1.5),xlab="Distance",ylab="Semi-variance",main="Gaussian variogram")
x=seq(0,200,,n+1);lines(x=x,y=1-exp(-((x/sqrt(-(a^2)/log(0.05)))^2)),col='red',lty=2)
legend("bottomright",lwd=c(1.5,1,1),lty=c(1,2,2),col=c("gray",1,"red"),legend=c("Experimental variogram","Average of the experimental variogram","Theoritical variogram"))


#Exponential:
n=50
a=20
mygrid <- db.create(flag.grid = T,nx=c(n+1,n+1),dx=c(4,4),x0=c(0,0))
#ptm<-proc.time()
mymodel <- model.create(vartype=2,range=a)
simexpo <- simtub(dbout=mygrid,model=mymodel,nbtuba=100)
#proc.time() -ptm
plot(simexpo, col=heat.colors(6))
#variogram
var = vario.grid(simexpo,dirincr=c(0,1),nlag=100)
plot(var,col="gray",lwd=1.5, ylim=c(0,1.5),xlab="Distance",ylab="Semi-variance",main="Exponential variogram")
x=seq(0,200,,n+1);lines(x=x,y=1-exp(-(x/(-a/log(0.05)))),col='red',lty=2)
legend("bottomright",lwd=c(1.5,1,1),lty=c(1,2,2),col=c("gray",1,"red"),legend=c("Experimental variogram","Average of the experimental variogram","Theoritical variogram"))

#ptm<-proc.time()
#proc.time() -ptm


hist(simexpo@items[,4],freq=F,ylim=c(0,0.4),xlab="z")
lines(x=seq(-3,3,,50), y=dnorm(seq(-3,3,,50)), col='red')




################# Average variograms #####################

#Average of variograms of the exponential model of range 28
n=50
a=28
mygrid <- db.create(flag.grid = T,nx=c(n+1,n+1),dx=c(4,4),x0=c(0,0))
mymodel <- model.create(vartype=2,range=a)
result=NA
for(i in 1:20){
  simexpo <- simtub(dbout=mygrid,model=mymodel,nbtuba=100,seed=i)
  geo = variog(as.geodata(simexpo@items[,2:4], coords.col = 1:2, data.col = 3))
  result = cbind(result,geo$v)
}
result=result[,-1]
matplot(x=geo$u,y=result, type = "l",lty=3,ylim=c(0,1.5),col=1,xlab="Distance",ylab="Semi-variance",main="Exponential variogram of range 28")
quant = apply(result,1,quantile,probs=c(0.05,0.95)) #confidence intervals
polygon(x=c(geo$u,geo$u[13:1]),y=c(quant[1,],quant[2,13:1]),col=rgb(0,0,0.8,0.3),lwd=0.3,border=NA)
lines(x=geo$u,y=rowMeans(result), type = "l",lty=2,lwd=1.5,col="blue")
x=seq(0,270,,n+1);lines(x=x,y=1-exp(-(x/(-a/log(0.05)))),col='red',lwd=2)
abline(v=28, col="gray", lwd=0.5)
legend("bottomright", legend=c("Theoritical variogram","Average variogram","Experimental variograms","Confidence interval"), lty=c(1:3,1),lwd=c(1,1,1,10),col=c("red","blue",1,rgb(0,0,0.8,0.3)))

#Average of variograms of the exponential model of range 80
n=50
a=80
mygrid <- db.create(flag.grid = T,nx=c(n+1,n+1),dx=c(4,4),x0=c(0,0))
mymodel <- model.create(vartype=2,range=a)
result=NA
for(i in 1:20){
  simexpo <- simtub(dbout=mygrid,model=mymodel,nbtuba=100,seed=i)
  geo = variog(as.geodata(simexpo@items[,2:4], coords.col = 1:2, data.col = 3))
  result = cbind(result,geo$v)
}
result=result[,-1]
matplot(x=geo$u,y=result, type = "l",lty=3,ylim=c(0,1.5),col=1,xlab="Distance",ylab="Semi-variance",main="Exponential variogram of range 80")
quant = apply(result,1,quantile,probs=c(0.05,0.95)) #confidence intervals
polygon(x=c(geo$u,geo$u[13:1]),y=c(quant[1,],quant[2,13:1]),col=rgb(0,0,0.8,0.3),lwd=0.3,border=NA)
lines(x=geo$u,y=rowMeans(result), type = "l",lty=2,lwd=1.5,col="blue")
x=seq(0,270,,n+1);lines(x=x,y=1-exp(-(x/(-a/log(0.05)))),col='red',lwd=2)
abline(v=80, col="gray", lwd=0.5)
legend("bottomright", legend=c("Theoritical variogram","Average variogram","Experimental variograms","Confidence interval"), lty=c(1:3,1),lwd=c(1,1,1,10),col=c("red","blue",1,rgb(0,0,0.8,0.3)))

#Average of variograms of the exponential model of range 160
n=50
a=160
mygrid <- db.create(flag.grid = T,nx=c(n+1,n+1),dx=c(4,4),x0=c(0,0))
mymodel <- model.create(vartype=2,range=a)
result=NA
for(i in 1:20){
  simexpo <- simtub(dbout=mygrid,model=mymodel,nbtuba=100,seed=i)
  geo = variog(as.geodata(simexpo@items[,2:4], coords.col = 1:2, data.col = 3))
  result = cbind(result,geo$v)
}
result=result[,-1]
matplot(x=geo$u,y=result, type = "l",lty=3,ylim=c(0,1.5),col=1,xlab="Distance",ylab="Semi-variance",main="Exponential variogram of range 160")
quant = apply(result,1,quantile,probs=c(0.05,0.95)) #confidence intervals
polygon(x=c(geo$u,geo$u[13:1]),y=c(quant[1,],quant[2,13:1]),col=rgb(0,0,0.8,0.3),lwd=0.3,border=NA)
lines(x=geo$u,y=rowMeans(result), type = "l",lty=2,lwd=1.5,col="blue")
x=seq(0,270,,n+1);lines(x=x,y=1-exp(-(x/(-a/log(0.05)))),col='red',lwd=2)
abline(v=160, col="gray", lwd=0.5)
legend("bottomright", legend=c("Theoritical variogram","Average variogram","Experimental variograms","Confidence interval"), lty=c(1:3,1),lwd=c(1,1,1,10),col=c("red","blue",1,rgb(0,0,0.8,0.3)))





################### Particular cases #####################

melem.name()

  #Matern variogram:
a=c(0.25,0.188,0.14,0.117,0.025)
nu=c(0.5,1,2,3,100)
x=seq(0,1,,100);plot(x=x,y=1-cov.spatial(x,cov.model="matern",cov.pars=c(1,a[1]),kappa=nu[1]),col=2,lty=1,type="l",ylab=expression(paste(gamma,"(h)")),xlab="h")
lines(x=x,y=1-cov.spatial(x,cov.model="matern",cov.pars=c(1,a[2]),kappa=nu[2]),col=3,lty=3)
lines(x=x,y=1-cov.spatial(x,cov.model="matern",cov.pars=c(1,a[3]),kappa=nu[3]),col=4,lty=4)
lines(x=x,y=1-cov.spatial(x,cov.model="matern",cov.pars=c(1,a[4]),kappa=nu[4]),col=5,lty=5)
lines(x=x,y=1-cov.spatial(x,cov.model="matern",cov.pars=c(1,a[5]),kappa=nu[5]),col=6,lty=6)
legend("bottomright",lty=c(NA,1,3,4,5,6),col=c(NA,2,3,4,5,6),
       legend=c("Matern theoritical variograms with :",expression(paste("a=0.25 and ",nu,"=0.5")),
                expression(paste("a=0.188 and ",nu,"=1")),expression(paste("a=0.14 and ",nu,"=2")),
                expression(paste("a=0.117 and ",nu,"=3")),expression(paste("a=0.025 and ",nu,"=100"))))

  #Power variogram:
#parameters
b=c(1,4,30)
c=c(1,1.5,1.9) #c always between 0 et 2
x=seq(0,1,,100);plot(x=x,y=b[1]-cov.spatial(x,cov.model="power",cov.pars=c(b[1],c[1])),col=2,lty=2,type="l",ylim=c(0,1.5),ylab=expression(paste(gamma,"(h)")),xlab="h") #(of a brownian movement)
lines(x=x,y=b[2]-cov.spatial(x,cov.model="power",cov.pars=c(b[2],c[2])),col=3,lty=2)
lines(x=x,y=b[3]-cov.spatial(x,cov.model="power",cov.pars=c(b[3],c[3])),col=4,lty=2)
legend("bottomright",lty=c(NA,2,2,2),col=c(NA,2,3,4),
       legend=c("Power theoritical variograms with :",paste("b=1 and c=1"),
                paste("b=4 and c=1.5"),paste("b=30 and c=1.9")))






