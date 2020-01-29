        # Script SSPM

#   Pour connaitre le temps d'execution :
#ptm<-proc.time()
#proc.time() -ptm


library(gstat)
library(sp)

# Exploring spatial data

    # Spatial data frames
data(meuse)
head(meuse)
coordinates(meuse) = ~x+y
bubble(meuse, "zinc", col=c("#00ff0088", "#00ff0088"), main = "zinc concentrations (ppm)")

    # Spatial data on a regular grid
data(meuse.grid)
head(meuse.grid)
coordinates(meuse.grid) = ~x+y
gridded(meuse.grid) = TRUE
image(meuse.grid["dist"])
title("distance to river (red = 0)")

zinc.idw = idw(zinc~1, meuse, meuse.grid)
spplot(zinc.idw["var1.pred"], main = "zinc inverse distance weighted interpolations")

#To linearized the relationship between zinc and distance to the river
plot(log(zinc)~sqrt(dist), meuse)
abline(lm(log(zinc)~sqrt(dist), meuse))


    # Variograms
#log(zinc)~1 means that we assume a constant trend for the variable log(zinc)
lzn.vgm = variogram(log(zinc)~1, meuse)
lzn.vgm
lzn.fit = fit.variogram(lzn.vgm, model = vgm(1, "Sph", 900, 1))
lzn.fit
plot(lzn.vgm, lzn.fit, main="Variogram") # Variogramme 





        # LU Decomposition

# If the covariance matrix is an identity matrix, it means that the grid (n by n)
# will only be a grid with N (=n*n) gaussian variables.

# If the covariance is different, 

    #Non-conditional simulation
library(Matrix)
#Create a covariance matrix
n=4
K=matrix(c(1,0.5,0.1,0.1,0.5,1,0.5,0.1,0.1,0.5,1,0.5,0.1,0.1,0.5,1), nrow=n)
K
#Decompose K matrix
dec <- expand(lu(K))
dec$L %*% dec$U #We verify the good decompostion of K.
L <- dec$L
L
#Create a standard gaussian sample
Y <- matrix(rnorm(n))
Y
#Calculate Z=LY
Z <- as.matrix( L %*% Y )
Z




    #Non-conditional simulation



library(geoR)





        # Simulating stationary Gaussian random ﬁelds

rm(list=ls())
library(RandomFields)
set.seed(1)
param <- c(0, 1, 0, 1)
model <- "exponential"
#RFparameters(PracticalRange=FALSE)
krige.method <- "O" ## random field assumption corresponding to
## those of ordinary kriging
length<-100
x <-  seq(0, 1, length=length)
y <-  seq(0, 1, length=length)
xy<-expand.grid(x,y)
data <- GaussRF(x=xy[,1], y=xy[,2], grid=FALSE, model=model, param=param)
# another grid, where values are to be simulated
# standardisation of the output
lim <- range( x)
zlim <- c(min(data)-0.1, max(data)+0.1)
colour <- gray( seq(0,1,l= 16))
data.mat<-matrix(data,length(x),length(x))
## visualise generated spatial data
par(pty='s')
image(x, y, data.mat, xlim=lim, ylim=lim, zlim=zlim, col=colour)

set.seed(19)
par(pty='s')
sub.pos<-sample(1:length(data),100)
plot(xy[sub.pos[1:25],1],xy[sub.pos[1:25],2],xlab='x',ylab='y',type='n')
points(xy[sub.pos[1:25],1],xy[sub.pos[1:25],2],pch=20)

kdata <-  Kriging(krige.method=krige.method,
                  x=xy[,1], y=xy[,2],  grid=FALSE,
                  model=model, param=param,
                  given=xy[sub.pos[1:25],], data=data[sub.pos[1:25]])
kdata.mat<-matrix(kdata,length(x),length(x))
par(pty='s')

image(x, y, kdata.mat, xlim=lim, ylim=lim, zlim=zlim, col=colour)
points(xy[sub.pos[1:25],1],xy[sub.pos[1:25],2])

    #conditional simulation n=25
cdata <- CondSimu(krige.method, x=xy[,1], y=xy[,2],  grid=FALSE,
                  model=model, param=param,
                  given=xy[sub.pos[1:25],],data=data[sub.pos[1:25]])

cdata.mat<-matrix(cdata,length(x),length(x))
par(pty='s')

image(x, y, cdata.mat, xlim=lim, ylim=lim, zlim=zlim, col=colour)
points(xy[sub.pos[1:25],1],xy[sub.pos[1:25],2])

    #conditional simulation n=50

cdata <- CondSimu(krige.method, x=xy[,1], y=xy[,2],  grid=FALSE,
                  model=model, param=param,given=xy[sub.pos[1:50],],data=data[sub.pos[1:50]])


cdata.mat<-matrix(cdata,length(x),length(x))
par(pty='s')

image(x, y, cdata.mat, xlim=lim, ylim=lim, zlim=zlim, col=colour)
points(xy[sub.pos[1:50],1],xy[sub.pos[1:50],2])

    #conditional simulation n=100
cdata <- CondSimu(krige.method, x=xy[,1], y=xy[,2],  grid=FALSE,
                  model=model, param=param,given=xy[sub.pos[1:100],],data=data[sub.pos[1:100]])

cdata.mat<-matrix(cdata,length(x),length(x))
par(pty='s')

image(x, y, cdata.mat, xlim=lim, ylim=lim, zlim=zlim, col=colour)
points(xy[sub.pos[1:100],1],xy[sub.pos[1:100],2])















