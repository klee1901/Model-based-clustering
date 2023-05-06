# Kieran Lee - Feb 2023
# Performs k-means clustering for some 2D data from a Gaussian mixture model,
# also creates plot to indicate shape of clusters.

# data generation
N <- 300
xs <- matrix(NA,1,N)
ys <- matrix(NA,1,N)
xs[1:100] <- rnorm(100,-2,0.8)
ys[1:100] <- rnorm(100,-2,0.8)
xs[101:200] <- rnorm(100,2,0.8)
ys[101:200] <- rnorm(100,-1,0.3)
xs[201:300] <- rnorm(100,1,0.2)
ys[201:300] <- rnorm(100,1,0.2)
cl <- c(matrix(1,1,100),matrix(2,1,100),matrix(3,1,100))
plot(xs,ys,col=cl+1)
#kmeans
init <- sample(1:100,3)
muxs <- xs[init]
muys <- ys[init]
for (i in 1:10) {
  alloc <- matrix(1,1,300)
  for (j in 1:300) {
    dist <- (xs[j]-muxs[1])^2+(ys[j]-muys[1])^2
    for (k in 2:3) {
      if ((xs[j]-muxs[k])^2+(ys[j]-muys[k])^2<dist) {
        alloc[j] <- k
        if (k==2) {
          dist <- (xs[j]-muxs[k])^2+(ys[j]-muys[k])^2
        }
      }
    }
  }
  for (k in 1:3) {
    if (!is.na(sum(alloc==k))){
      muxs[k] <- sum(xs[alloc==k])/sum(alloc==k)
      muys[k] <- sum(ys[alloc==k])/sum(alloc==k)
    }
  }
}
# plot(xs,ys,col=alloc+1)
# points(muxs,muys,pch=4)
# abline(0.5*(muys[2]^2-muys[1]^2+muxs[2]^2-muxs[1]^2)/(muys[2]-muys[1]),(muxs[1]-muxs[2])/(muys[2]-muys[1]))
# abline(0.5*(muys[3]^2-muys[1]^2+muxs[3]^2-muxs[1]^2)/(muys[3]-muys[1]),(muxs[1]-muxs[3])/(muys[3]-muys[1]))
# abline(0.5*(muys[2]^2-muys[3]^2+muxs[2]^2-muxs[3]^2)/(muys[2]-muys[3]),(muxs[3]-muxs[2])/(muys[2]-muys[3]))
xgrid <- seq(-5,5,0.025)
ygrid <- seq(-5,2,0.025)
grid <- matrix(1,length(xgrid),length(ygrid))
for (i in 1:length(xgrid)) {
  for (j in 1:length(ygrid)) {
    dist <- (xgrid[i]-muxs[1])^2+(ygrid[j]-muys[1])^2
    for (k in 2:3) {
      if ((xgrid[i]-muxs[k])^2+(ygrid[j]-muys[k])^2<dist) {
        grid[i,j] <- k
        if (k==2) {
          dist <- (xgrid[i]-muxs[k])^2+(ygrid[j]-muys[k])^2
        }
      }
    }
  }
}
plot(xs,ys,pch=19)
for (i in 1:length(xgrid)) {
  points(matrix(xgrid[i],1,length(ygrid)),ygrid,col=grid[i,]+1,pch=19,cex=0.1)
}
# Plot confidence intervals
library(car)
plot(xs,ys,pch=19,ann=F)
ellipse(c(-2,-2),matrix(c(0.8,0,0,0.8),2,2)^2,1,col=3,lty=2)
ellipse(c(-2,-2),matrix(c(1.6,0,0,1.6),2,2)^2,1,col=3,lty=2)
ellipse(c(2,-1),matrix(c(0.8,0,0,0.3),2,2)^2,1,lty=2)
ellipse(c(2,-1),matrix(c(1.6,0,0,0.6),2,2)^2,1,lty=2)
ellipse(c(1,1),matrix(c(0.2,0,0,0.2),2,2)^2,1,col=2,lty=2)
ellipse(c(1,1),matrix(c(0.4,0,0,0.4),2,2)^2,1,col=2,lty=2)

## Connected clusters
N <- 300
xs <- matrix(NA,1,N)
ys <- matrix(NA,1,N)
xs[1:100]  <- runif(100,0,2)
xs[101:200]  <- runif(100,0,1)
ys[201:300]  <- runif(100,0,2)

# Cluster 1
for (i in 1:100) {
  ys[i] <- (xs[i]-1)^2
}
# Cluster 2
for (i in 101:200) {
  ys[i] <- 2-xs[i]
}
# Cluster 3
for (i in 201:300) {
  xs[i] <- 3-(ys[i]-1)^2
}
# Apply noise
ys <- ys + rnorm(300,0,0.1)
plot(xs,ys,axes=F,xlab='',ylab='')
for (k in 1:3) {
  points(xs[(100*(k-1)+1):(100*k)],ys[(100*(k-1)+1):(100*k)],col=k+1)
}
