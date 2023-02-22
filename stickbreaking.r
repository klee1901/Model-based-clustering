alpha <- 1
n <- 100
reps <- 3
for (r in 1:reps) {
  vs <- rbeta(n,1,alpha)
  eta <- matrix(vs[1],1,n)
  for (i in 2:n) {
    eta[i] <- eta[i-1]*(1-vs[i-1])*vs[i]/vs[i-1]
  }
  esum <- matrix(0,1,n+1)
  esum[2] <- eta[1]
  for (i in 2:n) {
    esum[i+1] <- esum[i]+eta[i]
  }
  if (r==1) {
    plot(c(0,1),c(1,1),type='l',ylim=c(1,reps),main = "Stick breaking for alpha = 1",xlab=NA,ylab=NA,yaxt='n',col=7)
    points(esum,rep(1,n+1),pch=124,col=7)
  } else {
    points(c(0,1),c(r,r),type='l',col=7)
    points(esum,rep(r,n+1),pch=124,col=7)
  }
}
