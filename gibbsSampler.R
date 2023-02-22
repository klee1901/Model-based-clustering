# required to sample from dirichlet
library(gtools)
# required to sample from invgamma
library(nimble)
# required for finding 'best' cluster
# loads mcclust which enables calculation of the psm
library(mcclust.ext)
# required to calculate ESS
library(coda)

for (rand in 1:3) {
  N <- 1000
  iterations <- 2000
  prob <- runif(N)
  sample <- matrix(NA,1,N)
  # Draw sample for selected dist (according to prob)
  for (i in 1:N) {
    if (prob[i] < 0.5) {
      sample[i] <- rnorm(1,2,1)
    } else {
      sample[i] <- rnorm(1,-2,1)
    }
  }
  trueCl <- ifelse(prob<0.5,2,1)
  hist(sample,breaks=30)
  #var(sample)
  
  # set hyperparameters
  mu0 <- 0
  s0 <- 4
  a <- 8
  b <- 10
  alpha <- c(3,2)
  # initialize parameters
  
  for (run in 1:1){
    weights <- rdirichlet(1,alpha)
    mus <- rnorm(2,mu0,s0)
    sigs <- rinvgamma(2,a,b)
    # iterate a number of times algorithm 5.1
    clusters <- matrix(NA,iterations-1000,N)
    entropy <- matrix(0,1,iterations-1000)
    for (i in 1:iterations) {
      #print(mus)
      zs <- matrix(1,1,N)
      prob <- runif(N)
      # for each observation
      for (ind in 1:N) {
        # calculate likelihood it came from each group
        g1lik <- weights[1]*dnorm(sample[ind],mus[1],sigs[1]^0.5)#/(sigs[1]*sqrt(2*pi)) * exp(-0.5*(sample[ind]-mus[1])^2*sigs[1]^(-2))
        g2lik <- weights[2]*dnorm(sample[ind],mus[2],sigs[2]^0.5)#/(sigs[2]*sqrt(2*pi)) * exp(-0.5*(sample[ind]-mus[2])^2*sigs[2]^(-2))
        if (prob[ind] > g1lik/(g1lik+g2lik)) {
          zs[ind] <- 2
        }
      }
      # compute n, ybar for both groups
      ns <- c(sum(zs==1),sum(zs==2))
      ybars <- c(sum(sample[zs==1])/ns[1],sum(sample[zs==2])/ns[2])
      ybars2 <- c(sum(sample[zs==1]^2)/ns[1],sum(sample[zs==2]^2)/ns[2])
      # For each group
      for (g in 1:2) {
        mus[g] <- rnorm(1,(mu0*sigs[g]+ns[g]*ybars[g]*s0^2)/(sigs[g]+s0^2*ns[g]),(sigs[g]*s0^2/(sigs[g]+s0^2*ns[g]))^0.5)
        sigs[g] <- rinvgamma(1,a+ns[g]/2,b+0.5*ns[g]*(mus[g]^2-2*mus[g]*ybars[g]+ybars2[g]))
      }
      weights <- rdirichlet(1,alpha+ns)
      # discount first 1000 iterations as 'burn in'
      if (i>1000){
        clusters[i-1000,] <- zs
        for (g in 1:2) {
          if (ns[g] != 0) {
            entropy[i-1000] <- entropy[i-1000]-(ns[g]/N)*log(ns[g]/N)
          }
        }
      }
    }
    # calculate posterier similarity matrix
    psm <- comp.psm(clusters)
    print(sprintf("Sample: %i; Iteration: %d",rand,run))
    print(sprintf("C1: weight = %.2f; mean = %.2f; st.dev = %.2f",weights[1],mus[1],sigs[1]^0.5))
    print(sprintf("C2: weight = %.2f; mean = %.2f; st.dev = %.2f",weights[2],mus[2],sigs[2]^0.5))
    print(sprintf("Entropy: %.2f",effectiveSize(as.vector(entropy))))
    # find optimal clustering and plot allocation
    opt <- minVI(psm,cls.draw=clusters,method='draws')
    print(sprintf("Min VI: %.2f",opt$value))
    hist(sample,freq=FALSE,breaks=30)
    points(sample,weights[1]*dnorm(sample,mus[1],sigs[1]^0.5)+weights[2]*dnorm(sample,mus[2],sigs[2]^0.5),col=clusters[iterations-1000,]+1)
    curve(0.5*dnorm(x,-2,1)+0.5*dnorm(x,2,1),add=TRUE,col='blue',lty=3,lwd=2)
    print(min(sum(apply(table(opt$cl, trueCl),2,max)),sum(apply(table(opt$cl, trueCl),1,max)))/N)
  }
}
# # find index of optimal cluster
# for (i in 1:10000) {
#   if (sum(clusters[i,]==opt$cl)==N) {
#     print(i)
#   }
# }
# Need to fix clusters by mean (label switching)
# Is it reasonable to base VI off whole cluster