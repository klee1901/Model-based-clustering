# Kieran Lee - Feb 2023
# Copy of gibbsNIG to form potential multidimensional version.

# required to sample from dirichlet
library(gtools)
# required to sample from invgamma
library(nimble)
# required for location-scale tdist
library(wiqid)
# required to calculate ESS
library(coda)
# required for finding 'best' cluster
# loads mcclust which enables calculation of the psm
library(mcclust.ext)

gibbsNIG <- function(sample,hyperparameters = list(mu0=0,k0=0.2,a=8,b=10,alpha=c(3,2))) {
  # sample is a list containing values ($sample) and true partition ($clustLabels)
  # hyperparameters is a list containing: mu0, k0, a, b, alpha.
  # set length of chain and burn-in
  iterations <- 2000
  burnIn <- 200
  # determine number of observations
  N <- length(sample$sample)
  # set hyperparameters
  mu0 <- hyperparameters$mu0
  k0 <- hyperparameters$k0
  a <- hyperparameters$a
  b <- hyperparameters$b
  alpha <- hyperparameters$alpha
  # initialize parameters
  weights <- rdirichlet(1,alpha)
  # use prior marginals
  sigs <- riwish(v0,S)
  mus <- rt2(2,mu0,(b/(a*k0))^0.5,2*a)
  # iterate a number of times algorithm 5.1
  clusters <- matrix(NA,iterations-burnIn,N)
  entropies <- matrix(0,1,iterations-burnIn)
  for (i in 1:iterations) {
    zs <- matrix(1,1,N)
    prob <- runif(N)
    # for each observation
    for (ind in 1:N) {
      # calculate likelihood it came from each group
      g1lik <- weights[1]*dnorm(sample$sample[ind],mus[1],sigs[1]^0.5)
      g2lik <- weights[2]*dnorm(sample$sample[ind],mus[2],sigs[2]^0.5)
      zs[ind] <- sample(c(1,2),1,prob=c(g1lik,g2lik))
    }
    # compute n, ybar for both groups
    ns <- c(sum(zs==1),sum(zs==2))
    ybars <- c(sum(sample$sample[zs==1])/ns[1],sum(sample$sample[zs==2])/ns[2])
    ybars2 <- c(sum(sample$sample[zs==1]^2)/ns[1],sum(sample$sample[zs==2]^2)/ns[2])
    # For each group
    for (g in 1:2) {
      if (ns[g] != 0) {
        # update using posterior full conditionals
        sigs[g] <- rinvgamma(1,a+ns[g]/2+1/2,b+0.5*k0*(mus[g]-mu0)^2
                             +0.5*ns[g]*(mus[g]^2-2*mus[g]*ybars[g]+ybars2[g]))
        mus[g] <- rnorm(1,(k0*mu0+ns[g]*ybars[g])/(k0+ns[g]),
                        (sigs[g]/(k0+ns[g]))^0.5)
      } else {
        # re-initalise if no observations have been allocated to the group
        sigs[g] <- rinvgamma(1,a,b)
        mus[g] <- rt2(1,mu0,(b/(a*k0))^0.5,2*a)
      }
    }
    weights <- rdirichlet(1,alpha+ns)
    # discount first 1000 iterations as 'burn in'
    if (i>burnIn){
      clusters[i-burnIn,] <- zs
      temp_tb <- table(zs) / N
      entropies[i-burnIn] <- -sum(temp_tb * log(temp_tb))
      #for (g in 1:2) {
      #  if (ns[g] != 0) { 
      #    entropies[i-burnIn] <- entropies[i-burnIn]-(ns[g]/N)*log(ns[g]/N)
      #  }
      #}
    }
  }
  # calculate posterior similarity matrix
  psm <- comp.psm(clusters)
  # find optimal clustering
  opt <- minVI(psm,cls.draw=clusters,method='draws')
  list('bestPartition'=opt,'entropies'=entropies,'mc'=clusters)
}