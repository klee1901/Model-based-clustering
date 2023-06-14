# Kieran Lee - May 2023
# Implements a Gibbs sampler for a finite mixture model with supposed G
# components.

# required to sample from dirichlet
library(gtools)
# required to sample from invgamma
library(nimble)
# required for location-scale tdist
library(wiqid)
# required to calculate ESS
library(coda)
# required for finding 'best' partition
# loads mcclust which enables calculation of the psm
library(mcclust.ext)

gibbsNIG <- function(sample,hyperparameters = list(mu0=0,k0=0.2,a=8,b=10,alpha=1)) {
  # sample is a list containing values ($sample) and true partition ($clustLabels)
  # hyperparameters is a list containing: mu0, k0, a, b, alpha.
  # set length of chain and burn-in
  iterations <- 2000
  burnIn <- 500
  # calculate true number of components
  G <- length(table(sample$clustLabels))
  # label components
  labs <- 1:G
  # determine number of observations
  N <- length(sample$sample)
  # set hyperparameters
  mu0 <- hyperparameters$mu0
  k0 <- hyperparameters$k0
  a <- hyperparameters$a
  b <- hyperparameters$b
  alpha <- rep(hyperparameters$alpha,G)
  # initialize parameters
  weights <- rdirichlet(1,alpha)
  # use prior marginals
  sigs <- rinvgamma(G,a,b)
  mus <- rt2(G,mu0,(b/(a*k0))^0.5,2*a)
  # set up matrix to store partition generated at each iteration (after burn-in)
  partitions <- matrix(NA,iterations-burnIn,N)
  # set up matrix to store partitions entropy at each iteration
  entropies <- matrix(0,1,iterations-burnIn)
  # iterate a number of times algorithm 5.1
  for (i in 1:iterations) {
    zs <- matrix(1,1,N)
    prob <- runif(N)
    # for each observation
    for (ind in 1:N) {
      # calculate likelihood it belongs to each component
      liks <- apply(cbind(weights[1,],mus,sigs),1,function(x) {x[1]*dnorm(sample$sample[ind],x[2],x[3]^0.5)})
      # sample an allocated with probability proportional to the likelihood
      # calculated above
      zs[ind] <- sample(labs,1,prob=liks)
    }
    # compute n, ybar for each group
    ns <- sapply(labs,function(lab) {sum(zs==lab)})
    ybars <- sapply(labs,function(lab) {sum(sample$sample[zs==lab])/ns[lab]})
    ybars2 <- sapply(labs,function(lab) {sum(sample$sample[zs==lab]^2)/ns[lab]})
    # For each group
    for (g in 1:G) {
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
    # update weight parameter using conditionals
    weights <- rdirichlet(1,alpha+ns)
    # discount first 500 iterations as 'burn in'
    if (i>burnIn){
      partitions[i-burnIn,] <- zs
      temp_tb <- table(zs) / N
      entropies[i-burnIn] <- -sum(temp_tb * log(temp_tb))
    }
  }
  # calculate posterior similarity matrix
  psm <- comp.psm(partitions)
  # find optimal partition, depending on VI distance
  opt <- minVI(psm,cls.draw=partitions,method='draws')
  list('bestPartition'=opt,'entropies'=entropies)
}