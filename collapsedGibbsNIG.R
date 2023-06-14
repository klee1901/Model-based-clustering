# Kieran Lee - Feb 2023
# Runs a collapsed Gibbs sampler to approximate a partition for supplied data.
# Assumes data follows an infinite mixture model. Will return optimal partition,
# trace of number of components and entropies.

# required to sample from invgamma
library(nimble)
# required to calculate ESS
library(coda)
# required for finding 'best' partition
# loads mcclust which enables calculation of the psm
library(mcclust.ext)
# required for location-scale tdist
library(wiqid)

collapsedGibbsNIG <- function(sample,hyperparameters = list(mu0=0,k0=0.2,a=8,b=10,alpha=1)) {
  # sample is a list containing values ($sample) and true partition ($clustLabels)
  # hyperparameters is a list containing: mu0, k0, a, b, alpha.
  # set length of chain and burn-in
  iterations <- 2000
  burnIn <- 500
  # determine number of observations
  N <- length(sample$sample)
  # set hyperparameters
  mu0 <- hyperparameters$mu0
  k0 <- hyperparameters$k0
  a <- hyperparameters$a
  b <- hyperparameters$b
  alpha <- hyperparameters$alpha
  # Create arrays to hold a maximum N mus and sigmas.
  mus <- matrix(NA,1,N)
  sigs <- matrix(NA,1,N)
  # Initialise using prior marginals
  sigs[1] <- rinvgamma(1,a,b)
  mus[1] <- rt2(1,mu0,(b/(a*k0))^0.5,2*a)
  # Create an array to track the current number of components, each element is a 
  # pointer to the mus and sigma arrays
  labels <- c(1)
  # partitions will store the allocation after each iteration
  partitions <- matrix(NA,iterations,N)
  # entropy will track both the entropy and no. of components
  entropy <- matrix(0,2,iterations)
  for (iter in 1:iterations) {
    # intialise partitions for this iteration using the final solution from the
    # last iteration, or single component initially
    if (iter>1){
      partitions[iter,] <- partitions[iter-1,]
    } else {
      partitions[1,] <- matrix(1,1,N)
    }
    # STEP 1: update partition
    for (obs in 1:N) {
      # determine current allocation
      alloc <- matrix(0,1,length(labels))
      for (labelInd in 1:length(labels)) {
        # count current obs (not including the one of interest) currently 
        #   allocated to the labelInd component
        alloc[labelInd] <- sum(partitions[iter,-obs]==labels[labelInd])
      }
      # remove any empty components
      # (if current obs is in a singleton)
      if (sum(alloc==0)==1){
        remove <- which(alloc==0)
        mus[labels[remove]] <- NA
        sigs[labels[remove]] <- NA
        # remove empty component from label list (special handling required if 
        #   component is first or last in list)
        if (remove > 1){
          if (remove < length(labels)){
            labels <- c(labels[1:(remove-1)],labels[(remove+1):length(labels)])
          } else {
            labels <- labels[1:(length(labels)-1)]
          }
        } else {
          labels <- labels[2:length(labels)]
        }
      }
      # determine current component allocation
      # (recalculated in case a component has been destroyed)
      alloc <- matrix(0,1,length(labels))
      for (labelInd in 1:length(labels)) {
        # count current obs (not including the one of interest) currently 
        #   allocated to the labelInd component
        alloc[labelInd] <- sum(partitions[iter,-obs]==labels[labelInd])
      }
      # calculate allocation probabilities
      temp_probs <- c(apply(cbind(alloc[1,], mus[labels], sigs[labels]), 1,
                            function(x) x[1] / (alpha + N - 1) *
                              dnorm(sample$sample[obs], x[2], x[3]^0.5)), 
                      alpha / (alpha + N - 1) * dt2(sample$sample[obs],mu0,
                                                    (b*(1+1/k0)/a)^0.5,2*a))
      # sample an allocation for the observation
      partitions[iter, obs] <- sample(x = c(labels, which(is.na(mus))[1]),
                                    size = 1, prob = temp_probs)
      # if we choose to assign the obs to a new component
      if(partitions[iter, obs] == which(is.na(mus))[1]){
        # add the appropriate pointer to the labels list (the lowest value not
        # currently used as a pointer)
        labels <- c(labels,partitions[iter, obs])
        # initialise sigma and mu for the new component based using the
        # posterior marginals based off 1 observation (the current obs)
        sigs[partitions[iter, obs]] <- rinvgamma(1,a+0.5,b+0.5*k0*(mu0-sample$sample[obs])^2/(k0+1))
        mus[partitions[iter, obs]] <- rt2(1,(k0*mu0+sample$sample[obs])/(k0+1),((2*b*(k0+1)+k0*(mu0-sample$sample[obs])^2)/((2*a+1)*(k0+1)^2))^0.5,2*a+1)
      }
    }
    # STEP 2: Update parameters
    for (labelInd in 1:length(labels)) {
      # count current obs (not including the one of interest) currently 
      #   allocated to the labelInd component
      n <- sum(partitions[iter,]==labels[labelInd])
      ybar <- sum(sample$sample[partitions[iter,]==labels[labelInd]])/n
      ybar2 <- sum(sample$sample[partitions[iter,]==labels[labelInd]]^2)/n
      # use the posterior full conditionals to update mu and sigma
      sigs[labels[labelInd]] <- rinvgamma(1,a+n/2+1/2,b+0.5*k0*
                                            (mus[labels[labelInd]]-mu0)^2+0.5*n*
                                            (mus[labels[labelInd]]^2
                                                 -2*mus[labels[labelInd]]*ybar
                                                        +ybar2))
      mus[labels[labelInd]] <- rnorm(1,(k0*mu0+n*ybar)
                                     /(k0+n),
                                     (sigs[labels[labelInd]]/
                                        (k0+n))^0.5)
    }
    temp_tb <- table(partitions[iter, ]) / N
    entropy[1,iter] <- length(labels)
    entropy[2,iter] <- -sum(temp_tb * log(temp_tb))
  }
  psm <- comp.psm(partitions[(burnIn+1):iterations,])
  # find optimal partition
  opt <- minVI(psm,cls.draw=partitions[(burnIn+1):iterations,],method='draws')
  output <- list("bestPartition"=opt,"components"=entropy[1,(burnIn+1):iterations],"entropies"=entropy[2,(burnIn+1):iterations])
}