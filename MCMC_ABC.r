# Kieran Lee - Mar 2023
# Runs MCMC ABC algorithm to generate an approximate partition for supplied
# data. Will return optimal partition, trace of number of clusters and 
# entropies, and Markov Chain.

# required to sample from invgamma
library(nimble)
# required to calculate ESS
library(coda)
# required for finding 'best' cluster
# loads mcclust which enables calculation of the psm
library(mcclust.ext)
# required for location-scale tdist
library(wiqid)

MCMC_ABC_NIG <- function(sample, tol = 1.7, hyperparameters = list(mu0=0,k0=0.2,a=8,b=10,alpha=1)) {
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
  # Create an array to track the current number of clusters, each element is a 
  # pointer to the mus and sigma arrays
  n <- length(table(sample$clustLabels))
  labels <- c(1:n)
  # Initialise using prior marginals
  sigs[1:n] <- rinvgamma(n,a,b)
  mus[1:n] <- rt2(n,mu0,(b/(a*k0))^0.5,2*a)
  # clusters will store the allocation after each iteration
  clusters <- matrix(NA,iterations+1,N)
  clusters[1,] <- sample$clustLabels
  # entropy will track both the entropy and no. of clusters
  entropy <- matrix(0,2,iterations)
  #print(c(mus[1],sigs[1]))
  for (iter in 1:iterations) {
    if (iter==1) {print('iter')}
    syntheticData <- matrix(NA,1,N)
    dist <- tol+1
    while (dist > tol) {
      #print(iter)
      # reset clusters (in case of previous iterations)
      # determine original allocation
      alloc <- matrix(0,1,length(labels))
      for (labelInd in 1:length(labels)) {
        # count current obs currently allocated to the labelInd cluster
        alloc[labelInd] <- sum(clusters[iter,]==labels[labelInd])
      }
      
      # remove any empty clusters (those that've been generated in previous
      # iterations)
      if (sum(alloc==0)>0){
        remove <- which(alloc==0)
        mus[labels[remove]] <- NA
        sigs[labels[remove]] <- NA
        # remove empty cluster from label list (special handling required if 
        #   cluster is first or last in list)
        labels <- labels[!(labels %in% labels[remove])]
      }
      
      # Generate cluster proposal
      newClusters <- matrix(0,1,N)
      for (obs in 1:N) {
        # determine current allocation
        alloc <- matrix(0,1,length(labels))
        for (labelInd in 1:length(labels)) {
          # count current obs (not including the one of interest) currently 
          #   allocated to the labelInd cluster
          alloc[labelInd] <- (sum(clusters[iter,]==labels[labelInd]) + 
                                sum(newClusters==labels[labelInd]))
        }
        # allocate to component using predictive density
        temp_probs <- c(apply(cbind(alloc[1,], mus[labels], sigs[labels]), 1,
                              function(x) x[1] / (alpha + N - 1) *
                                dnorm(sample$sample[obs], x[2], x[3]^0.5)), 
                        alpha / (alpha + N - 1) * dt2(sample$sample[obs],mu0,
                                                      (b*(1+1/k0)/a)^0.5,2*a))
        newClusters[obs] <- sample(x = c(labels, which(is.na(mus))[1]),
                                   size = 1, prob = temp_probs)
        # if we choose to generate a new cluster
        if(newClusters[obs] == which(is.na(mus))[1]){
          # add the appropriate pointer to the labels list (the lowest value not
          # currently used as a pointer)
          labels <- c(labels,newClusters[obs])
          # initialise sigma and mu for the new cluster based using the posterior 
          # marginals based off 1 observation (the current obs)
          sigs[newClusters[obs]] <- rinvgamma(1,a+0.5,b+0.5*k0*
                                                   (mu0-sample$sample[obs])^2/(k0+1))
          mus[newClusters[obs]] <- rt2(1,(k0*mu0+sample$sample[obs])/(k0+1),
                                          ((2*b*(k0+1)+k0*(mu0-sample$sample[obs])^2)/
                                             ((2*a+1)*(k0+1)^2))^0.5,2*a+1)
        }
      }
      
      # Generate synthetic data
      for (obs in 1:N) {
        syntheticData[obs] <- rnorm(1,mus[newClusters[obs]],sigs[newClusters[obs]]^0.5)
      }
      # Calculate distance between synthetic and actual data
      dist <- sum(((sort(syntheticData)-sort(sample$sample))^2)/N)^0.5
      #print(length(table(newClusters)))
      if (iter==1) {print(dist)}
    }
    
    
    # STEP 2: Update parameters
    # Set clusters allocation
    clusters[iter+1,] <- newClusters
    for (labelInd in 1:length(labels)) {
      # count current obs (not including the one of interest) currently 
      #   allocated to the labelInd cluster
      n <- sum(clusters[iter+1,]==labels[labelInd])
      ybar <- sum(sample$sample[clusters[iter+1,]==labels[labelInd]])/n
      ybar2 <- sum(sample$sample[clusters[iter+1,]==labels[labelInd]]^2)/n
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
    temp_tb <- table(clusters[iter, ]) / N
    entropy[1,iter] <- length(labels)
    entropy[2,iter] <- -sum(temp_tb * log(temp_tb))
  }
  psm <- comp.psm(clusters[(burnIn+2):(iterations+1),])
  # find optimal clustering
  opt <- minVI(psm,cls.draw=clusters[(burnIn+2):(iterations+1),],method='draws')
  output <- list("bestPartition"=opt,"clusters"=entropy[1,(burnIn+1):iterations],"entropies"=entropy[2,(burnIn+1):iterations])#,"mc"=clusters[(burnIn+2):(iterations+1),])
}
  