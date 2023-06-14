# Kieran Lee - Mar 2023
# Runs MCMC ABC algorithm to generate an approximate partition for supplied
# data. Will return optimal partition, trace of number of clusters and 
# entropies.

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
        #temp_probs <- c(alloc[1,], alpha)
        newClusters[obs] <- sample(x = c(labels, which(is.na(mus))[1]),
                                   size = 1, prob = c(alloc[1,], alpha))
        # if we choose to generate a new cluster
        if(newClusters[obs] == which(is.na(mus))[1]){
          # add the appropriate pointer to the labels list (the lowest value not
          # currently used as a pointer)
          labels <- c(labels,newClusters[obs])
          # initialise sigma and mu for the new cluster based using the prior 
          # marginals
          sigs[newClusters[obs]] <- rinvgamma(1,a,b)
          mus[newClusters[obs]] <- rt2(1,mu0,(b/(a*k0))^0.5,2*a)
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
    for (obs in 1:N) {
      # determine rank of current observation
      rank <- which(sort(sample$sample)==sample$sample[obs])[1]
      # determine position of synthetic obs with same rank
      synPlace <- which(syntheticData==sort(syntheticData)[rank])[1]
      # allocate obs to component of paired synthetic obs
      clusters[iter+1,obs] <- newClusters[synPlace]
    }
    #clusters[iter+1,] <- newClusters
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
    temp_tb <- table(clusters[iter+1, ]) / N
    entropy[1,iter] <- length(labels)
    entropy[2,iter] <- -sum(temp_tb * log(temp_tb))
  }
  psm <- comp.psm(clusters[(burnIn+2):(iterations+1),])
  # find optimal clustering
  opt <- minVI(psm,cls.draw=clusters[(burnIn+2):(iterations+1),],method='draws')
  output <- list("bestPartition"=opt,"clusters"=entropy[1,(burnIn+1):iterations],"entropies"=entropy[2,(burnIn+1):iterations],"mc"=clusters[(burnIn+2):(iterations+1),])
}

MCMC_ABC_NIG_rep <- function(sample, tol = 1.7, hyperparameters = list(mu0=0,k0=0.2,a=8,b=10,alpha=1)) {
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
  syntheticData <- matrix(NA,1,N)
  iter <- 1
  while (is.na(clusters[iterations+1,1])) {
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
      #temp_probs <- c(alloc[1,], alpha)
      newClusters[obs] <- sample(x = c(labels, which(is.na(mus))[1]),
                                 size = 1, prob = c(alloc[1,], alpha))
      # if we choose to generate a new cluster
      if(newClusters[obs] == which(is.na(mus))[1]){
        # add the appropriate pointer to the labels list (the lowest value not
        # currently used as a pointer)
        labels <- c(labels,newClusters[obs])
        # initialise sigma and mu for the new cluster based using the prior 
        # marginals
        sigs[newClusters[obs]] <- rinvgamma(1,a,b)
        mus[newClusters[obs]] <- rt2(1,mu0,(b/(a*k0))^0.5,2*a)
      }
    }
    
    # Generate synthetic data
    for (obs in 1:N) {
      syntheticData[obs] <- rnorm(1,mus[newClusters[obs]],sigs[newClusters[obs]]^0.5)
    }
    # Calculate distance between synthetic and actual data
    dist <- sum(((sort(syntheticData)-sort(sample$sample))^2)/N)^0.5
    
    # STEP 2: Update parameters
    if (dist<tol) {
      # Set clusters allocation
      for (obs in 1:N) {
        # determine rank of current observation
        rank <- which(sort(sample$sample)==sample$sample[obs])[1]
        # determine position of synthetic obs with same rank
        synPlace <- which(syntheticData==sort(syntheticData)[rank])[1]
        # allocate obs to component of paired synthetic obs
        clusters[iter+1,obs] <- newClusters[synPlace]
      }
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
      temp_tb <- table(clusters[iter+1, ]) / N
      entropy[1,iter] <- length(labels)
      entropy[2,iter] <- -sum(temp_tb * log(temp_tb))
      iter <- iter + 1
    } else if (iter==1) {
      n <- length(labels)
      sigs[labels] <- rinvgamma(n,a,b)
      mus[labels] <- rt2(n,mu0,(b/(a*k0))^0.5,2*a)
    } else {
      for (labelInd in 1:length(labels)) {
        # count current obs (not including the one of interest) currently 
        #   allocated to the labelInd cluster
        n <- sum(clusters[iter,]==labels[labelInd])
        if (n>0) {
          ybar <- sum(sample$sample[clusters[iter,]==labels[labelInd]])/n
          ybar2 <- sum(sample$sample[clusters[iter,]==labels[labelInd]]^2)/n
          # use the posterior full conditionals to update mu and sigma
          sigs[labels[labelInd]] <- rinvgamma(1,a+n/2+1/2,b+0.5*k0*
                                                (mus[labels[labelInd]]-mu0)^2+0.5*n*
                                                (mus[labels[labelInd]]^2
                                                 -2*mus[labels[labelInd]]*ybar
                                                 +ybar2))
          mus[labels[labelInd]] <- rnorm(1,(k0*mu0+n*ybar)/(k0+n),
                                         (sigs[labels[labelInd]]/(k0+n))^0.5)
        }
      }
    }
  }
  psm <- comp.psm(clusters[(burnIn+2):(iterations+1),])
  # find optimal clustering
  opt <- minVI(psm,cls.draw=clusters[(burnIn+2):(iterations+1),],method='draws')
  output <- list("bestPartition"=opt,"clusters"=entropy[1,(burnIn+1):iterations],"entropies"=entropy[2,(burnIn+1):iterations],"mc"=clusters[(burnIn+2):(iterations+1),])
}

MCMC_ABC_NIG_1percent <- function(sample, hyperparameters = list(mu0=0,k0=0.2,a=8,b=10,alpha=1)) {
  # sample is a list containing values ($sample) and true partition ($clustLabels)
  # hyperparameters is a list containing: mu0, k0, a, b, alpha.
  # set length of chain and burn-in
  iterations <- 4000
  burnIn <- 0
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
  dists <- matrix(NA,1,iterations+1)
  clusters[1,] <- sample$clustLabels
  # entropy will track both the entropy and no. of clusters
  entropy <- matrix(0,2,iterations)
  #print(c(mus[1],sigs[1]))
  for (iter in 1:iterations) {
    syntheticData <- matrix(NA,1,N)
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
      #temp_probs <- c(alloc[1,], alpha)
      newClusters[obs] <- sample(x = c(labels, which(is.na(mus))[1]),
                                 size = 1, prob = c(alloc[1,], alpha))
      # if we choose to generate a new cluster
      if(newClusters[obs] == which(is.na(mus))[1]){
        # add the appropriate pointer to the labels list (the lowest value not
        # currently used as a pointer)
        labels <- c(labels,newClusters[obs])
        # initialise sigma and mu for the new cluster based using the prior 
        # marginals
        sigs[newClusters[obs]] <- rinvgamma(1,a,b)
        mus[newClusters[obs]] <- rt2(1,mu0,(b/(a*k0))^0.5,2*a)
      }
    }
      
    # Generate synthetic data
    for (obs in 1:N) {
      syntheticData[obs] <- rnorm(1,mus[newClusters[obs]],sigs[newClusters[obs]]^0.5)
    }
    # Calculate distance between synthetic and actual data
    dists[iter+1] <- sum(((sort(syntheticData)-sort(sample$sample))^2)/N)^0.5
    
    # STEP 2: Update parameters
    # Set clusters allocation
    for (obs in 1:N) {
      # determine rank of current observation
      rank <- which(sort(sample$sample)==sample$sample[obs])[1]
      # determine position of synthetic obs with same rank
      synPlace <- which(syntheticData==sort(syntheticData)[rank])[1]
      # allocate obs to component of paired synthetic obs
      clusters[iter+1,obs] <- newClusters[synPlace]
    }
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
    temp_tb <- table(clusters[iter+1, ]) / N
    entropy[1,iter] <- length(labels)
    entropy[2,iter] <- -sum(temp_tb * log(temp_tb))
  }
  dists[1] <- max(dists[2:(iterations+1)])
  qrt <- quantile(dists[2:(iterations+1)],probs=0.5)
  partitionsInd <- which(dists < qrt)
  MCout <- clusters[partitionsInd,]
  psm <- comp.psm(MCout)
  # find optimal clustering
  opt <- minVI(psm,cls.draw=MCout,method='draws')
  output <- list("bestPartition"=opt,"clusters"=entropy[1,(burnIn+1):iterations],"entropies"=entropy[2,(burnIn+1):iterations])
}

adaptive_MCMC_ABC_NIG <- function(sample, tol = 1.7, hyperparameters = list(mu0=0,k0=0.2,a=8,b=10,alpha=1)) {
  # sample is a list containing values ($sample) and true partition ($clustLabels)
  # hyperparameters is a list containing: mu0, k0, a, b, alpha.
  # set length of chain and burn-in
  iterations <- 2000
  burnIn <- 500
  # determine number of observations
  N <- length(sample$sample)
  # calculate initial and limit values for the tolerance
  tol0 <- 2#(N*log(N))**0.5
  tolStar <- 0.1
  tols <- c()
  tol <- tol0
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
    if (iter%%200==0) {print('iter')}
    syntheticData <- matrix(NA,1,N)
    # Update tolerance
    if (iter>1) {tol <- exp(-(iter-min(burnIn,iter))/1000)*tol0+(1-exp(-(iter-min(burnIn,iter))/1000))*exp(-(iter-min(burnIn,iter))/1000)*quantile(tols,probs=0.1)+(1-exp(-(iter-min(burnIn,iter))/1000)-(1-exp(-(iter-min(burnIn,iter))/1000))*exp(-(iter-min(burnIn,iter))/1000))*tolStar}
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
        #temp_probs <- c(alloc[1,], alpha)
        newClusters[obs] <- sample(x = c(labels, which(is.na(mus))[1]),
                                   size = 1, prob = c(alloc[1,], alpha))
        # if we choose to generate a new cluster
        if(newClusters[obs] == which(is.na(mus))[1]){
          # add the appropriate pointer to the labels list (the lowest value not
          # currently used as a pointer)
          labels <- c(labels,newClusters[obs])
          # initialise sigma and mu for the new cluster based using the prior 
          # marginals
          sigs[newClusters[obs]] <- rinvgamma(1,a,b)
          mus[newClusters[obs]] <- rt2(1,mu0,(b/(a*k0))^0.5,2*a)
        }
      }
      
      # Generate synthetic data
      for (obs in 1:N) {
        syntheticData[obs] <- rnorm(1,mus[newClusters[obs]],sigs[newClusters[obs]]^0.5)
      }
      # Update tolerance?
      # Calculate distance between synthetic and actual data
      dist <- sum(((sort(syntheticData)-sort(sample$sample))^2)/N)^0.5
      #print(length(table(newClusters)))
      if (iter%%200==0) {print(tol)}
    }
    # Update list of distances
    if (length(tols) < 100) {
      tols <- c(tols, dist)
    } else {
      tols <- c(tols[2:100],dist)
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
