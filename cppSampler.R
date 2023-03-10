# Kieran Lee - Dec 2022
# Legacy version of Collapsed Gibbs Sampler, distributions for sampling require
# updating.

#setwd('C:\\Users\\kiera\\Documents\\UoN\\Year 4\\MATH4001')
#library(coda)
#source('createSample.r')

timestamp()
mu0 <- 0
k0 <- 0.2
a <- 8
b <- 10
alpha <- 1
N <- 100

for (sample in 1:3) {#each (sample = 1:2, .packages = c('coda','gtools','nimble'), .combine=cbind) %dopar% {
  sampleClusts <- gaussian2Sample(N)
  
  # initialize parameters
  
  for (run in 1:1){
    mus <- matrix(NA,1,N)
    sigs <- matrix(NA,1,N)
    sigs[1] <- rinvgamma(1,a,b)
    mus[1] <- rnorm(1,mu0,(sigs[1]/k0)^0.5)
    labels <- c(1)
    clusters <- matrix(NA,11000,N)
    entropy <- matrix(0,3,11000)
    for (iter in 1:11000) {
      # take clusters for final solution on last iteration.
      # or single cluster initially
      if (iter>1){
        clusters[iter,] <- clusters[iter-1,]
      } else {
        clusters[1,] <- matrix(1,1,N)
      }
      # STEP 1: update clustering
      for (obs in 1:N) {
        # determine current allocation
        alloc <- matrix(0,1,length(labels))
        for (labelInd in 1:length(labels)) {
          # count current obs (not including the one of interest) currently 
          #   allocated to the labelInd cluster
          alloc[labelInd] <- sum(clusters[iter,-obs]==labels[labelInd])
        }
        # remove any empty clusters
        # (if current obs is in a singleton)
        if (sum(alloc==0)==1){
          remove <- which(alloc==0)
          mus[labels[remove]] <- NA
          sigs[labels[remove]] <- NA
          # remove empty cluster from label list (special handling required if 
          #   cluster is first or last in list)
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
        alloc <- matrix(0,1,length(labels))
        for (labelInd in 1:length(labels)) {
          # count current obs (not including the one of interest) currently 
          #   allocated to the labelInd cluster
          alloc[labelInd] <- sum(clusters[iter,-obs]==labels[labelInd])
        }
        # allocate to existing cluster if allocation probability is correct
        temp_probs <- c(apply(cbind(alloc[1,], mus[labels], sigs[labels]), 1,
                              function(x) x[1] / (alpha + N - 1) *
                                dnorm(sampleClusts$sample[obs], x[2], x[3]^0.5)), 
                        alpha / (alpha + N - 1) * dt((sampleClusts$sample[obs]-mu0)/(b*(1+1/k0)/a)^0.5,2*a))/(b*(1+1/k0)/a)^0.5
        clusters[iter, obs] <- sample(x = c(labels, which(is.na(mus))[1]),
                                      size = 1, prob = temp_probs)
        if(clusters[iter, obs] == which(is.na(mus))[1]){
          labels <- c(labels,clusters[iter, obs])
          sigs[clusters[iter, obs]] <- rinvgamma(1,a,b)
          mus[clusters[iter, obs]] <- rnorm(1,mu0,(sigs[clusters[iter, obs]]/k0)^0.5)
        }
      }
      # STEP 2: Update parameters
      for (labelInd in 1:length(labels)) {
        # count current obs (not including the one of interest) currently 
        #   allocated to the labelInd cluster
        n <- sum(clusters[iter,]==labels[labelInd])
        ybar <- sum(sampleClusts$sample[clusters[iter,]==labels[labelInd]])/n
        ybar2 <- sum(sampleClusts$sample[clusters[iter,]==labels[labelInd]]^2)/n
        sigs[labels[labelInd]] <- rinvgamma(1,a+n/2+1/2,b+0.5*n*
                                              (mus[labels[labelInd]]^2-2*
                                                 mus[labels[labelInd]]*ybar+ybar2)
                                               +0.5*k0*(mus[labels[labelInd]]-mu0)^2)
        mus[labels[labelInd]] <- rnorm(1,(k0*mu0+n*ybar)
                                       /(k0+n),
                                       (sigs[labels[labelInd]]/
                                          (k0+n))^0.5)
        entropy[1,iter] <- entropy[1,iter] + 1
        entropy[2,iter] <- entropy[2,iter] - (n/N)*log(n/N)
        entropy[3,iter] <- entropy[3,iter] - abs(mus[labels[labelInd]]/N)*log(abs(mus[labels[labelInd]])/N)
      }
    }
    # calculate posterior similarity matrix
    print(sprintf("Sample: %i; Iteration: %d",sample,run))
    #psm <- comp.psm(clusters[1001:11000,])
    prob <- matrix(0,1,N)
    for (cluster in 1:N) {
      if (!is.na(mus[cluster])) {
        print(sprintf("%i observations in %i, mean = %.2f, st.dev = %.2f",sum(clusters[11000,]==cluster),cluster,mus[cluster],sigs[cluster]^0.5))
        prob <- prob + (sum(clusters[11000,]==cluster)/N)*dnorm(sampleClusts$sample,mus[cluster],sigs[cluster]^0.5)
      }
    }
    # find optimal clustering and plot allocation
    #opt <- minVI(psm,cls.draw=clusters[1001:11000,],method='draws')
    #hist(sampleClusts$sample,freq=FALSE,breaks=30)
    #points(sampleClusts$sample,prob,col=clusters[11000,]+1)
    plot(1:11000,entropy[1,],type='l',lty=1,col=2,ylim=c(0,20))
    lapply(2:3,function(i) lines(1:11000,entropy[i,],lty=1,col=i+1))
    print(apply(entropy[,1001:11000],1,function(x) effectiveSize(as.vector(x))))
    # print(table(opt$cl, sampleClusts$clustLabels))
    # print(sum(apply(table(opt$cl, sampleClusts$clustLabels),2,max))/N)
    # print(psm[1:10,1:7])
  }
  #effectiveSize(as.vector(entropy))
}
# print(apply(entropy,1,function(x) effectiveSize(as.vector(x))))
# plot(1:11000,entropy[1,],type='l',lty=1,col=2,ylim=c(0,20))
# lapply(2:3,function(i) lines(1:11000,entropy[i,],lty=1,col=i+1))
timestamp()