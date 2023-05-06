# Kieran Lee - Mar 2023
# Will run the Gibbs sampler, Collapsed Gibbs sampler and MCMC-ABC algorithm.
# Compare's performance of the 3 approaches. Requires the initialization of 
# parallel clusters.

source('createSample.r')
source('gibbsNIG.r')
source('collapsedGibbsNIG.r')
source('VItwoPartitions.r')
source('MCMC_ABC.r')

N <- 1000
numOfSamples <- 6
# adjust
skip <- 0 #already processed samples

# Create samples
results <- c(results[1:skip],sapply(1:numOfSamples,function(samNum){list('res'=list('sample'=gaussian4Sample(N),'gibbs'='','cGibbs'='','ABC'=''))}))

# Run Gibbs
timestamp()
results <- c(results[1:skip],foreach(sampleNum = (skip+1):(numOfSamples+skip), .packages = c('coda','nimble','gtools','wiqid','mcclust.ext'), .combine=c) %dopar% {
  start <- Sys.time()
  results[sampleNum]$res$gibbs <- gibbsNIG(results[sampleNum]$res$sample)
  end <- Sys.time()
  results[sampleNum]$res$gibbs <- c(results[sampleNum]$res$gibbs,'time'=end-start)
  results[sampleNum]
})
timestamp()
save(results,file="compOne.Rdata")

# Run Collapsed Gibbs
skip <- 36
breakP <- 9 #number of runs of algorithm
numOfSamples <- 14
timestamp()
results <- c(results[1:skip],foreach(sampleNum = (skip+1):(skip+breakP), .packages = c('coda','nimble','gtools','wiqid','mcclust.ext'), .combine=c) %dopar% {
  start <- Sys.time()
  results[sampleNum]$res$cGibbs <- collapsedGibbsNIG(results[sampleNum]$res$sample)
  end <- Sys.time()
  results[sampleNum]$res$cGibbs <- c(results[sampleNum]$res$cGibbs,'time'=end-start)
  results[sampleNum]
},results[(skip+breakP+1):(skip+numOfSamples)])
timestamp()
save(results,file="compOne.Rdata")

# Reload data
load("compOne.Rdata")
# Create break point
breakP <- floor(numOfSamples/2)
# Run ABC-MCMC
timestamp()
results <- c(results[1:skip],foreach(sampleNum = (skip+1):(skip+breakP), .packages = c('coda','nimble','gtools','wiqid','mcclust.ext'), .combine=c) %dopar% {
  start <- Sys.time()
  results[sampleNum]$res$ABC <- MCMC_ABC_NIG(results[sampleNum]$res$sample)
  end <- Sys.time()
  results[sampleNum]$res$ABC <- c(results[sampleNum]$res$ABC,'time'=end-start)
  results[sampleNum]
},results[(skip+breakP+1):(numOfSamples+skip)])
timestamp()
save(results,file="compOne.Rdata")

# Part 2
load("compOne.Rdata")
breakP <- 0
skip <- 18
numOfSamples <- 12
timestamp()
results <- c(results[1:(skip+breakP)],foreach(sampleNum = (skip+breakP+1):(skip+breakP+numOfSamples), .packages = c('coda','nimble','gtools','wiqid','mcclust.ext'), .combine=c) %dopar% {
  start <- Sys.time()
  results[sampleNum]$res$ABC <- MCMC_ABC_NIG(results[sampleNum]$res$sample)
  end <- Sys.time()
  results[sampleNum]$res$ABC <- c(results[sampleNum]$res$ABC,'time'=end-start)
  results[sampleNum]
})
timestamp()
save(results,file="compOne.Rdata")

for (result in results) {
  print("Result")
  print(sprintf("Gibbs ESS: %f",effectiveSize(as.vector(result$gibbs$entropies))))
  print(sprintf("Gibbs VI: %f",VItwoPartitions(result$sample$clustLabels,result$gibbs$bestPartition$cl)))
  print(sprintf("Collapsed gibbs ESS: %f",effectiveSize(as.vector(result$cGibbs$entropies))))
  print(sprintf("Collapsed gibbs VI: %f",VItwoPartitions(result$sample$clustLabels,result$cGibbs$bestPartition$cl)))
  print(sprintf("MCMC-ABC ESS: %f",effectiveSize(as.vector(result$ABC$entropies))))
  print(sprintf("MCMC-ABC VI: %f",VItwoPartitions(result$sample$clustLabels,result$ABC$bestPartition$cl)))
}
