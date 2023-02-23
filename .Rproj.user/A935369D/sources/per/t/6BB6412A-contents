# Kieran Lee - Feb 2023
# Will run the Gibbs sampler and Collapsed Gibbs sampler, comparing performance.
# Requires the initialization of parallel clusters.

source('createSample.r')
source('gibbsNIG.r')
source('collapsedGibbsNIG.r')
source('VItwoPartitions.r')

N <- 1000

timestamp()
results <- foreach(sample = 1:3, .packages = c('coda','nimble','gtools','wiqid','mcclust.ext'), .combine=c) %dopar% {
  sampleClusts <- gaussian2Sample(N)
  outGibbs <- gibbsNIG(sampleClusts)
  outCGibbs <- collapsedGibbsNIG(sampleClusts)#,list(mu0=0,k0=0.2,a=8,b=10,alpha=1))
  output <- list('out'=list("sample"=sampleClusts,"gibbs"=outGibbs,"cGibbs"=outCGibbs))
  names(output) <- c(paste('out',as.character(sample),sep=''))
  output
}
save(results,file="compT.Rdata")

for (result in results) {
  print("Result")
  print(sprintf("Gibbs ESS: %f",effectiveSize(as.vector(result$gibbs$entropies))))
  print(sprintf("Gibbs VI: %f",VItwoPartitions(result$sample$clustLabels,result$gibbs$bestPartition$cl)))
  # print(sprintf("Gibbs VI: %f",result$gibbs$bestPartition$value))
  print(sprintf("Collapsed gibbs ESS: %f",effectiveSize(as.vector(result$cGibbs$entropies))))
  print(sprintf("Collapsed gibbs VI: %f",VItwoPartitions(result$sample$clustLabels,result$cGibbs$bestPartition$cl)))
  # print(sprintf("Collapsed gibbs VI: %f",result$cGibbs$bestPartition$value))
  # print(sprintf("Collapsed gibbs ESS (Clusters): %f",effectiveSize(as.vector(result$cGibbs$clusters))))
  # print(sprintf("Collapsed gibbs ESS (Avg entropies): %f",effectiveSize(as.vector(result$cGibbs$entropies/result$cGibbs$clusters))))
}
timestamp()