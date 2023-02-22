library(doParallel)
# set-up
cl <- makeCluster(3)
registerDoParallel(cl)
#foreach (sample = 1:2, .packages = c('coda','gtools','nimble')) %dopar% {
# check
getDoParWorkers()
# end
stopCluster(cl)
