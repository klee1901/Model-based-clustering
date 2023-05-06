# Kieran Lee - Feb 2022
# Contains code to plot the output of Gibbs and Collapsed Gibbs samplers to
# enable comparison. Will visualize VI distance to true partition and Effective
# Sample Size (ESS). Can plot sample by sample or boxplot overview.

source('VItwoPartitions.r')
library(coda)

n <- length(results)
vis <- data.frame(rep(NA,3*n),rep(NA,3*n),c(rep('gibbs',n),rep('collapsedGibbs',n),rep('ABC-MCMC',n)),c(1:n,1:n,1:n),rep(NA,3*n))
colnames(vis) <- c('dist','ESS','type','sample','time')
vis$sample <- factor(vis$sample, levels=1:n)

i <- 1
for (result in results) {
  vis$dist[i] <- VItwoPartitions(result$sample$clustLabels,result$gibbs$bestPartition$cl)
  vis$ESS[i] <- effectiveSize(as.vector(result$gibbs$entropies))
  # times 60 to convert to secs, times 0.75 to account for 25% burnin
  vis$time[i] <- 45*result$gibbs$time/length(result$gibbs$entropies)
  vis$dist[i+n] <- VItwoPartitions(result$sample$clustLabels,result$cGibbs$bestPartition$cl)
  vis$ESS[i+n] <- effectiveSize(as.vector(result$cGibbs$entropies))
  vis$time[i+n] <- 45*result$cGibbs$time/length(result$cGibbs$entropies)
  vis$dist[i+2*n] <- VItwoPartitions(result$sample$clustLabels,result$ABC$bestPartition$cl)
  vis$ESS[i+2*n] <- effectiveSize(as.vector(result$ABC$entropies))
  vis$time[i+2*n] <- 45*result$ABC$time/length(result$ABC$entropies)
  i <- i+1
}
vis$ESS <- log(vis$ESS)

library(ggplot2)
library(dplyr)
# show 5 at a time
for (start in seq(1,n,5)) {
  subData <- filter(vis,sample %in% c(start:(start+4)))
  ggplot(subData,aes(sample))+geom_bar(aes(weight=dist,fill=type),position='dodge')+labs(y='VI distance from true partition',title='Comparison of gibbs and collapsed gibbs samplers \n in approximating the true partition of 5 random samples')
}
#ggplot(subData,aes(sample,ESS))+geom_point(aes(color=type))+labs(y='log(ESS)',title='Comparison of Effective Sample Size (ESS) \n on 5 random samples')

# summarise
ggplot(vis,aes(type,dist))+geom_boxplot(aes(color=type))+labs(x='algorithm',y='VI distance',title='Comparison of VI distances')+theme(legend.position='none')
ggplot(vis,aes(type,ESS))+geom_boxplot(aes(color=type))+labs(x='Sampler type',y='log(ESS)',title='Comparison of ESS')+theme(legend.position='none')
ggplot(vis,aes(type,time))+geom_boxplot(aes(color=type))+labs(x='Sampler type',y='Average time (secs) per iteration',title='Comparison of time per iteration')+theme(legend.position='none')

# rug plot
s <- gaussian2Sample(2000,c(-0.5,0.5))
out <- gibbsNIG(s)
dt <- data.frame(s$clustLabels,s$sample[1,],out$bestPartition$cl)
names(dt) <- c('true','vals','pred')
dt$true <- factor(dt$true)
dt$pred <- factor(dt$pred)
ggplot(dt)+geom_point(aes(vals,true,col=pred),shape=4,size=5)+labs(x='Value',y='True component (used to generate)',title='Comparision of components used to generate data\nand component allocation in optimal partition',color='Component\nallocated')
#geom_jitter