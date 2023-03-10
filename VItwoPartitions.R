# Kieran Lee - Jan 2023
# Calculates VI distance between 2 partitions.

# section 3 Wade-Ghahramani
VItwoPartitions <- function(c1,c2) {
  N <- length(c1)
  nijs <- table(c1,c2)
  nis <- apply(nijs,1,sum)
  k1 <- length(nis)
  njs <- apply(nijs,2,sum)
  k2 <- length(njs)
  Hc1 <- -sum((nis/N)*log2(nis/N))
  Hc2 <- -sum((njs/N)*log2(njs/N))
  Ic1c2 <- 0
  for (i in 1:k1) {
    for (j in 1:k2) {
      if (nijs[i,j] != 0){
        Ic1c2 <- Ic1c2 + (nijs[i,j]/N)*log2((nijs[i,j]*N)/(nis[i]*njs[j]))
      }
    }
  }
  (Hc1 + Hc2 - 2*Ic1c2) / log2(N) 
}
