gaussian2Sample <- function (N=1000,mus=c(-2,2),sigs=c(1,1)){
  prob <- runif(N)
  sample <- matrix(NA,1,N)
  # Draw sample for selected dist (according to prob)
  for (i in 1:N) {
    if (prob[i] < 0.5) {
      sample[i] <- rnorm(1,mus[1],sigs[1])
    } else {
      sample[i] <- rnorm(1,mus[2],sigs[2])
    }
  }
  trueCl <- ifelse(prob<0.5,1,2)
  hist(sample,breaks=30)
  return(list("sample"=sample,"clustLabels"=trueCl))
}

gaussian4Sample <- function (N=1000,mus=c(-6,-2,2,6),sigs=rep(1,4)){
  prob <- runif(N)
  sample <- matrix(NA,1,N)
  trueCl <- ceiling(4*prob)
  # Draw sample for selected dist (according to prob)
  for (i in 1:N) {
    sample[i] <- rnorm(1,mus[trueCl[i]],sigs[trueCl[i]])
  }
  hist(sample,breaks=30)
  return(list("sample"=sample,"clustLabels"=trueCl))
}
