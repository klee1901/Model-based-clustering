# Kieran Lee - Oct 2022
# Generates sample from an example mixture model and visulises.

# Store the dist means and sdevs
msdev <- matrix(c(2,3,3.4,1,0.5,1.3),3,2)
# Generate probabilities for selecting dist
prob <- runif(1000)
sample <- matrix(NA,1,1000)
# Draw sample for selected dist (according to prob)
for (i in 1:1000) {
  if (prob[i] < 0.3) {
    sample[i] <- rnorm(1,msdev[1,1],msdev[1,2])
  } else if (prob[i] < 0.8) {
    sample[i] <- rnorm(1,msdev[2,1],msdev[2,2])
  } else {
    sample[i] <- rnorm(1,msdev[3,1],msdev[3,2])
  }
}
# Plot!
hist(sample,freq=FALSE,ylim=c(0,1),breaks = 30)
curve(0.3*dnorm(x,msdev[1,1],msdev[1,2]),add=TRUE,col='red')
curve(0.5*dnorm(x,msdev[2,1],msdev[2,2]),add=TRUE,col='red')
curve(0.2*dnorm(x,msdev[3,1],msdev[3,2]),add=TRUE,col='red')

