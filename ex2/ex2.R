  ### ex 2 ###


  ### preamble ###
library(xtable)
library("tidyverse")



  ### (a) ###
  #Switch the default random number generator in R to Wichmann-Hill. 

?Random

RNGkind('Wichmann-Hill')
set.seed(143, 'Wichmann-Hill')

#Simulate N = 1000 binomial random variables B(n = 10; p = 0:4) using three approaches: 

N <- 1000
n <- 10
p <- 0.4



##method 1: by inversion method, 
# Let us generate [0,1]-uniform r.v.
uni_rv <- runif(N)
binom_inversion <- function(X,n,p) {
  y= NULL
  for (k in 1:length(X)){
    x = X[k]
    for (l in 10:0) {if (x<=pbinom(l,n,p)) y[k]=l}
  }
  return(y)
}
bin1 <- binom_inversion(uni_rv,n,p)

## method 2: by simulating corresponding Bernoulli random variables by inversion method
# Let us produce n*N uniform r.v. and transform them
uni_rv2 <- runif(n*N)
Bern <- matrix(as.numeric(uni_rv2<p), nrow = N)
# then add up the independent bernoullis to get binomials
bin2 <- rowSums(Bern)

## method 3: by using R function rbinom
bin3 <- rbinom(N, n, p)




#Plot the empirical probability density functions of all three samples on one panel. 

#Comment on the results. 
#Switch the random number generator back to its default.