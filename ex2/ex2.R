### ex 2 ###


### preamble ###
library(xtable)
library("tidyverse")

### (a) ###
#Switch the default random number generator in R to Wichmann-Hill. 

#?Random

RNGkind('Wichmann-Hill')
set.seed(143, 'Wichmann-Hill')

#Simulate N = 1000 binomial random variables B(n = 10; p = 0:4) using three approaches: 

N <- 1000
n <- 10
p <- 0.4



#inversion method, 
#by simulating corresponding Bernoulli random variables 
#by inversion method
#by using R function rbinom. 



#Plot the empirical probability density functions of all three samples on one panel. 
#Comment on the results. 
#Switch the random number generator back to its default.