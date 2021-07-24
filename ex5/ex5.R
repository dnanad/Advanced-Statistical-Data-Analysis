####      Exercise 5      #####

### Preamble ####



library(tidyverse)
library("ggplot2")
library("dplyr")


#####





em_mixture <- function(params,X,clusters = 2,tol = .00001,maxits  = 100,showits = TRUE){
  
  # Arguments are starting parameters (means, covariances, cluster probability),
  # data, number of clusters desired, tolerance, maximum iterations, and whether
  # to show iterations
  
  # Starting points
  N     = nrow(X)
  nams  = names(params)
  mu    = params$mu
  var   = params$var
  probs = params$probs
  
  # Other initializations
  # initialize cluster 'responsibilities', i.e. probability of cluster
  # membership for each observation i
  ri = matrix(0, ncol = clusters, nrow = N) 
  it = 0
  converged = FALSE
  
  if (showits)                                  # Show iterations
    cat(paste("Iterations of EM:", "\n"))
  
  while ((!converged) & (it < maxits)) {
    probsOld = probs
    muOld = mu
    varOld = var
    riOld = ri
    
    # E
    # Compute responsibilities
    for (k in 1:clusters){
      ri[, k] = probs[k] * dnorm(X, mu[k], sd = sqrt(var[k]), log = FALSE)
    }
    
    ri = ri/rowSums(ri)
    
    # M
    rk = colSums(ri)           # rk is the weighted average cluster membership size
    probs = rk/N
    mu = (t(X) %*% ri) / rk    
    var = (t(X^2) %*% ri) / rk - mu^2
    
    # could do mu and var via log likelihood here, but this is more straightforward
    
    parmlistold     = rbind(probsOld, muOld, varOld)
    parmlistcurrent = rbind(probs, mu, var)
    
    it = it + 1
    
    # if showits true, & it =1 or divisible by 5 print message
    if (showits & it == 1 | it%%5 == 0)        
      cat(paste(format(it), "...", "\n", sep = ""))
    
    #converged = max(abs(parmlistold - parmlistcurrent)) <= tol
  }
  
  clust = which(round(ri) == 1, arr.ind = TRUE) # create cluster membership
  clust = clust[order(clust[, 1]), 2]           # order according to row rather than cluster
  
  out = list(
    probs   = probs,
    mu      = mu,
    var     = var,
    resp    = ri,
    cluster = clust
  )
  
  out
}



em_mixture <- function(params, X, clusters = 2, tol = .00001, maxits   = 100, showits  = TRUE){
  
  # Arguments are 
  # params: starting parameters (means, covariances, cluster probability)
  # X: data 
  # clusters: number of clusters desired
  # tol: tolerance
  # maxits: maximum iterations
  # showits: whether to show iterations
  
  require(mvtnorm)
  # Starting points
  N     = nrow(X)
  mu    = params$mu
  var   = params$var
  probs = params$probs
  
  # initializations
  
  # cluster 'responsibilities', i.e. probability of cluster membership for each
  # observation i
  ri = matrix(0, ncol=clusters, nrow=N)       
  ll = 0                                        # log likelihood
  it = 0                                        # iteration count
  converged = FALSE                             # convergence
  
  # Show iterations if showits == true
  if (showits)                                  
    cat(paste("Iterations of EM:", "\n"))
  
  while (!converged & it < maxits) { 
    probsOld = probs
    # muOld = mu                # Use direct values or loglike for convergence check
    # varOld = var
    llOld = ll
    riOld = ri
    
    ### E
    # Compute responsibilities
    for (k in 1:clusters){
      ri[,k] = probs[k] * dmvnorm(X, mu[k, ], sigma = var[[k]], log = FALSE)
    }
    
    ri = ri/rowSums(ri)
    
    ### M
    rk = colSums(ri)            # rk is weighted average cluster membership size
    probs = rk/N
    
    for (k in 1:clusters){
      varmat = matrix(0, ncol = ncol(X), nrow = ncol(X))    # initialize to sum matrices
      
      for (i in 1:N){
        varmat = varmat + ri[i,k] * X[i,]%*%t(X[i,])
      }
      
      mu[k,]   = (t(X) %*% ri[,k]) / rk[k]
      var[[k]] =  varmat/rk[k] - mu[k,]%*%t(mu[k,])
      
      ll[k] = -.5*sum( ri[,k] * dmvnorm(X, mu[k,], sigma = var[[k]], log = TRUE) )
    }
    
    ll = sum(ll)
    
    # compare old to current for convergence
    parmlistold =  c(llOld, probsOld)           # c(muOld, unlist(varOld), probsOld)
    parmlistcurrent = c(ll, probs)              # c(mu, unlist(var), probs)
    it = it + 1
    
    # if showits true, & it =1 or modulo of 5 print message
    if (showits & it == 1 | it%%5 == 0)         
      cat(paste(format(it), "...", "\n", sep = ""))
    
    converged = min(abs(parmlistold - parmlistcurrent)) <= tol
  }
  
  clust = which(round(ri) == 1, arr.ind = TRUE)        # create cluster membership
  clust = clust[order(clust[,1]), 2]            # order accoring to row rather than cluster
  
  
  list(
    probs   = probs,
    mu      = mu,
    var     = var,
    resp    = ri,
    cluster = clust,
    ll      = ll
  )
}


plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}



# quakes data set is in base R
data(quakes)

head(quakes)
depth = as.matrix(quakes[, 3, drop = FALSE])
options(scipen = 999)

p1 <- ggplot(quakes, aes(x = depth)) +
  geom_density()
p1
p1 + 
  geom_vline(xintercept = 80, col = "red", size = 2) + 
  geom_vline(xintercept = 565, col = "blue", size = 2)


#########starting values 
d_depth = sd(quakes$depth)
m1_depth = mean(quakes$depth) - (d_depth/2)
m2_depth=  mean(quakes$depth) + (d_depth/2)

start_values_depth = list(mu = c(m1_depth, m2_depth),
                          var = c(d_depth, d_depth),
                          probs = c(.5, .5))


#start_values_1 = list(mu = c(50, 90),
#                  var = c(1, 115),
#                  probs = c(.5, .5))
#mix_erupt   = em_mixture(start_values_1, X = eruptions,  tol = 1e-8)  

mix_waiting_depth = em_mixture(start_values_depth, X = depth, tol = 1e-3)




#data.frame(x = waiting) %>%
ggplot(quakes, aes(x = depth)) +
  geom_histogram(aes(x= depth, ..density..), binwidth = 5, colour = "black", 
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mix_waiting_depth$mu[1], sqrt(mix_waiting_depth$var[1]), lam = mix_waiting_depth$probs[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mix_waiting_depth$mu[2], sqrt(mix_waiting_depth$var[2]), lam = mix_waiting_depth$probs[2]),
                colour = "blue", lwd = 1.5) +
  ylab("Density")







# faithful data set is in base R
data(faithful)

head(faithful)
waiting = as.matrix(faithful[, 2, drop = FALSE])
options(scipen = 999)

p2 <- ggplot(faithful, aes(x = waiting)) +
  geom_density()
p2
p2 + 
  geom_vline(xintercept = 53, col = "red", size = 2) + 
  geom_vline(xintercept = 80, col = "blue", size = 2)





#########starting values 
d_faith = sd(faithful$waiting)
m1_faith = mean(faithful$waiting) - (d_faith/2)
m2_faith=  mean(faithful$waiting) + (d_faith/2)

start_values_faith = list(mu = c(m1_faith, m2_faith),
                      var = c(d_faith, d_faith),
                      probs = c(.5, .5))


mustart  = rbind(c(3, 60), c(3, 60.1))    # must be at least slightly different
covstart = list(cov(faithful), cov(faithful))
probs    = c(.01, .99)

# params is a list of mu, var, and probs 
starts = list(mu = mustart, var = covstart, probs = probs)  


#start_values_1 = list(mu = c(50, 90),
#                  var = c(1, 115),
#                  probs = c(.5, .5))
#mix_erupt   = em_mixture(start_values_1, X = eruptions,  tol = 1e-8)  

mix_waiting_faith = em_mixture(start_values_faith, X = waiting, tol = 1e-3)

mix_faithful = em_mixture(
  params = starts,
  X = as.matrix(faithful),
  clusters = 2,
  tol      = 1e-12,
  maxits   = 1500,
  showits  = TRUE
)

#data.frame(x = waiting) %>%
ggplot(faithful, aes(x = waiting)) +
  geom_histogram(aes(x= waiting, ..density..), binwidth = 1, colour = "black", 
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mix_waiting_faith$mu[1], sqrt(mix_waiting_faith$var[1]), lam = mix_waiting_faith$probs[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mix_waiting_faith$mu[2], sqrt(mix_waiting_faith$var[2]), lam = mix_waiting_faith$probs[2]),
                colour = "blue", lwd = 1.5) +
  ylab("Density")

