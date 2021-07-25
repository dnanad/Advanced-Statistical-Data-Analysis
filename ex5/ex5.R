####      Exercise 5      #####

### Preamble ####



library(tidyverse)
library("ggplot2")
library("dplyr")


#####EM-Algorithm to estimate the parameters of a Gaussian mixture distribution
###includes the explicit formulas for updating the means and variances of the mixture components as well as the weights
######a stopping criterion for the EM-Algorithm

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
  ll = 0                                        # log likelihood
  ll_list <-list()
  it = 0
  it_list <- list()
  converged = FALSE
  
  if (showits)                                  # Show iterations
    cat(paste("Iterations of EM:", "\n"))
  
  while ((!converged) & (it < maxits)) {
    probsOld = probs
    #muOld = mu
    #varOld = var
    llOld = ll
    riOld = ri
    # E
    for (k in 1:clusters){
      ri[, k] = probs[k] * dnorm(X, mu[k], sd = sqrt(var[k]), log = FALSE)
    }
    
    ri = ri/rowSums(ri)
    
    
    ### M
    rk = colSums(ri)            # rk is weighted average cluster membership size
    probs = rk/N
      
    #log likelihood 
    
    for (k in 1:clusters){
      varmat = matrix(0, ncol = ncol(X), nrow = ncol(X))    # initialize to sum matrices
      
      for (i in 1:N){
        varmat = varmat + ri[i,k] * X[i,]%*%t(X[i,])
      }
      
      mu[k]   = (t(X) %*% ri[,k]) / rk[k]
      var[[k]] =  varmat/rk[k] - mu[k]%*%t(mu[k])
      
      ll[k] = sum( ri[,k] * dnorm(X, mu[k], sqrt(var[k]), log = FALSE))
    }
    
    ll = sum(ll) 
    
    # compare old to current for convergence
    parmlistold =  c(llOld, probsOld)           
    parmlistcurrent = c(ll, probs)              
    it = it + 1
    ll_list[[it]]<-ll
    it_list[[it]]<-it
    # if showits true, & it =1 or modulo of 5 print message
    if (showits)# & it == 1 )#| it%%5 == 0)         
      cat(paste(format(it), "...", "\n", sep = ""))
    
    converged = min(abs(parmlistold - parmlistcurrent)) <= tol ##### a stopping criterion for the EM-Algorithm
  }
  
  clust = which(round(ri) == 1, arr.ind = TRUE)        # create cluster membership
  clust = clust[order(clust[,1]), 2]            # order accoring to row rather than cluster
  
  
  out = list(
    probs   = probs,
    mu      = mu,
    var     = var,
    resp    = ri,
    cluster = clust,
    ll = ll,
    ll_list = ll_list,
    it_list = it_list
  )
  
  out

}





plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}



###### quakes data set is in base R######
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
ggsave('quakes_depth_dist.png')

#########starting values 
d_depth = sd(quakes$depth)
m1_depth = mean(quakes$depth) - (d_depth/2)
m2_depth=  mean(quakes$depth) + (d_depth/2)

start_values_depth = list(mu = c(m1_depth, m2_depth),
                          var = c(d_depth, d_depth),
                          probs = c(.5, .5))


mix_depth = em_mixture(start_values_depth, X = depth, tol = 1e-3)

str(mix_depth)


#data.frame(x = waiting) %>%
ggplot(quakes, aes(x = depth)) +
  geom_histogram(aes(x= depth, ..density..), binwidth = 5, colour = "black", 
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mix_depth$mu[1], sqrt(mix_depth$var[1]), lam = mix_depth$probs[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mix_depth$mu[2], sqrt(mix_depth$var[2]), lam = mix_depth$probs[2]),
                colour = "blue", lwd = 1.5) +
  ylab("Density")
ggsave('quakes_depth_mix.png')

df_depth<-do.call(rbind, Map(data.frame, Iteration=mix_depth$it_list, Value=mix_depth$ll_list))

ggplot(df_depth, aes(x=Iteration, y=Value)) + 
  geom_line()+geom_point()
ggsave('quakes_depth_iter.png')




####### faithful data set is in base R########
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
ggsave('faith_wait_dist.png')




#########starting values 
d_faith = sd(faithful$waiting)
m1_faith = mean(faithful$waiting) - (d_faith/2)
m2_faith=  mean(faithful$waiting) + (d_faith/2)

start_values_faith = list(mu = c(m1_faith, m2_faith),
                      var = c(d_faith, d_faith),
                      probs = c(.5, .5))


mix_waiting_faith = em_mixture(start_values_faith, X = waiting, tol = 1e-3, maxits = 1000)

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
ggsave('faith_wait_mix.png')

df_faith<-do.call(rbind, Map(data.frame, Iteration=mix_waiting_faith$it_list, Value=mix_waiting_faith$ll_list))

ggplot(df_faith, aes(x=Iteration, y=Value)) + 
  geom_line()+geom_point()
ggsave('faith_wait_iter.png')
