###############################
### EM for gaussian mixture ###
###############################
gaussmixEM = function(params, X, clusters = 2, tol=.00001, maxits=100, showits=T){     
  # Arguments are starting parameters (means, covariances, cluster probability), data, number of clusters desired, tolerance,
  # maximum iterations, and whether to show iterations
  
  # Starting points
  N = nrow(X)
  nams = names(params)
  mu = params$mu
  var = params$var
  probs = params$probs
  
  # Other initializations
  ri = matrix(0, ncol=clusters, nrow=N)  #initialize cluster 'responsibilities', i.e. probability of cluster membership for each observation i
  it = 0
  converged = FALSE
  
  if (showits)                                  # Show iterations
    cat(paste("Iterations of EM:", "\n"))
  
  while ((!converged) & (it < maxits)) { 
    probsOld = probs
    muOld = mu
    varOld = var
    riOld = ri
    
    ### E
    # Compute responsibilities
    for (k in 1:clusters){
      ri[,k] = probs[k] * dnorm(X, mu[k], sd = sqrt(var[k]), log=F)
    }
    ri = ri/rowSums(ri)
    
    ### M
    rk = colSums(ri)                           # rk is the weighted average cluster membership size
    probs = rk/N
    mu = (t(X) %*% ri) / rk                      # could do mu and var via log likelihood here but this is more straightforward
    var = (t(X^2) %*% ri) / rk - mu^2
    
    parmlistold = rbind(probsOld, muOld, varOld)
    parmlistcurrent = rbind(probs, mu, var)
    it = it + 1
    if (showits & it == 1 | it%%5 == 0)        # if showits true, & it =1 or divisible by 5 print message
      cat(paste(format(it), "...", "\n", sep = ""))
    converged = max(abs(parmlistold - parmlistcurrent)) <= tol
  }
  
  clust = which(round(ri)==1, arr.ind=T)       # create cluster membership
  clust = clust[order(clust[,1]), 2]           # order accoring to row rather than cluster
  out = list(probs=probs, mu=mu, var=var, resp=ri, cluster=clust)
} 


### This example uses Old Faithful geyser eruptions.  This is only a univariate mixture for either eruption time or wait time.
### The next will be doing both variables, i.e. multivariate normal.  'Geyser' is supposedly more accurate, though seems to have 
### arbitrarily assigned some duration values.  See also http://www.geyserstudy.org/geyser.aspx?pGeyserNo=OLDFAITHFUL, but that only has
### intervals. Some July 1995 data is available

### faithful data set
data(faithful)
head(faithful)

# starting parameters; requires mean, variance and class probabilitiy
params1 = list(mu=c(2, 5), var=c(1, 1), probs=c(.5, .5))  # note that starts from mean must be in data range or it will break.  
params2 = list(mu=c(0, 90), var=c(1, 15), probs=c(.5, .5))  

X1 = matrix(faithful[,1])
X2 = matrix(faithful[,2])

test = gaussmixEM(params2, X=X2, tol=1e-3)

plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}
test$probs[1]

df <- cbind(data.frame(test), x = X2)
df %>%
  ggplot() +
    geom_histogram(aes(x= x, ..density..), binwidth = 1, colour = "black", 
                   fill = "white") +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(mu = df$mu.1[1],sigma = df$var.1[1], lam = test$probs[1]),
                  colour = "red", lwd = 1.5) +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(mu=df$mu.2[2], sigma = df$var.2[2], lam = test$probs[2]),
                  colour = "blue", lwd = 1.5) +
    ylab("Density")



ggplot(faithful, aes(x = X2)) +
  geom_histogram(aes(x= X2), binwidth = 1, colour = "black", 
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(test$mu[1], test$var[1], lam = test$probs[1]),
                colour = "red", lwd = 1.5)+
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(test$mu[2], test$var[2], lam = test$probs[2]),
                colour = "blue", lwd = 1.5) +
  ylab("Density") 



plot_mix_comps(X2, test$mu[1], test$var[1], lam = test$probs[1])

test$resp[1]*dnorm(X2, test$mu[1], test$var[1])






# NOT RUN {
##Analyzing the Old Faithful geyser data with a 2-component mixture of normals.
library("mixtools")
data(faithful)
attach(faithful)
set.seed(100)
system.time(mixmdl<-normalmixEM(waiting, arbvar = FALSE, epsilon = 1e-03))

data.frame(x = mixmdl$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                colour = "blue", lwd = 1.5) +
  ylab("Density")
