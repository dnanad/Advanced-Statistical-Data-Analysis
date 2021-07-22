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




##Plot the empirical probability density functions of all three samples on one panel. 
# build a dataframe with all samples 
df <- data.frame(bins = c(bin1, bin2, bin3), method=as.factor(c(rep(1, N), rep(2,N), rep(3,N))))

# plot all histograms in one plot
ggplot(data=df, alpha=0.4, aes(bins,stat(density), group=method, color=method)) +
  geom_histogram(aes(fill=method, alpha=0.6),binwidth = 1,center=1, position='dodge2')+
  geom_density(bw=0.6)+
  xlab('')+ylab('empirical probability')+
  scale_x_continuous(breaks=0:10)+
  guides(alpha='none')+
  theme_minimal()
ggsave('emp_all.png', width=8, height = 5.2)

# plot three histogramm in a column
ggplot(data=df, alpha=0.4, aes(bins,stat(density), group=method, color=method)) +
  geom_histogram(aes(fill=method, alpha=0.5), binwidth = 1,center=1)+
  geom_density(bw=0.6)+
  scale_x_continuous(breaks=0:10)+
  theme_minimal()+
  facet_wrap(facets = ~method, nrow = 3)+
  theme(strip.background = element_blank(),strip.text.x = element_blank())+
  guides(alpha='none')+
  xlab('')+ylab('empirical probability')
ggsave('emp_all_col.png', width=4.5, height = 6)
#Switch the random number generator back to its default.
RNGkind('default')
#Comment on the results. 


  ###(b) Aaccept-reject method###
set.seed(123)
N = 10000
c <- sqrt(2*pi)*exp(-1/2)
sample <- numeric()
Acc <- c(FALSE, FALSE)
k<- 0


#standard Cauchy distribution
#first we generate the normal variables and we keep the produced Cauchy variables for a later comparison of theoretical and empirical acceptance probabilities
while (sum(Acc)< N){
  k <- k+1
  u <- runif(2)
  cau <- tan(pi*u[1]-pi/2)
  sample[k] <- cau
  if (u[2]*c*dcauchy(cau)<=dnorm(cau)) {
    Acc[k] <- TRUE
  }
  else{
    Acc[k] <- FALSE
  }
}

## plot histogram and QQ-plot for simulated normals
NORM <- data.frame(x =sample[Acc]) #this is our resulting sample

ggplot() +
  geom_histogram(data=NORM,aes(x,stat(density),alpha=0.6), binwidth = 0.3, center=0)+       #,binwidth = 1,center=1)+
  xlab('')+ylab('density')+
  stat_function(geom='line', aes(x = -4:4), n = 200, fun =
                  dnorm, args = list(mean=0,sd=1), colour='orange')+
  guides(alpha='none')+
  theme_minimal()
ggsave('normal_hist.png', width = 5, height=2.3)

ggplot(data=NORM)+geom_qq(aes(sample=x), size=0.3)+theme_minimal()+geom_abline(slope=1, intercept=0)
ggsave('normal_qq.png', width=3.5, height =3.2)


## examine the acceptance probability (AP) ##
sum(Acc)/length(Acc) # empirical overall AP
1/c # theoretical overall AP


## Now look at more local APs
# we just use the range of accepted values since the Cauchy-distr. has very high tail probs, 
# and the min and max of a sample would probably avoid any reasonable illustration
# We bin and count the data using function hist to compute empirical APs
x=seq(min(sample[Acc]),max(sample[Acc]), length.out = 100)
H1 <- hist(sample[Acc], breaks = x)
Hges <- hist(sample[x[1]<sample & sample<x[100]], breaks=x)
H1$counts
Hges$counts
emp_acc_prob <- H1$counts/Hges$counts
bin_centers<- (x[1:99]+x[2:100])/2
# build a dataframe for plotting
local_acc_prob <- data.frame(bin_centers=c(bin_centers, bin_centers),
                             acc_prob=c(emp_acc_prob,dnorm(bin_centers)/(dcauchy(bin_centers)*c)),
                             color=c(rep('empirical', 99), rep('theoretical',99)))
ggplot(data=local_acc_prob)+
  geom_line(aes(bin_centers, acc_prob, colour=color))+
  scale_color_manual(name='',labels=c('empricial', 'theoretical'), values = c('red', 'black'))+
  xlab('')+ylab('Acceptance probability')+
  theme_minimal()+
  scale_x_continuous(limits=c(-4,4))
ggsave('acc_probs.png', width=4, height = 2.6)