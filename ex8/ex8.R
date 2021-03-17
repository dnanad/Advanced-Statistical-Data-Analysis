####      Ex 8     #####

### Preamble ####


## install/load required packages

#install.packages("xtable")
library(xtable)
#install.packages('tidyverse')
library("tidyverse")
#install.packages('gridExtra')
library(gridExtra)


## Exercises ##
### (a) ####
#Implement  kernel  density  estimation  in  an  R  function  that  depends  on  the  sample,bandwidth  and  a  kernel.
## Read  the  data  into  R

sp_data <- read_csv('StudentsPerformance.csv')
sp_data <- sp_data %>% select(`test preparation course`,`math score`, `reading score`, `writing score`)
xtable(head(sp_data,10), auto=TRUE, digits = 2)



### Different types of kernels
#indicator function
I <- function(x) {as.numeric((x<1) & (x > -1))}
#the rectangular or uniform kernel
rect_ker <- Vectorize(function(x){0.5*I(x)})
#the triangular kernel
tri_ker <- Vectorize(function(x) {(1-abs(x))*I(x)})
#the Epanechnikov kernel
epa_ker <- Vectorize(function(x) {0.75*(1-x^2)*I(x)})
#the Gaussian kernel
gauss_ker <- Vectorize(function(x) {1/(sqrt(2*pi)*exp(x^2/2))})

ggplot(data.frame(x = c(-3,3)), aes(x))+
  stat_function(aes(color='rectangular'), n=500, fun=rect_ker, geom='step')+
  stat_function(aes(color='triangular'), n=500, fun=tri_ker, geom='line')+
  stat_function(aes(color='Epanechnikov'), n=500, fun=epa_ker, geom='line')+
  stat_function(aes(color='Gaussian'), n=500, fun=gauss_ker, geom='line')+
  xlab('distance')+theme_minimal()+ylab('weight')+
  scale_color_discrete(name='Kernel')+theme(line = element_line(size=1))
ggsave('types_kernels.png', width=5, height=2.3)


## function for density estimation ##
densitiy_estimator <- function(X, h, method){
  if (h<=0) {warning('h is less equal zero, which is not allowed. Bandwidth has to be positive!')}
  else { 
    n <- length(X)
    if (method =='rect_ker') {
      k <- rect_ker
    } else if (method == 'tri_ker') {
      k <- tri_ker
    } else if (method == 'epa_ker') {
      k <- epa_ker
    } else if (method == 'gauss_ker') {
      k <- gauss_ker
    }
    f_hat <- Vectorize(function(x) {1/(n*h)*sum(k((x-X)/h))})
    kern_smooth <- list('f_hat'= f_hat, 'method'=method, 'bandwidth'= h, 'kernel'= k)
    return(kern_smooth)
  }
}

#Estimate  and  plot  the  density  of math.score with the Epanechnikov kernel and 4 different choices of the bandwidth putting them onto one plot.
## now we use different bandwidths (2,5, 12 and 30) and visualize

h1 <- densitiy_estimator(sp_data$`math score`, 2, 'epa_ker')#$f_hat
h2 <- densitiy_estimator(sp_data$`math score`, 5, 'epa_ker')
h3 <- densitiy_estimator(sp_data$`math score`, 12, 'epa_ker')
h4 <- densitiy_estimator(sp_data$`math score`, 30, 'epa_ker')

ggplot(data.frame(x = c(0,100)))+
  geom_histogram(data=sp_data, aes(`math score`, stat(density)), alpha=0.1, linetype=1, fill='yellow', colour='black', bins=30)+
  stat_function(aes(color='2'), n=2000,fun=h1$f_hat, geom='line')+
  stat_function(aes(color='5'), n=2000,fun=h2$f_hat, geom='line')+
  stat_function(aes(color='12'), n=2000,fun=h3$f_hat, geom='line')+  
  stat_function(aes(color='30'), n=2000,fun=h4$f_hat, geom='line')+
  theme_minimal()+xlab('math score')+
  scale_color_discrete(name='bandwidth', breaks=c('2','5','12','30'),theme(line = element_line(size=1)))
ggsave('diff_bandwidths.png', width=5.5, height=2.3)


#Fix the bandwidth you find most reasonable and plot kernel density estimators for 4 different choices of kernel functions (Gauss, Epanech-nikov, uniform, tridiagonal),
#putting them onto one plot.

#h between 5 and 12 seems a good choice h<-8

bw <- 8 

kde1 <- densitiy_estimator(sp_data$`math score`, bw, 'rect_ker')$f_hat
kde2 <- densitiy_estimator(sp_data$`math score`, bw, 'tri_ker')$f_hat
kde3 <- densitiy_estimator(sp_data$`math score`, bw, 'epa_ker')$f_hat
kde4 <- densitiy_estimator(sp_data$`math score`, bw, 'gauss_ker')$f_hat



ggplot(data.frame(x = c(0,100)))+
  geom_histogram(data=sp_data, aes(`math score`, stat(density)), alpha=0.1, linetype=1, fill='yellow', colour='black', bins=30)+
  stat_function(aes(color='rectangular'), n=2000,fun=kde1, geom='line')+
  stat_function(aes(color='triangular'), n=2000,fun=kde2, geom='line')+
  stat_function(aes(color='Epanechnikov'), n=2000,fun=kde3, geom='line')+  
  stat_function(aes(color='Gaussian'), n=2000,fun=kde4, geom='line')+
  theme_minimal()+xlab('math score')+
  scale_color_discrete(name='kernel', breaks=c('rectangular','triangular', 'Epanechnikov','Gaussian'))
ggsave('diff_ker.png', width=5.5, height=2.3)


### (b) ###
#Implement the cross-validation criterion to find the optimal bandwidth.


densitiy_estimator(sp_data$`math score`, bw, 'gauss_ker')$kernel

## Cross Validation function
CV <- function(sample, h, method){
  n <- length(sample)
  de <- densitiy_estimator(sample, h, method)
  K <- de$kernel
  integrand <- Vectorize(function(x)(de$f_hat(x)^2))
  a <- integrate(integrand, lower = 0, upper = 100)$value
  b <- 0
  for (i in 1:n) {
    b <- b+sum(k((sample[-i]-sample[i])/h))
  }
  CrossVal <- a-b*2/((n-1)*n*h)
  return(CrossVal)
}



#Compare your resulting bandwidths with the ones obtained by built-in R functions bw.ucv and bw.bcv used indensity for all three scores math.score, reading.score and writing.score.
