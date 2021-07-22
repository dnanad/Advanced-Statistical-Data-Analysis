####      Ex 8     #####

### Preamble ####

library(xtable)
library("tidyverse")
library(gridExtra)


###  a) density estimation ####

###     load data set
dat <- read_csv('StudentsPerformance.csv')
dat <- dat %>% select(`test preparation course`,`math score`, `reading score`, `writing score`)
xtable(head(dat), auto=TRUE, digits = 2)


### first define different kernels for plotting and lateron in the functions
I <- function(x) {as.numeric((x<1) & (x > -1))}
rect <- Vectorize(function(x){0.5*I(x)})
tri <- Vectorize(function(x) {(1-abs(x))*I(x)})
epa <- Vectorize(function(x) {0.75*(1-x^2)*I(x)})
gau <- Vectorize(function(x) {1/(sqrt(2*pi)*exp(x^2/2))})
ggplot()+
  stat_function(aes(x=-3:3, color='rectangular'), n=200, fun=rect, geom='step')+
  stat_function(aes(x=-3:3, color='triangular'), n=200, fun=tri, geom='line')+
  stat_function(aes(x=-3:3, color='Epanechnikov'), n=200, fun=epa, geom='line')+
  stat_function(aes(x=-3:3, color='Gaussian'), n=200, fun=gau, geom='line')+
  xlab('distance')+theme_minimal()+ylab('weight')+
  scale_color_discrete(name='kernel')+theme(line = element_line(size=1))
ggsave('kernels.png', width=5, height=2.3)

## create function for density estimation ##
densitiy_est <- function(X, h, method){
  if (h<=0) {warning('h has to be positive!')}
  else { 
    n <- length(X)
    if (method =='rect') {k <- rect
    } else if (method == 'tri') {k <- tri
    } else if (method == 'epa') {k <- epa
    } else if (method == 'gau') {k <- gau}
    f_hat <- Vectorize(function(x) {1/(n*h)*sum(k((x-X)/h))})
    kern_smooth <- list('f_hat'=f_hat, 'method'=method, 'bandwidth'=h, 'kernel'=k)
    return(kern_smooth)
  }
}


## now use different bandwidths (2,6, 15 and 30) and visualise
ggplot()+
  geom_histogram(data=dat, aes(`math score`, stat(density)), alpha=0.5, linetype=1, fill='gray', colour='black', bins=30)+
  stat_function(aes(x=0:100, color='2'), n=2000,fun=densitiy_est(dat$`math score`, 2, 'epa')$f_hat, geom='line')+
  stat_function(aes(x=0:100, color='6'), n=2000,fun=densitiy_est(dat$`math score`, 6, 'epa')$f_hat, geom='line')+
  stat_function(aes(x=0:100, color='15'), n=2000,fun=densitiy_est(dat$`math score`, 15, 'epa')$f_hat, geom='line')+  
  stat_function(aes(x=0:100, color='30'), n=2000,fun=densitiy_est(dat$`math score`, 30, 'epa')$f_hat, geom='line')+
  theme_minimal()+xlab('math score')+
  scale_color_discrete(name='bandwidth', breaks=c('2','6','15','30'))
ggsave('bandwidths.png', width=5.5, height=2.3)


### this is the first try. Another posibility to visualise without ggplot - nevermind  
# hist(dat$`math score`, freq=FALSE, xlim=c(0, 105))
# points(table(dat$`math score`)/(sum(dat$`math score`)/2))
# Est <- tibble()
# for (i in c(1,2, 6, 15, 30)) {
#   h <- i
#   sm <- densitiy_est(dat$`math score`, h, 'epa')
#   X <- sort(dat$`math score`)
#   Est <- sm$f_hat(X)
#   points(X, Est, type='l', col=i)
# }



## use and visualise different kernels
# it seemed for the epanechnikov kernel that h in (6,15) may be a good choice --> set h= 10
ggplot()+
  geom_histogram(data=dat, aes(`math score`, stat(density)), alpha=0.5, linetype=1, fill='gray', colour='black', bins=30)+
  stat_function(aes(x=0:100, color='rectangular'), n=2000,fun=densitiy_est(dat$`math score`, 10, 'rect')$f_hat, geom='line')+
  stat_function(aes(x=0:100, color='triangular'), n=2000,fun=densitiy_est(dat$`math score`, 10, 'tri')$f_hat, geom='line')+
  stat_function(aes(x=0:100, color='Epanechnikov'), n=2000,fun=densitiy_est(dat$`math score`, 10, 'epa')$f_hat, geom='line')+  
  stat_function(aes(x=0:100, color='Gaussian'), n=2000,fun=densitiy_est(dat$`math score`, 10, 'gau')$f_hat, geom='line')+
  theme_minimal()+xlab('math score')+
  scale_color_discrete(name='kernel', breaks=c('rectangular','triangular', 'Epanechnikov','Gaussian'))
ggsave('kernelfits.png', width=6, height=2.3)



###   b) cross-validation method ####

## define CV-function
CV <- function(data, h, method){
  n <- length(data)
  sm <- densitiy_est(data, h, method)
  k <- sm$kernel
  integrand <- Vectorize(function(x)(sm$f_hat(x)^2))
  a <- integrate(integrand, 0, 100)$value
  B <- 0
  for (i in 1:n) {
    B <- B+sum(k((data[-i]-data[i])/h))
  }
  CrossV <- a-B*2/((n-1)*n*h)
  return(CrossV)
}
CV <- Vectorize(CV, vectorize.args = 'h')


## a first visualisation of the CV:
#CV(dat$`math score`, h, 'epa')
H <- seq(1, 25, length.out = 40) 
CV_H <- CV(dat$`math score`,H, 'epa')
CVdf <- tibble(h=H, CV=CV_H)
which(CV_H==min(CV_H))
minCV<- CVdf[which(CV_H==min(CV_H)),]
minCV$CV
min_math_epa <- minCV

# now use optimize method to search for best bandwidth
# use different limits to show we get very different results
opt_epa_math<-as.data.frame(optimize(CV, lower=0, upper=30, data=dat$`math score`, method = 'epa'))
opt_epa_math2<-as.data.frame(optimize(CV, lower=0, upper=13, data=dat$`math score`, method = 'epa'))
opt_epa_math3<-as.data.frame(optimize(CV, lower=1, upper=5, data=dat$`math score`, method = 'epa'))

## see in the plot that the results vary for different limits (even when the first optimum is still contained in the limits of the second try)
ggplot()+
  geom_line(data=CVdf, aes(x=h, y=CV))+
  geom_point(data=minCV, aes(x=h, y=CV), color='red')+
  geom_point(data=opt_epa_math, aes(x=minimum, y=objective, color='0,30'))+
  geom_point(data=opt_epa_math2, aes(x=minimum, y=objective,color='0,13'))+
  geom_point(data=opt_epa_math3, aes(x=minimum, y=objective,color='1,5'))+
  scale_color_discrete(name='Limits for optimize')
#annotate('text', x=4.6+minCV$h, y=minCV$CV, label=paste('h =', round(minCV$h,2)), size=7)
ggsave('CV1.png', width=6, height=2.3)



#### now we do a unified stragety with a combination of scan and optimize
####  with a Gaussian kernel and Epanechnikov 
####  afterwards we use the R implemented functions bw.ucv and bw.bcv and bring all the results together
####  Do this for all scores and bring the data together in a final data frame:

## math score
# Guassian kernel and my own CV function:
H <- seq(1, 20, length.out = 40) 
CV_H <- CV(dat$`math score`,H, 'gau')
CVdf <- tibble(h=H, CV=CV_H)
minCV <- CVdf[which(CV_H==min(CV_H)),]
minCV_math_gau<- minCV
opt_math_gau<-as.data.frame(optimize(CV, lower=4.5, upper=5.5, data=dat$`math score`, method = 'gau'))

# do a plot
g_math_gau<- ggplot()+
  geom_line(data=CVdf, aes(x=h, y=CV))+
  geom_point(data=minCV, aes(x=h, y=CV), color='red')+
  geom_point(data=opt_math_gau, aes(x=minimum, y=objective), color='blue')+
  theme_minimal()

# now use the R methods provided
bwucv_math <- bw.ucv(dat$`math score`, lower=1, upper=50) # I restrict the search bounds to a reasonable area
bwbcv_math <- bw.bcv(dat$`math score`, lower=1, upper=50) # the results are different, since they use the Gauss-kernel!


## reading score
# Epanechnikov kernel and my own CV function:
CV_H <- CV(dat$`reading score`,H, 'epa')
CVdf <- tibble(h=H, CV=CV_H)
minCV <- CVdf[which(CV_H==min(CV_H)),]
minCV_read_epa<- minCV
opt_read_epa<-as.data.frame(optimize(CV, lower=1, upper=2, data=dat$`reading score`, method = 'epa'))
# do a plot
g_read_epa<- ggplot()+
  geom_line(data=CVdf, aes(x=h, y=CV))+
  geom_point(data=minCV, aes(x=h, y=CV), color='red')+
  geom_point(data=opt_read_epa, aes(x=minimum, y=objective), color='blue')+
  theme_minimal()
# Guassian kernel and my own CV function:
CV_H <- CV(dat$`reading score`,H, 'gau')
CVdf <- tibble(h=H, CV=CV_H)
minCV <- CVdf[which(CV_H==min(CV_H)),]
minCV_read_gau<- minCV
opt_read_gau<-as.data.frame(optimize(CV, lower=4, upper=5, data=dat$`reading score`, method = 'gau'))
# do a plot
g_read_gau<- ggplot()+
  geom_line(data=CVdf, aes(x=h, y=CV))+
  geom_point(data=minCV, aes(x=h, y=CV), color='red')+
  geom_point(data=opt_read_gau, aes(x=minimum, y=objective), color='blue')+
  theme_minimal()
g_read_gau
# now use the R methods provided
bwucv_read <- bw.ucv(dat$`reading score`, lower=1, upper=50) 
bwbcv_read <- bw.bcv(dat$`reading score`, lower=1, upper=50) 


## writing score
# Epanechnikov kernel and my own CV function:
CV_H <- CV(dat$`writing score`,H, 'epa')
CVdf <- tibble(h=H, CV=CV_H)
minCV <- CVdf[which(CV_H==min(CV_H)),]
minCV_write_epa<- minCV
opt_write_epa<-as.data.frame(optimize(CV, lower=1, upper=2, data=dat$`writing score`, method = 'epa'))
# do a plot
g_write_epa<- ggplot()+
  geom_line(data=CVdf, aes(x=h, y=CV))+
  geom_point(data=minCV, aes(x=h, y=CV), color='red')+
  geom_point(data=opt_write_epa, aes(x=minimum, y=objective), color='blue')+
  theme_minimal()
# Guassian kernel and my own CV function:
CV_H <- CV(dat$`writing score`,H, 'gau')
CVdf <- tibble(h=H, CV=CV_H)
minCV <- CVdf[which(CV_H==min(CV_H)),]
minCV_write_gau<- minCV
opt_write_gau<-as.data.frame(optimize(CV, lower=4, upper=5, data=dat$`writing score`, method = 'gau'))
# do a plot
g_write_gau<- ggplot()+
  geom_line(data=CVdf, aes(x=h, y=CV))+
  geom_point(data=minCV, aes(x=h, y=CV), color='red')+
  geom_point(data=opt_write_gau, aes(x=minimum, y=objective), color='blue')+
  theme_minimal()
g_write_gau
# now use the R methods provided
bwucv_write <- bw.ucv(dat$`writing score`, lower=1, upper=50) 
bwbcv_write <- bw.bcv(dat$`writing score`, lower=1, upper=50) 


### now bring together the results
## first the data frame
res_math <- as.numeric(c(min_math_epa, opt_math_gau, 
                         bwbcv_math, CV(dat$`math score`,bwbcv_math, 'gau'), 
                         bwucv_math, CV(dat$`math score`,bwucv_math, 'gau')))
res_read <- as.numeric(c(opt_read_epa, opt_read_gau, 
                         bwbcv_read, CV(dat$`reading score`,bwbcv_read, 'gau'), 
                         bwucv_read, CV(dat$`reading score`,bwucv_read, 'gau')))
res_write <- as.numeric(c(opt_write_epa, opt_write_gau, 
                          bwbcv_write, CV(dat$`writing score`,bwucv_write, 'gau'),
                          bwucv_write, CV(dat$`writing score`,bwucv_write, 'gau')))
res_df <- as.data.frame(rbind(res_math, res_read, res_write))
row.names(res_df)<- c('math', 'reading', 'writing')
colnames(res_df) <-c('hepa','CVepa','hgau','CVgau','bwbcvh','bwbcvCV',
                     'bwucvh','bwucvCV')
res_df <- as_tibble(res_df)
res_df
xtable(res_df, auto=TRUE, digits=c(1,2,4,2,4,2,4,2,4))
## now bring all the plots together
p <- arrangeGrob(g_math_gau,
                 g_read_epa,g_read_gau,
                 g_write_epa,g_write_gau, layout_matrix = matrix(c(1,1,2,3,4,5), nrow = 3, byrow=TRUE))
ggsave('CVplots.png', width=10, height = 12,p)



###   c) comparison of the two groups ####

prep_stud <- dat %>% filter(`test preparation course`=='completed')
nonprep_stud <- dat %>% filter(`test preparation course`=='none')
dens_mat_prep <- densitiy_est(prep_stud$`math score`, min_math_epa$h, 'epa')
dens_mat_nonprep <- densitiy_est(nonprep_stud$`math score`, min_math_epa$h, 'epa')

# following makes it short and beautiful, but pls dont ask how long it took me to come to this
test <- dat%>% gather(`math score`,`reading score`,`writing score`, key='domain', value='score')
test
h_CV <- mean(c(min_math_epa$h, opt_read_epa$minimum, opt_write_epa$minimum))
x11()
ggplot(data=test, aes(score,stat(density), group=`test preparation course`))+
  facet_wrap(~domain, nrow=1)+
  geom_histogram(aes(fill=`test preparation course`, alpha=0.2),position='identity', binwidth = 5)+
  geom_density(aes(col=`test preparation course`), bw=h_CV)+ 
  theme(legend.position="top")+ 
  guides(alpha='none')
ggsave('groups.png', width=7, height=3)

## this was my first (well actually my 201st) try
## therefore we would need the gridExtra pckg
# p1 <- ggplot(data=dat, aes(`math score`,stat(density), group=`test preparation course`))+
#   geom_histogram(aes(fill=`test preparation course`, alpha=0.2),position='dodge2', binwidth = 5)+
#   geom_density(aes(col=`test preparation course`), bw=h_CV_math) +
#   theme(legend.position="none")+ xlim(c(0,100))+ylim(c(0,0.036))
# 
# p2 <- ggplot(data=dat, aes(`reading score`,stat(density), group=`test preparation course`))+
#   geom_histogram(aes(fill=`test preparation course`, alpha=0.2),position='dodge2', binwidth = 5)+
#   geom_density(aes(col=`test preparation course`), bw=h_CV_reading) +
#   guides(alpha='none')+
#   theme(legend.position = c(0.4, 0.98), legend.direction = 'horizontal', 
#         legend.text.align = 0.5, axis.title.y = element_blank())+ 
#   xlim(c(0,100))+ylim(c(0,0.036))
# 
# p3 <- ggplot(data=dat, aes(`writing score`,stat(density), group=`test preparation course`))+
#   geom_histogram(aes(fill=`test preparation course`, alpha=0.2),position='dodge2', binwidth = 5)+
#   geom_density(aes(col=`test preparation course`), bw=h_CV_writing) +
#   guides(alpha='none', size='none')+#stat_bin(geom='bar', position='dodge2')
#   theme(legend.position="none", axis.title.y = element_blank()) + 
#   xlim(c(0,100))+ylim(c(0,0.036))
# x11()
# grid.arrange(p1, p2,p3, widths=c(1,1,1), nrow = 1)