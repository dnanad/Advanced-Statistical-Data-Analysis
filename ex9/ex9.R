####      Exercise 9      #####

### Preamble ####

library(xtable)
library("tidyverse")
library(gridExtra)


#available.packages() %>%
#  as_tibble() %>%
#  filter(Package == 'NonpModelCheck')

#readr_desc = available.packages() %>%
#  as_tibble() %>%
#  filter(Package == 'NonpModelCheck')

#readr_desc

#install.packages('NonpModelCheck')
library(NonpModelCheck)


###  a) function to fit spline ####

###     load data set
dat <- read_csv2('children_ex9.csv')
dat <- dat %>% select(hypage, zwast)
hypage <- dat$hypage
zwast <- dat$zwast

## create function for estimation ####
LocPolyReg <- function(Y, X, h, deg, method, der=deg, compute_weights=FALSE) {
  # first get right kernel
  I <- function(x) {as.numeric((x<=1) & (x >= -1))}
  if (method =='rect') {k <- function(x){0.5*I(x)}
  } else if (method == 'tri') {k <- Vectorize(function(x) {(1-abs(x))*I(x)})
  } else if (method == 'epa') {k <- Vectorize(function(x) {0.75*(1-x^2)*I(x)})
  } else if (method == 'gau') {k <- Vectorize(function(x) {1/(sqrt(2*pi)*exp(x^2/2))})}
  # determine values in X
  # this is because we dont want to fit a lm for every observation but just for every different value
  # X attains --> this saves a lot of work! Later on we can pick the desired values for the observations
  # (see e.g. function GCV)
  Range <- sort(unique(X))
  R <- length(Range)
  # then create matrix f_hat with the estimated coefficients/derivates
  # and weights (which we need later on for the GCV)
  f_hat <- matrix(nrow=deg+1, ncol=R)
  colnames(f_hat) <- Range
  Weights <- NULL
  for (j in 1:R){
    x <- Range[j]
    X_x <- X-x
    Xmat <- rep(1, length(X_x))
    for (l in 1:deg) {Xmat <- cbind(Xmat, X_x^l)}
    V <- k((X_x)/h)
    M <- lm(Y~Xmat-1, weights = V)
    coeff <- M$coefficients
    f_hat[,j] <-  coeff*factorial(seq(0,deg))
    if (compute_weights==TRUE) {
      obs <- which(X==x)[1]
      if (obs%in% labels(influence(M)$hat)){
        Weights[j] <- influence(M)$hat[as.character(obs)]
      } else {Weights[j] <- 0}
    }
  }
  # if more derivates than degree wished, fill up with 0s
  if (der > deg) {for (j in (deg+1):der) f_hat<-rbind(f_hat, rep(0, R))}
  if (compute_weights==TRUE) {output <- list(f_hat=f_hat, Weights=Weights)}
  else {output <- list(f_hat=f_hat)}
  return(output)
}


## plot regression for different bandwidths
bandwidths <- tibble()
h <- c(1,2,10,40)
for (j in 1:4){
  A <- LocPolyReg(zwast, hypage, h[j], 1, 'epa', 3)$f_hat
  bandwidths <- rbind(bandwidths, cbind(seq(0,59), rep(h[j],60), A[1,]))
}
names(bandwidths)<- c('age', 'bandwidth', 'reg')

bandwidths <- bandwidths %>% group_by(bandwidth)
bandwidths$bandwidth <- as.character(bandwidths$bandwidth)

g1 <- ggplot()+
  geom_point(data=dat, aes(x=hypage, y=zwast), alpha=0.2, size=1)
g2<- ggplot()+
  geom_line(data=bandwidths, mapping=aes(x=age, y=reg, group=bandwidth, color=bandwidth))+
  theme(legend.position = 'top')+
  #alpha=0.8, size=1.01
  scale_color_discrete(name='bandwidth', breaks=c('1','2','10','40'))
#ggsave('bandwidths.png', width=7, height=6)

## plot regression for different kernels
kernels <- data.frame()
h <- 10
method <- c('rect','tri','epa','gau')
for (j in 1:4){
  A <- LocPolyReg(zwast, hypage, h, 1, method[j])$f_hat
  kernels <- rbind(kernels, data.frame(seq(0,59), rep(method[j],60), A[1,]))
}
names(kernels)<- c('age', 'kernel', 'reg')
kernels <- kernels %>% as_tibble()
kernels <- kernels %>% group_by(kernel)

g3 <- ggplot()+
  geom_line(data=kernels, mapping=aes(x=age, y=reg, group=kernel, color=kernel))+
  theme(legend.position = 'top')
p <- arrangeGrob(g1, g2, g3, nrow=3)#, widths=c(0.7, 0.9,1))
ggsave('polyfits.png', width=7, height=8, p)  



###   b) find best bandwidth with GCV ####
## we write the function for GCV with arguments bandwidth (h) and degree (deg) only
## since all the other stuff is fixed

GCV <- function(h,deg){
  fit <- LocPolyReg(zwast,hypage,h, deg, 'epa', deg, compute_weights = TRUE)
  Y_hat <- fit$f_hat[1,][hypage+1]
  num<- sum((zwast-Y_hat)^2)
  W <- fit$Weights[hypage+1]
  den <- 1-(1/length(zwast))*sum(W)
  return(num/den)
}
GCV <- Vectorize(GCV)
#GCV(2,2)


## now search for optimal h for degree 1:
h_opt <- as.data.frame(optimize(GCV, c(0,60),deg=1))
#h_opt <- list(minimum=2.111064, objective=6517.513)
H <- seq(0.1,4, length.out = 30)
CV_H <- GCV(H, 1)
CVdf <- tibble(h=x, CV=CV_H)
which(CV_H==min(CV_H))
minCV_1<- CVdf[which(CVdf$CV ==min(CVdf$CV)),]

ggplot()+
  geom_line(data=CVdf, aes(x=h, y=CV))+
  geom_point(data=minCV_1, aes(x=h, y=CV,color='scanning') )+
  geom_point(data=h_opt, aes(x=minimum, y=objective, color='optimize'))+
  scale_color_discrete(name='Method')
ggsave('opth.png', width=5.5, height=2.3)

## search for deg=2
CV_H <- GCV(H, 2)
CVdf <- tibble(h=x, CV=CV_H)
which(CV_H==min(CV_H))
minCV_2<- CVdf[which(CVdf$CV ==min(CVdf$CV)),]

## search for deg=3
CV_H <- GCV(H, 3)
CVdf <- tibble(h=x, CV=CV_H)
which(CV_H==min(CV_H))
minCV_3<- CVdf[which(CVdf$CV ==min(CVdf$CV)),]

## search for deg=4
CV_H <- GCV(H, 4)
CVdf <- tibble(h=x, CV=CV_H)
which(CV_H==min(CV_H))
minCV_4<- CVdf[which(CVdf$CV ==min(CVdf$CV)),]

## bring results together in a table
bdwdth<- rbind(minCV_1,minCV_2,minCV_3, minCV_4)
xtable(bdwdth, digits = 2, auto=TRUE)

## now fit the optimal bandwidth curves for all degrees and plot results
degrees <- data.frame()
h <-bdwdth$h
Deg <- c('1','2','3', '4')
for (j in 1:4){
  A <- LocPolyReg(zwast, hypage, h[j], j, 'epa')$f_hat
  reg_fun <- reg(A)
  degrees <- rbind(degrees, data.frame(seq(0,59), rep(Deg[j],60), A[1,]))
}
names(degrees)<- c('age', 'degree', 'reg')
degrees <- degrees %>% as_tibble()
degrees <- degrees %>% group_by(degree)
ggplot()+
  geom_line(data=degrees, mapping=aes(x=age, y=reg, group=degree, color=degree))+
  scale_color_discrete()
ggsave('fitsh.png', width=5.5, height=4)


###   c) derivative regression using implemented function localpoly.reg ####
# use the optimal bandwidth vector h defined above
LPR1 <- localpoly.reg(hypage, zwast, bandwidth=h[1], degree.pol = 1, deriv = 1)
LPR2 <- localpoly.reg(hypage, zwast, bandwidth=h[2], degree.pol = 2, deriv = 1)
LPR3 <- localpoly.reg(hypage, zwast, bandwidth=h[3], degree.pol = 3, deriv = 1)
LPR4 <- localpoly.reg(hypage, zwast, bandwidth=h[4], degree.pol = 4, deriv = 1)
# build a tibble for easy plotting:
df <- tibble(age=LPR1$x, '1'=LPR1$predicted, '2'=LPR2$predicted,
             '3'=LPR3$predicted, '4'=LPR4$predicted)
df <- df %>%gather('1','2','3', '4', key='degree', value='predicted')%>% group_by(degree)
g1 <- ggplot(data=df)+
  geom_line(aes(x=age, y=predicted, color=degree))+
  theme(legend.position='top')+ylim(-0.55,0.9)


## because bandwidths dont seem reasonable to me, I decided to fit the derivative again
## I used h=10 since it was an apparantly good bandwidth in the first part of Ex7
# first fit the actual curve
LPR1 <- localpoly.reg(hypage, zwast, bandwidth=10, degree.pol = 1, deriv = 0)
LPR2 <- localpoly.reg(hypage, zwast, bandwidth=10, degree.pol = 2, deriv = 0)
LPR3 <- localpoly.reg(hypage, zwast, bandwidth=10, degree.pol = 3, deriv = 0)
LPR4 <- localpoly.reg(hypage, zwast, bandwidth=10, degree.pol = 4, deriv = 0)
df <- tibble(age=LPR1$x, '1'=LPR1$predicted, '2'=LPR2$predicted,
             '3'=LPR3$predicted, '4'=LPR4$predicted)
df <- df %>%gather('1','2','3', '4', key='degree', value='predicted')%>% group_by(degree)
g2 <- ggplot(data=df)+
  geom_line(aes(x=age, y=predicted, color=degree))+
  guides(color='none')+ylim(-0.55,0.9)
# then the derivative
LPR1 <- localpoly.reg(hypage, zwast, bandwidth=10, degree.pol = 1, deriv = 1)
LPR2 <- localpoly.reg(hypage, zwast, bandwidth=10, degree.pol = 2, deriv = 1)
LPR3 <- localpoly.reg(hypage, zwast, bandwidth=10, degree.pol = 3, deriv = 1)
LPR4 <- localpoly.reg(hypage, zwast, bandwidth=10, degree.pol = 4, deriv = 1)
df <- tibble(age=LPR1$x, '1'=LPR1$predicted, '2'=LPR2$predicted,
             '3'=LPR3$predicted, '4'=LPR4$predicted)
df <- df %>%gather('1','2','3', '4', key='degree', value='predicted')%>% group_by(degree)
g3 <- ggplot(data=df)+
  geom_line(aes(x=age, y=predicted, color=degree))+
  guides(color='none')+ylim(-0.55,0.9)
ggsave('deriv.png', width=5.5, height=9, arrangeGrob(g1, g2, g3, nrow=3, heights = c(1,0.85,0.85)))
## in the last plot we see slightly positive derivatives --> improvement? Yes, a little bit. 