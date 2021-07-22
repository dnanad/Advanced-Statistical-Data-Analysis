####      Exercise 7      #####

### Preamble ####
library(xtable)
library("tidyverse")

###     load data set
dat <- read.csv('student-mat.csv', sep=',', header=TRUE)

###    a)  Distribution of G1, G2, G3  ####

x <- seq(0,20,by=1)
xa <- seq(-0.5, 20.5, by=1)
grades <- select(dat, c('G1','G2','G3'))
#grades <- grades %>%
#grades
tablesummarise <- grades %>% gather('G1', 'G2', 'G3',key='timepoint', value = 'grade') %>%
  group_by(timepoint) %>% summarise(min= min(grade),
                                    q25 = quantile(grade, 0.25),
                                    med = median(grade),
                                    M =mean(grade),
                                    q75 = quantile(grade, 0.75),
                                    max=max(grade),
                                    IQR = q75-q25,
                                    SD = sd(grade),
                                    F = var(grade)/mean(grade))
xtable(tablesummarise, auto=TRUE, digits=2)


x11()
#par(mfrow=c(3,2))
layout(matrix(c(1,4,7,2,5,8,3,6,9), nrow = 3, byrow = TRUE))
disp <- NULL
for (k in 1:3){
  G <- grades[,k]
  hist(G, breaks=xa, freq=FALSE, ylim=c(0,0.15), main=paste('G', k), xlab='Grade', cex.main=1.8)
  curve(dnorm(x, mean(G), sd(G)), col='red', add=TRUE)
  points(x, dpois(x, mean(G)), col='blue', type='p')
  qqnorm(G, main='Normal QQ-Plot')
  qqline(G)
  thpois <- qpois(seq(0,1,length.out = length(G)), lambda = mean(G))
  qqplot(thpois,G,  xlab='Theoretical Quantiles', ylab='Sample Quantiles', main='Poisson QQ-Plot')
  abline(0,1)
  disp[k] <- var(G)/mean(G)
}
disp

## drop 0s from G3 and plot the same stuff
G3 <- subset(grades, G3!=0)$G3
G3
x11()
par(mfrow=c(1,3))
hist(G3, breaks=xa, freq=FALSE, ylim=c(0,0.2), main='G3 withouth zeros', xlab='Grade', cex.main=1.8)
curve(dnorm(x, mean(G3), sd(G3)), col='red', add=TRUE)
points(x, dpois(x, mean(G3)), col='blue', type='p')
qqnorm(G3, main='Normal QQ-Plot')
qqline(G3)
thpois <- qpois(seq(0,1,length.out = length(G3)), lambda = mean(G3))
qqplot(thpois,G3,  xlab='Theoretical Quantiles', ylab='Sample Quantiles', main='Poisson QQ-Plot')
abline(0,1)
disp[4] <- var(G3)/mean(G3)
disp

### b) Fit GLMs for G1     ####

X <- dat %>% select(-G1,-G2,-G3) #extract explanatory variables
fmla <- as.formula(paste('G1 ~ ',paste(names(X), collapse = '+')))
model1 <- glm(fmla, data=dat,family='poisson' )
summary(model1)
model0 <- glm(G1~1, data=dat, family='poisson')
anova(model1, model0, test='Chisq')

plot(model1$residuals)
x11()
par(mfrow=c(2,2))
plot(model1)

## examine goodness-of-fit with residuals ####

## compute Pearson residuals
res_pears<-residuals(model1, type = 'pearson')
# compute chisq statistic
sum(res_pears^2)

##compute Ascombe's residuals
res_Asc <- 3*((dat$G1)^(2/3)-model1$fitted.values^(2/3))/(2*model1$fitted.values^(1/6))


x11()
par(mfrow=c(2, 3))

plot(model1$fitted.values, res_pears, main='Fitted against residuals', ylab='Pearson residuals', xlab='', cex.lab=1.3)
qqnorm(res_pears, ylab='Quantiles of resdiuals', xlab='')
hist(res_pears, freq=FALSE, main='Histogram', xlab='Pearson residuals', ylim=c(0,0.52))
curve(dnorm(x, 0, sd(res_pears)), add=TRUE, col='red')

plot(model1$fitted.values,res_Asc, ylab = 'Anscombe residuals', xlab ='Fitted values', cex.lab =1.3)
qqnorm(res_Asc, main='', ylab='Quantiles of resdiuals')
hist(res_Asc, freq=FALSE, main='', xlab='Anscombe residuals', ylim=c(0,0.52))
curve(dnorm(x, 0, sd(res_Asc)), add=TRUE, col='red')

# ## this was just a little playing around for examining the stationarity of residuals
# s=10
# x <- seq(min(model1$fitted.values), max(model1$fitted.values), length.out = s)
# cent <- (x[1:(s-1)]+x[2:s])/2
# SD <- NULL
# for (i in 1:(length(x)-1)){
#   temp <- data.frame(values = model1$fitted.values) %>% filter(((values<x[i+1]) & (values >x[i])))
#   SD[i] <- sd(temp$values)
# }
# plot(model1$fitted.values,res_pears, ylim=c(-0.5,0.6))
# points(cent, SD,type='l', col='red' )


###  c) Comparison of nested models ####

model2 <- update(model1, .~sex+Fedu+studytime+failures+schoolsup+famsup+goout)
summary(model2)  
model2$model
# look at https://www.kaggle.com/uciml/student-alcohol-consumption for interpretation

an_o_dev <- anova(model1, model2, test='Chisq')
an_o_dev


## examine goodness-of-fit for model2 with residuals ####
## compute Pearson residuals
res_pears2<-residuals(model2, type = 'pearson')
##compute Ascombe's residuals
res_Asc2 <- 3*((dat$G1)^(2/3)-model2$fitted.values^(2/3))/(2*model2$fitted.values^(1/6))

#x11()
par(mfrow=c(2, 3))

plot(model2$fitted.values, res_pears2, main='Fitted against residuals', ylab='Pearson residuals', xlab='', cex.lab=1.3)
qqnorm(res_pears2, ylab='Quantiles of resdiuals', xlab='')
hist(res_pears2, freq=FALSE, main='Histogram', xlab='Pearson residuals', ylim=c(0,0.52))
curve(dnorm(x, 0, sd(res_pears2)), add=TRUE, col='red')

plot(model2$fitted.values,res_Asc2, ylab = 'Anscombe residuals', xlab ='Fitted values', cex.lab =1.3)
qqnorm(res_Asc2, main='', ylab='Quantiles of resdiuals')
hist(res_Asc2, freq=FALSE, main='', xlab='Anscombe residuals', ylim=c(0,0.52))
curve(dnorm(x, 0, sd(res_Asc2)), add=TRUE, col='red')



## now fit a third model3     ####
model3 <- update(model2, .~.-goout+Walc)
summary(model3)

## examine goodness-of-fit for model3 with residuals ####
## compute Pearson residuals
res_pears3<-residuals(model3, type = 'pearson')
##compute Ascombe's residuals
res_Asc3 <- 3*((dat$G1)^(2/3)-model3$fitted.values^(2/3))/(2*model3$fitted.values^(1/6))
x11()
par(mfrow=c(2, 3))
plot(model3$fitted.values, res_pears3, main='Fitted against residuals', ylab='Pearson residuals', xlab='', cex.lab=1.3)
qqnorm(res_pears3, ylab='Quantiles of resdiuals', xlab='')
hist(res_pears3, freq=FALSE, main='Histogram', xlab='Pearson residuals', ylim=c(0,0.53))
curve(dnorm(x, 0, sd(res_pears3)), add=TRUE, col='red')
plot(model3$fitted.values,res_Asc3, ylab = 'Anscombe residuals', xlab ='Fitted values', cex.lab =1.3)
qqnorm(res_Asc3, main='', ylab='Quantiles of resdiuals')
hist(res_Asc3, freq=FALSE, main='', xlab='Anscombe residuals', ylim=c(0,0.53))
curve(dnorm(x, 0, sd(res_Asc3)), add=TRUE, col='red')


## comparison by a Chisq-distr not possible but one compare the deviances:
c(model2$deviance, model3$deviance)
# or the likelihood
logLik(model2)
logLik(model3)
