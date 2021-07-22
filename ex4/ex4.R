####      Exercise 4      #####

### Preamble ####

library(xtable)
library("tidyverse")
library(survival)
library(ggfortify)
library(gridExtra)

###   a) Survivor function for whole sample   ####

###     load data set and select relevant variables
dat <- read.table('Thoracic.txt', header=FALSE)
dat <- dat %>% select(14, 16, 17)
names(dat) <- c('PRE30', 'AGE','Risk1Y')  
xtable(head(dat), auto=TRUE,digits = 2)

# get to know survfit.object  
?survfit.object

###    a)  Survival analysis for   ####
dead_time <- dat$AGE+1
dead <- dat$Risk1Y

fit_KM <- survfit(Surv(AGE, Risk1Y)~1, data = dat, type='kaplan-meier')
fit_FH <- survfit(Surv(AGE, Risk1Y)~1, data = dat, type='fleming-harrington')
dfKM <- tibble(age=fit_KM$time, fit=fit_KM$surv, upper=fit_KM$upper, lower=fit_KM$lower, method=factor('KM'))
dfFH <- tibble(age=fit_KM$time, fit=fit_FH$surv, upper=fit_FH$upper, lower=fit_FH$lower, method=factor('FH'))
df <- bind_rows(dfKM, dfFH)
df <- df %>% gather('fit', 'upper', 'lower', key='curve', value='survival_probability')
df
ggplot(data=df)+
  geom_step(aes(x=age, y=survival_probability, color=method, linetype=curve))+
  scale_color_manual(values=c('KM'='blue', 'FH'='orange'))+
  labs(linetype='')+
  ylab('survival probability')+
  theme_minimal()
ggsave('survfits.png', width=5, height = 3.3)



## fit parametric models
reg_weibull <- survreg(Surv(AGE, Risk1Y)~1, data = dat, dist='weibull')
w_sh <- 1/reg_weibull$scale
w_sc <- exp(reg_weibull$coefficients[1])
dfweibull <- tibble(age=seq(21,87,length.out = 400 ), weibull_fit = 1-pweibull(age, w_sh,w_sc))

reg_exp <- survreg(Surv(AGE, Risk1Y)~1, data = dat, dist='exponential')
lambda <- exp(-reg_exp$coefficients)
dfexp <- tibble(age=seq(21,87,length.out = 400 ),exp_fit=1-pexp(age, exp(-reg_exp$coefficients)))

# gather dfKM for plotting
dfKM <- dfKM %>% gather('fit', 'upper','lower', key=curve, value=survival_probability)

ggplot()+
  geom_step(data=dfKM, aes(x=age, y=survival_probability, linetype=curve, color='KM'))+
  scale_linetype_manual(values=c('fit'=1, 'upper'=2, 'lower'=3))+
  geom_step(data=dfweibull, aes(age, weibull_fit, color='Weib'))+
  geom_step(data=dfexp, aes(age, exp_fit,color='exp'))+
  guides(linetype='none')+
  scale_color_manual(values = c('KM'='black','exp' = 'orange','Weib'='blue'), 
                     labels=c('KM','exp','Weibull'))+
  labs(color='fit')+
  ylab('survival probability')+
  theme_minimal()
ggsave('paramsurvfits.png', width=5, height = 3.3)



### b) Comparison of smokers and non-smokers   ####
smokers <- subset(dat, PRE30==TRUE)
nonsmokers <- subset(dat, PRE30==FALSE)
# fit group-slitted model
KMb <- survfit(Surv(AGE, Risk1Y)~PRE30, data = dat, type='kaplan-meier')

## plot two KM fits for the groups in one plot
autoplot(KMb, censor.shape ='.', censor.alpha = 0.6)+
  labs(x = "age", y = "survival probability") +
  labs(colour = "Smoker") +
  theme_minimal()+
  scale_fill_manual(labels = c("no", "yes"), values = c('cyan2', 2), name='Smoker')+
  scale_colour_manual(labels = c("no", "yes"), values = c('darkcyan', 2), name='Smoker')
ggsave('smokers.png', width=5, height = 3.3)  

## what is the proportion of smokers in the group?
length(smokers$AGE)/length(dat$AGE) # with 82.13% rather high
length(smokers$AGE)
length(nonsmokers$AGE)

## compare groups with log-rank test
survdiff(Surv(AGE, Risk1Y)~PRE30, data = dat) # with p=0.1 it indicates a difference which would be interpreted as not significant by most researchers
nonsmokers[order(nonsmokers$AGE),]

## plot KM and Weibull fit for two groups
## here it is easier to reproduce the above plots for each group seperately:
KMsmokers <- survfit(Surv(AGE, Risk1Y)~1, data = smokers, type='kaplan-meier')
KMnonsmokers <- survfit(Surv(AGE, Risk1Y)~1, data = nonsmokers, type='kaplan-meier')

dfKMsmokers <- tibble(age=KMsmokers$time, fit=KMsmokers$surv, upper=KMsmokers$upper, lower=KMsmokers$lower, group=factor('smokers'))
dfKMnonsmokers <- tibble(age=KMnonsmokers$time, fit=KMnonsmokers$surv, upper=KMnonsmokers$upper, lower=KMnonsmokers$lower, group=factor('nonsmokers'))
dfKMsmokers <- dfKMsmokers %>% gather('fit', 'upper','lower', key=curve, value=survival_probability)
dfKMnonsmokers <- dfKMnonsmokers %>% gather('fit', 'upper','lower', key=curve, value=survival_probability)


reg_weibull_smokers <- survreg(Surv(AGE, Risk1Y)~1, data = smokers, dist='weibull')
w_sh <- 1/reg_weibull_smokers$scale
w_sc <- exp(reg_weibull_smokers$coefficients[1])
dfweibull_smokers <- tibble(age=seq(21,87,length.out = 400 ), weibull_fit = 1-pweibull(age, w_sh,w_sc))
reg_weibull_nonsmokers <- survreg(Surv(AGE, Risk1Y)~1, data = nonsmokers, dist='weibull')
w_sh <- 1/reg_weibull_nonsmokers$scale
w_sc <- exp(reg_weibull_nonsmokers$coefficients[1])
dfweibull_nonsmokers <- tibble(age=seq(21,87,length.out = 400 ), weibull_fit = 1-pweibull(age, w_sh,w_sc))



gsmoke<- ggplot()+
  geom_step(data=dfKMsmokers, aes(x=age, y=survival_probability, linetype=curve, color='KM'))+
  scale_linetype_manual(values=c('fit'=1, 'upper'=2, 'lower'=3))+
  geom_step(data=dfweibull_smokers, aes(age, weibull_fit, color='Weib'))+
  guides(linetype='none', color='none')+
  scale_color_manual(values = c('KM'='black','Weib'='blue'), 
                     labels=c('KM','Weibull'))+
  labs(color='fit')+ylim(c(0,1))+
  ylab('survival probability') + ggtitle('Smokers')+
  theme_minimal()
gnonsmoke<- ggplot()+
  geom_step(data=dfKMnonsmokers, aes(x=age, y=survival_probability, linetype=curve, color='KM'))+
  scale_linetype_manual(values=c('fit'=1, 'upper'=2, 'lower'=3))+
  geom_step(data=dfweibull_nonsmokers, aes(age, weibull_fit, color='Weib'))+
  guides(linetype='none')+
  scale_color_manual(values = c('KM'='black','Weib'='blue'), 
                     labels=c('KM','Weibull'))+
  labs(color='fit')+
  ylim(c(0,1))+
  ylab('') + ggtitle('Nonsmokers')+
  theme_minimal()
grid.arrange(gsmoke, gnonsmoke, nrow=1, widths =c(0.78,1))
g <- arrangeGrob(gsmoke, gnonsmoke, nrow=1, widths =c(0.78,1))
ggsave(file='groupsparam.png', g, width=9, height = 4)
## Weibull seems quite good