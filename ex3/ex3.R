  ####      Exercise 3      #####
  
  ### Preamble ####


library(xtable)
library("tidyverse")
library(gridExtra)
library(bootstrap)


  ###    a)  Bootstrap-CI for StdDev and median  ####
  ###     fix size parameters for the Weibull distr. (and seed)
set.seed(2210)
lambda <- 13 # scale parameter
k <- 1 # shape parameter

# following are the true sd and true median
tsd <- lambda*sqrt(gamma(1+2/k)-gamma(1*1/k)^2)
tmed <- lambda*(log(2)^(1/k))


## we first use the simple quantile method (and do it by hand)
# therefore define  usefull functions to shorten the code:
Boots <- function (Dat, R,n) {      # this bootstraps the StdDev and median for a given sample 
  boot_sd <- NULL
  boot_med <- NULL
  for (l in 1:R){
    temp <- sample(Dat, n, replace=TRUE)   # we use an n out of n method with replacement
    boot_sd[l] <- sd(temp)
    boot_med[l] <- median(temp)
  }
  return(data.frame(boot_sd, boot_med))   # return bootstrap-values for sd and median in a df
}

MC_Simu_boot <- function(M,R,n) {   # this does a MC-simulation for the bootstrap
  MC_sd <- rep(0,M)
  MC_med <- rep(0,M)
  length_CI_sd <- rep(0,M)
  length_CI_med <- rep(0,M)
  cover_sd <- rep(0,M)
  cover_med <- rep(0, M)
  for (l in 1:M) {
    Temp <- rweibull(n,k,lambda)    # get a new sample of Weibull-r.v.
    MC_sd[l] <- sd(Temp)            # store sd and median
    MC_med[l] <- median(Temp)
    temp_Boots <- Boots(Temp, R,n)  # do the bootstrap
    temp_CI_sd <- quantile(temp_Boots$boot_sd, c(0.025, 0.975))
    temp_CI_med <- quantile(temp_Boots$boot_med, c(0.025,0.975))
    cover_sd[l] <- tsd < temp_CI_sd[2] & tsd > temp_CI_sd[1]     # save coverage 
    cover_med[l] <- tmed < temp_CI_med[2] & tmed > temp_CI_med[1]
    length_CI_sd[l] <- temp_CI_sd[2]-temp_CI_sd[1]       # store length of CI's
    length_CI_med[l] <- temp_CI_med[2]-temp_CI_med[1]
    
  }
  return(data.frame(MC_sd, length_CI_sd, MC_med, length_CI_med, cover_sd, cover_med))  # return all values in a df
}

# now start the first trial with R = 1000 and n=100
n <- 100
R <- 1000
M <- 1000

Dat <- rweibull(n,k,lambda)             # first we draw a Weibull sample
Boot_Dat <- Boots(Dat, R, n)            # draw bootstrap samples out of this
CI_sd <- quantile(Boot_Dat$boot_sd, c(0.025, 0.975))
CI_med <- quantile(Boot_Dat$boot_med, c(0.025,0.975))   # and save the bootstrap CIs for this data

# Simulate coverage with Monte Carlo and compute average CI-length
MC_Simu1 <- MC_Simu_boot(M,R,n)     # use MC simulation and save relevant information
cov_sd1 <- sum(MC_Simu1$cover_sd)/M  # estimate coverage probability
cov_med1 <- sum(MC_Simu1$cover_med)/M
avg_CIl_sd1 <- mean(MC_Simu1$length_CI_sd)   # and mean CI-lengths
avg_CIl_med1 <- mean(MC_Simu1$length_CI_med)


# some nice plots ####
# we mark the bootstrap CIs of the first drawn sample there in red
sum(MC_Simu1$MC_sd < CI_sd[2] & MC_Simu1$MC_sd > CI_sd[1])/M  # estimate coverage probability
sum(MC_Simu1$MC_med < CI_med[2] & MC_Simu1$MC_med > CI_med[1])/M

CIs <- data.frame(x1=CI_sd, x2=CI_sd, y1=c(0,0), y2=c(170,170))
gsd <- ggplot(data=MC_Simu1)+
  geom_histogram(aes(MC_sd), alpha=0.8, linetype=1, fill='gray', colour='black', bins=25)+
  geom_segment(data=CIs, aes(x = x1, xend=x2, y=y1, yend=y2, color='segment'))+
  guides(color='none')+ylim(c(0,170))+
  xlab('SD')+ylab('frequency')+
  theme_minimal()

CIs <- data.frame(x1=CI_med, x2=CI_med, y1=c(0,0), y2=c(170,170))
gmed<- ggplot(data=MC_Simu1)+
  geom_histogram(aes(MC_med), alpha=0.8, linetype=1, fill='gray', colour='black', bins=25)+
  geom_segment(data=CIs, aes(x = x1, xend=x2, y=y1, yend=y2, color='segment'))+
  guides(color='none')+ylim(c(0,170))+
  xlab('median')+
  list(theme(axis.text.y = element_blank(),
             axis.ticks.y = element_blank(),
             axis.title.y = element_blank()))+
  theme_minimal()

CIs <- data.frame(x1=(CI_sd[2]-CI_sd[1]), x2=(CI_sd[2]-CI_sd[1]), y1=0, y2=170)
gCILsd <- ggplot(data=MC_Simu1)+
  geom_histogram(aes(length_CI_sd), alpha=0.8, linetype=1, fill='gray', colour='black', bins=25)+
  geom_segment(data=CIs, aes(x = x1, xend=x2, y=y1, yend=y2, color='segment'))+
  guides(color='none')+ylim(c(0,170))+
  xlab('Length of CI for SD')+ylab('frequency')+
  theme_minimal()
CIs <- data.frame(x1=(CI_med[2]-CI_med[1]), x2=(CI_med[2]-CI_med[1]), y1=0, y2=170)
gCILmed<- ggplot(data=MC_Simu1)+
  geom_histogram(aes(length_CI_med), alpha=0.8, linetype=1, fill='gray', colour='black', bins=25)+
  geom_segment(data=CIs, aes(x = x1, xend=x2, y=y1, yend=y2, color='segment'))+
  guides(color='none')+ylim(c(0,170))+
  xlab('Length of CI for median')+
  list(theme(axis.text.y = element_blank(),
             axis.ticks.y = element_blank(),
             axis.title.y = element_blank()))+
  theme_minimal()
p <- arrangeGrob(gsd, gmed,gCILsd, gCILmed, nrow=2, widths=c(1, 0.85))
ggsave('MChist.png', width=5, height = 5,p)


## now check what changes with the size of R and n  ####

R <- 1000
n <- 1000
MC_Simu <- MC_Simu_boot(M,R,n)
cov_sd2 <- sum(MC_Simu$cover_sd)/M  
cov_med2 <- sum(MC_Simu$cover_med)/M
avg_CIl_sd2 <- mean(MC_Simu$length_CI_sd)
avg_CIl_med2 <- mean(MC_Simu$length_CI_med)

R <- 5000
n <- 100
MC_Simu <- MC_Simu_boot(M,R,n)
cov_sd3 <- sum(MC_Simu$cover_sd)/M  
cov_med3 <- sum(MC_Simu$cover_med)/M
avg_CIl_sd3 <- mean(MC_Simu$length_CI_sd)
avg_CIl_med3 <- mean(MC_Simu$length_CI_med)


r <- c(1000, 1000, 5000)
N <- c(100, 1000, 100)
Cover_med <- c(cov_med1,cov_med2,cov_med3)
Cover_sd <- c(cov_sd1,cov_sd2,cov_sd3)
AVG_CIL_med <- c(avg_CIl_med1, avg_CIl_med2, avg_CIl_med3)
AVG_CIL_sd <- c(avg_CIl_sd1, avg_CIl_sd2, avg_CIl_sd3)
results <- data.frame(R=r, n=N, Cover_med=Cover_med, Cover_sd=Cover_sd, CI_length_med = AVG_CIL_med, CI_length_sd = AVG_CIL_sd)
results
xtable(results, auto = TRUE, digits = 2)


## now using bcanon ####
R <- 1000
n <- 100 

cov_sd4 <- NULL
cov_med4 <- NULL
z0_med <- NULL
a_med <- NULL 
z0_sd <- NULL
a_sd <- NULL
length_CI_sd4 <- NULL
length_CI_med4 <- NULL
for (l in 1:M){                     # for the number of MC samples
  Temp <- rweibull(n, k, lambda)    # get a new sample of Weibull-r.v.
  BCa_med <- bcanon(Dat, R, median, alpha=c(0.025,0.975))  #bcanon for median
  z0_med[l] <- BCa_med$z0; a_med[l] <- BCa_med$acc         # store relevant variables
  length_CI_med4[l] <- BCa_med$confpoints[2,2]- BCa_med$confpoints[1,2]
  cov_med4[l] <- (tmed<BCa_med$confpoints[2,2])&(tmed > BCa_med$confpoints[1,2])
  BCa_sd <- bcanon(Dat, R, sd, alpha=c(0.025,0.975))      #bcanon for sd
  z0_sd[l] <- BCa_sd$z0; a_sd[l] <- BCa_sd$acc
  length_CI_sd4[l] <- BCa_sd$confpoints[2,2]- BCa_sd$confpoints[1,2]
  cov_sd4[l] <- (tsd<BCa_sd$confpoints[2,2])&(tsd > BCa_sd$confpoints[1,2])
}
MC_bcanon <- data.frame(cov_med = cov_med4, cov_sd = cov_sd4, z0med = z0_med, z0sd = z0_sd, 
                        amed=a_med, asd=a_sd, CIlength_med = length_CI_med4, CIlength_sd = length_CI_sd4)
cov_med4 <- mean(cov_med4)
cov_sd4 <- mean(cov_med4)
z0_med <- mean(z0_med)
z0_sd <- mean(z0_sd)
a_med <- mean(a_med)
a_sd <- mean(a_med)
length_CI_med4 <- mean(length_CI_med4)
length_CI_sd4 <- mean(length_CI_sd4)


#combine with previous results
results2 <- c(R, n, cov_med4, cov_sd4, length_CI_med4, length_CI_sd4)
results <- rbind(results, results2)
results

z0_med;z0_sd;a_med;a_sd
pnorm(z0_med)
pnorm(z0_sd)


### b) rdi4p of shhs1-data  ####

dat <- read.csv('shhs1.txt', sep='\t', header=TRUE)
rdi4p <- dat$rdi4p

# plot histogram of the data and fitted exponential density
ggplot()+
  geom_histogram(data=dat, aes(rdi4p, stat(density)), alpha=0.8, linetype=1, fill='gray', colour='black', bins=30)+
  stat_function(geom='line', aes(x = 0:150), n = 200, fun =
                  dexp, args = list(rate=1/mean(rdi4p)), colour='orange')+
  guides(color='none')+#ylim(c(0,170))+
  xlab('rdi4p')+ylab('frequency')+
  theme_minimal()
ggsave('hist_data.png', width=3, height =3)

#compute some sample statistics
mean(rdi4p)
1/mean(rdi4p)
length(rdi4p)
median(rdi4p); sd(rdi4p)

## now do bootstrap CIs using the simple percentile method and bcanon:
bot_b <- Boots(rdi4p, 1000, length(rdi4p))
CI_sd_b <- quantile(bot_b$boot_sd, c(0.025, 0.975))
CI_med_b <- quantile(bot_b$boot_med, c(0.025,0.975))
BCa_med_b <- bcanon(rdi4p, 1000, median)
BCa_sd_b <- bcanon(rdi4p, 1000, sd)
results_b <- data.frame(lower_med = c(CI_med_b[1],BCa_med_b$confpoints[1,2]),
                        upper_med=c(CI_med_b[2],BCa_med_b$confpoints[2,2]),
                        lower_sd = c(CI_sd_b[1],BCa_sd_b$confpoints[1,2]),
                        upper_sd = c(CI_sd_b[2],BCa_sd_b$confpoints[2,2]))
row.names(results_b) <- c('bootstrap percentile', 'bootstrap acc bias-corr')
results_b
xtable(results_b, auto = TRUE, digits = 2)

#bcanon produces NaNs for the median becaue the value occurs twice
sort(rdi4p)[(length(rdi4p)/2-4):(length(rdi4p)/2+4)]
median(rdi4p)
BCa_sd_b$acc