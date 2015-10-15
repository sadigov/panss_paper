##########################################################################################
#
#   IRT Analysis of Bitopertin PANSS NSFS
#   Trial n25310
#   Author: Shamil Sadikhov
#   Date 10 Aug 2014
#   GRM IRT model for visit 1 only
#   Using RSTAN package
#
##########################################################################################

rm(list = ls())
system("g++ --version")

require(Hmisc)
library("ltm")
library("mcmcplots")
require(plyr)
require(reshape2)
require(rstan)
require(abind)

setwd('u:/My Documents/Statistics_projects/IRT/bitopertin')
Sys.setlocale("LC_TIME", "English") # with incorrect locale as.Date does not recognize months, months appear in German
set.seed(1234567)

# bitodata <- sasxport.get("u:/My Documents/Bitopertin/MAAP/PPNS/Analysis/n25310a/data_21112013/csv_21112013", method='csv')
# names(bitodata)
# panss <- bitodata$panss
load('panss.Rdata')
panss.nsfs <- subset(panss, qsscat %in% c('NEGATIVE SCALE', 'GENERAL PSYCHOPATHOLOGY SCALE') & qstestcd %in% c('QS027N01', 'QS027N02', 'QS027N03', 'QS027N04', 'QS027N06', 'QS027G07', 'QS027G16'))

# rm(bitodata)
# save(panss, file = 'panss.Rdata')

lgrmar1 <- function(){
# Graded response longitudinal model
# Unconstrained - discrimination parameters alpha freely estimated
    
for (t in 1:T){
    for (i in 1:n){
        for (j in 1:p){
            Y[i, j, t] ~ dcat(prob[i, j, t, 1:K[j]])

            for (k in 1:(K[j] - 1)) {
            logit(P[i, j, t, k]) <- kappa[j, k] - alpha[j]*theta[i, t]
                                   
        }
            P[i, j, t, K[j]] <- 1.0
        } 
    
        for (j in 1:p){
            prob[i, j, t, 1] <- P[i, j, t, 1]
            for (k in 2:K[j]) {
                prob[i, j, t, k] <- P[i, j, t, k] - P[i, j, t, k - 1]
            }
        }
    }
}

for (i in 1:n){
    theta[i, 1:T] ~ dmnorm(mu.theta[], Pr.theta[,])
}

## Prior for mu.theta
mu.theta[1] <- 0.0
for (t in 2:T){
    mu.theta[t] ~ dnorm(m.mu.theta, pr.mu.theta)
}

pr.mu.theta <- pow(s.mu.theta, -2)

## AR(1) structure for Sigma.theta
sigsq.theta <- 1.0
Sigma.theta[1, 1] <- sigsq.theta
for (t in 2:T){
    Sigma.theta[t, t] <- sigsq.theta
    for (j in 1:(t - 1)){
        Sigma.theta[t, j] <- sigsq.theta*pow(rho, t - j)
        Sigma.theta[j, t] <- Sigma.theta[t, j]
    }
}

Pr.theta[1:T, 1:T] <- inverse(Sigma.theta[,])
rho ~ dunif(-1.0, 1.0)

## Priors on item parameters
for (j in 1:p){
    alpha[j] ~ dnorm(m.alpha, pr.alpha) %_% I(0, )
}
pr.alpha <- pow(s.alpha, -2)

# thresholds need to be ranked
for (j in 1:p){
    for (k in 1:(K[j] - 1)){
        kappa.star [j, k] ~ dnorm(m.kappa, pr.kappa)
        kappa[j, k] <- ranked(kappa.star[j, 1:(K[j] - 1)], k)
    }
}
pr.kappa <- pow(s.kappa, -2)       
}


## some temporary filename:
filename <- file.path(getwd(), "lgrmar1.txt")
## write model file:
write.model(lgrmar1, filename)
file.show(filename)

# https://github.com/stan-dev/example-models/blob/master/bugs_examples/vol1/bones/bones.stan

bito_stan1 <- '
data {
  int<lower=0> N;                // number of subjects
  int<lower=0> J;                // number of items
  int<lower=1, upper=7> K[J];    // max level for each item
  int Y[N, J];                   // data
}

parameters {
    matrix[J, 6] kappastar;                 // intercept
    real theta[N];                          // ability
    real<lower=0.0, upper=10.0> alpha[J];   // discrimination

}

model{
   matrix[J, 6] kappa; 
   real p[N, J, 7];
   real Q[N, J, 6];
   theta ~ normal(0.0, 1.0);
  //  alpha ~ uniform(0.0, 10.0);
   for (j in 1:J) {
      for (k in 1:(K[j] - 1)) {
         kappastar[j,k] ~ normal(0.0, 50.0);
      }
 kappa[j] <- sort_asc(kappastar[j]);
   }

  for (i in 1:N) {
    for (j in 1:J) {
      for (k in 1:(K[j] - 1)) {
        Q[i, j, k] <- inv_logit(kappa[j,k] - alpha[j]*theta[i]);
        }
        p[i, j, 1] <-  Q[i, j, 1];
      for (k in 2:(K[j] - 1)){
        p[i, j, k] <- Q[i, j, k] - Q[i, j, k-1];
      }
      p[i, j, K[j]] <- 1 - Q[i, j, K[j] - 1];

      // incement log probability directly because Y[i, j]
      // has categorical distribution with varying dimension
      increment_log_prob(log(p[i, j, Y[i, j]]));  
    }
  }
}
'

# try to get kappa as a parameter
bito_stan2 <- '
data {
  int<lower=0> N;                // number of subjects
  int<lower=0> J;                // number of items
  int<lower=1, upper=7> K[J];    // max level for each item
  int Y[N, J];                   // data
}

parameters {
    matrix[J, 6] kappastar;                 // intercept
    real theta[N];                          // ability
    real<lower=0.0, upper=10.0> alpha[J];   // discrimination

}

transformed parameters {
   matrix[J, 6] kappa;
   for (j in 1:J) {
        kappa[j] <- sort_asc(kappastar[j]);
   }

}

model{

   real p[N, J, 7];
   real Q[N, J, 6];
   theta ~ normal(0.0, 1.0);
  //  alpha ~ uniform(0.0, 10.0);
   for (j in 1:J) {
      for (k in 1:(K[j] - 1)) {
         kappastar[j,k] ~ normal(0.0, 50.0);
      }
   }

  for (i in 1:N) {
    for (j in 1:J) {
      for (k in 1:(K[j] - 1)) {
        Q[i, j, k] <- inv_logit(kappa[j,k] - alpha[j]*theta[i]);
        }
        p[i, j, 1] <-  Q[i, j, 1];
      for (k in 2:(K[j] - 1)){
        p[i, j, k] <- Q[i, j, k] - Q[i, j, k-1];
      }
      p[i, j, K[j]] <- 1 - Q[i, j, K[j] - 1];

      // incement log probability directly because Y[i, j]
      // has categorical distribution with varying dimension
      increment_log_prob(log(p[i, j, Y[i, j]]));  
    }
  }
}
'

# try using an array of ordered vectors
bito_stan3 <- '
data {
  int<lower=0> N;                // number of subjects
  int<lower=0> J;                // number of items
  int<lower=1, upper=7> K[J];    // max level for each item
  int Y[N, J];                   // data
}

parameters {
    ordered[6] kappa[J];                    // intercept
    real theta[N];                          // ability
    real<lower=0.0, upper=10.0> alpha[J];   // discrimination

}

model{

   real p[N, J, 7];
   real Q[N, J, 6];
   theta ~ normal(0.0, 1.0);
  //  alpha ~ uniform(0.0, 10.0);
   for (j in 1:J) {
      for (k in 1:(K[j] - 1)) {
         kappa[j,k] ~ normal(0.0, 5.0);
      }
   }

  for (i in 1:N) {
    for (j in 1:J) {
      for (k in 1:(K[j] - 1)) {
        Q[i, j, k] <- inv_logit(kappa[j,k] - alpha[j]*theta[i]);
        }
        p[i, j, 1] <-  Q[i, j, 1];
      for (k in 2:(K[j] - 1)){
        p[i, j, k] <- Q[i, j, k] - Q[i, j, k-1];
      }
      p[i, j, K[j]] <- 1 - Q[i, j, K[j] - 1];

      // incement log probability directly because Y[i, j]
      // has categorical distribution with varying dimension
      increment_log_prob(log(p[i, j, Y[i, j]]));  
    }
  }
}
'

#Prepare data
panss.data <- droplevels(na.omit(subset(panss.nsfs, select = c(usubjid, visitnum, qstestcd, qsstresn), visitnum <= 12 & visitnum >=4)))
panss.data <- panss.data[!duplicated(panss.data),]
panss.melt <- melt(panss.data, id.vars = c('usubjid', 'visitnum', 'qstestcd'), measure.vars = 'qsstresn')
panss.wide <- dcast(panss.melt, formula =  usubjid + qstestcd ~ visitnum + variable)
panss.wide <- na.omit(panss.wide)
panss.melt2 <- melt(panss.wide)
panss.wide2 <- dcast(panss.melt2, formula = usubjid + variable ~ qstestcd)
panss.wide2 <- na.omit(panss.wide2)

X <- abind(split(panss.wide2, panss.wide2$variable), along = 3)
Y <- apply(X[, -c(1,2), ], c(2,3), as.numeric) 

# add 1 to make range 1-7 -> for BUGS dcat() function to work
Y <- Y+1
## Y <- aperm(Y, c(1,3,2))

# Test a subset of the data


Y <- Y[1:50, , 1:2] # note: WinBUGS runs with this subset
Y <- Y[ , , 1:2]    # note: WinBUGS runs with this subset
Y <- Y[ , , 1:3]    # ran in BUGS but slow
Y <- Y[ , , 1:4]    # ran in BUGS but slow


# Y <- Y[ , , c(1,2,5)]
Y <- Y[ , , 1] # note: 

## # fake simulated data for error report
## Y <- rmultinom(n = 150, size = 7, prob = c(0.15,0.15,0.15,0.15,0.15,0.15,0.1))
## Y <- t(Y) + 1
## write.csv(Y, file = "simdata.txt", row.names=FALSE) 

# time
T <- dim(Y)[3]
# n.variables    
J <- dim(Y)[2]
# sample size
N <- dim(Y)[1]
# K - max level for ordinal variable for each item
# K <- apply(apply(Y, c(2,3), max), 1, max)
K <- apply(Y, c(2), max)

m.alpha <- 1.0
s.alpha <- 2.5
m.kappa <- 0.0
s.kappa <- 2.5
m.mu.theta <- 0
s.mu.theta <- 1

# data for stan()
bito_dat <- list(N = N,
                 J = J, 
                 Y = Y,
                 K = K)

# set_cppo(mode = "debug")
# set_cppo(mode = "fast") 
fit1 <- stan(model_code = bito_stan1, data = bito_dat, 
             iter = 5000, chains = 3, verbose = FALSE)

fit2 <- stan(model_code = bito_stan1, data = bito_dat, 
             iter = 5000, chains = 3, verbose = FALSE, init = '"0.5"' )

# function form 1 without arguments 
initf1 <- function() {
    list(theta = rnorm(N),  alpha = runif(J), kappastar = t(apply(array(rnorm(J*6), dim=c(J,6)), 1, sort)))
}

initf2 <- function() {
    list(theta = rnorm(N),  alpha = runif(J), kappa = t(apply(array(rnorm(J*6), dim=c(J,6)), 1, sort)))
}


fit3 <- stan(model_code = bito_stan1, data = bito_dat, 
             iter = 5000, chains = 3, verbose = FALSE, init = initf1 )

fit4 <- stan(model_code = bito_stan2, data = bito_dat, 
             iter = 5000, chains = 3, verbose = FALSE, init = initf1 )

fit5 <- stan(model_code = bito_stan3, data = bito_dat, 
             iter = 5000, chains = 3, verbose = FALSE, init = initf2 )

print(fit1, pars=c("theta", "alpha", "kappastar"), probs=c(.1,.5,.9))

fit1.res <- extract(fit1)

#trace plots
rstan::traceplot(fit1, c("alpha"), ncol=2, inc_warmup=F)
rstan::traceplot(fit1, c("kappastar"), ncol=12, inc_warmup=F)


print(fit3, pars=c("theta", "alpha", "kappastar"), probs=c(.1,.5,.9))

fit3.res <- extract(fit3)

#trace plots
rstan::traceplot(fit3, c("alpha"), ncol=2, inc_warmup=F)
rstan::traceplot(fit3, c("kappastar"), ncol=12, inc_warmup=F)


print(fit4, pars=c("theta", "alpha", "kappastar"), probs=c(.1,.5,.9))
#trace plots
rstan::traceplot(fit4, c("alpha"), ncol=2, inc_warmup=F)
rstan::traceplot(fit4, c("kappa"), ncol=12, inc_warmup=F)


print(fit5, pars=c("theta", "alpha"), probs=c(.1,.5,.9))
#trace plots
rstan::traceplot(fit4, c("alpha"), ncol=2, inc_warmup=F)
rstan::traceplot(fit4, c("kappa"), ncol=12, inc_warmup=F)
