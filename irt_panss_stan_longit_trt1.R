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

# rm(list = ls())
# system("g++ --version")

require(Hmisc)
library("ltm")
library("mcmcplots")
require(plyr)
require(reshape2)
require(inline)
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


##########################################################################################
# using an array of ordered vectors
# AR(1) population-averaged model
##########################################################################################
bito_stan_long1 <- '
data {
  real m_mu_theta;               // prior for mean of theta
  real var_mu_theta;             // prior for variance of theta
  int<lower=0> T;                // number of categorical time points
  int<lower=0> N;                // number of subjects
  int<lower=0> J;                // number of items
  int<lower=1, upper=7> K[J];    // max level for each item
  int Y[N, J, T];                // data, 3D array of integers
}

parameters {
    ordered[6] kappa[J];                    // intercept : array of ordered vectors
    vector[T] theta[N];                     // ability : theta[N,T] - array of vectors
    real<lower=0.0, upper=10.0> alpha[J];   // discrimination
    real<lower=0.0, upper=1.0> rho;         // AR(1) parameter
    vector[T-1] mu_theta_raw;
}

model{
   real sigsq_theta;
   vector[T] mu_theta;
   matrix[T, T] var_theta;
   matrix[T, T] L_theta;
   real p[N, J, T, 7];
   real Q[N, J, T, 6];

   theta ~ multi_normal(mu_theta, var_theta);

// prior for mu_theta
   mu_theta[1] <- 0.0;
   for (t in 2:T){
      mu_theta[t] ~ normal(m_mu_theta, var_mu_theta);
   }

// AR(1) structure for the covariance of thetas
   sigsq_theta <- 1.0;
   var_theta[1,1] <- sigsq_theta;
   for (t in 2:T) {
      var_theta[t,t] <- sigsq_theta;
      for (j in 1:(t-1)) {
        var_theta[t,j] <- sigsq_theta*pow(rho, t-j);
        var_theta[j,t] <- var_theta[t,j];
      }
   }

  //  alpha ~ uniform(0.0, 10.0);
   for (j in 1:J) {
      for (k in 1:(K[j] - 1)) {
         kappa[j,k] ~ normal(0.0, 5.0);
      }
   }
  for (t in 1:T){
    for (i in 1:N) {
      for (j in 1:J) {
        for (k in 1:(K[j] - 1)) {
          Q[i, j, t, k] <- inv_logit(kappa[j,k] - alpha[j]*theta[i,t]);
          }
          p[i, j, t, 1] <-  Q[i, j, t, 1];
        for (k in 2:(K[j] - 1)){
          p[i, j, t, k] <- Q[i, j, t, k] - Q[i, j, t, k-1];
        }
        p[i, j, t, K[j]] <- 1 - Q[i, j, t, K[j] - 1];

      // incement log probability directly because Y[i, j]
      // has categorical distribution with varying dimension
        increment_log_prob(log(p[i, j, t, Y[i, j, t]]));
      }
    }
  }
}
'

##########################################################################################
# Compound symmetry, population-averaged model
# Treatment effect
##########################################################################################
bito_stan_long_trt1 <- '
data {
  int n_trt;                     // number of arms
  real m_mu_theta;        // prior for mean of theta
  real var_mu_theta;             // prior for variance of mean of theta
  int<lower=0> T;                // number of categorical time points
  int<lower=0> N;                // number of subjects
  int<lower=0> J;                // number of items
  int<lower=1, upper=7> K[J];    // max level for each item
  int Y[N, J, T];                // data, 3D array of integers
  int trt[N];                    // treatment arm
}

parameters {
    ordered[6] kappa[J];                    // intercept : array of ordered vectors
    vector[T] theta[N];                     // ability : theta[N,T] - array of vectors
    real<lower=0.0, upper=10.0> alpha[J];   // discrimination
    real<lower=0.0, upper=1.0> rho;         // Correlation parameter
    vector[T] mu_theta_raw[n_trt];
}

model{
   real sigsq_theta;
   vector[T] mu_theta[n_trt];
   matrix[T, T] var_theta;
   matrix[T, T] L_theta;
   real p[N, J, T, 7];
   real Q[N, J, T, 6];

   alpha ~ normal(1.0, 3.0);

   for (t in 1:T)
      for (j in 1:n_trt)
          mu_theta_raw[j, t] ~ normal(m_mu_theta, var_mu_theta);

   for (j in 1:J) {
     for (k in 1:(K[j] - 1)) {
       kappa[j,k] ~ normal(0.0, 2.5);
     }
   }

// prior for mu_theta
   mu_theta[1, 1] <- 0.0;  // fix to 0 for identifiability

// set the values of mu_theta before sampling - stan complains otherwise
for (j in 2:n_trt)
   mu_theta[j, 1] <- mu_theta_raw[j, 1];

   for (t in 2:T)
     for (j in 1:n_trt)
        mu_theta[j, t] <- mu_theta_raw[j, t];


// Compound symmetry structure for the covariance of thetas
   sigsq_theta <- 1.0;  // fix to 1 for identifiability
   var_theta[1,1] <- sigsq_theta;
   for (t in 2:T) {
      var_theta[t,t] <- sigsq_theta;
      for (j in 1:(t-1)) {
        var_theta[t,j] <- sigsq_theta*rho;
        var_theta[j,t] <- var_theta[t,j];
      }
   }

   L_theta <- cholesky_decompose(var_theta);
for (i in 1:N)
   theta[i] ~ multi_normal_cholesky(mu_theta[trt[i]], L_theta);


  for (t in 1:T){
    for (i in 1:N) {
      for (j in 1:J) {
        for (k in 1:(K[j] - 1)) {
          Q[i, j, t, k] <- inv_logit(kappa[j,k] - alpha[j]*theta[i,t]);
          }
          p[i, j, t, 1] <-  Q[i, j, t, 1];
        for (k in 2:(K[j] - 1)){
          p[i, j, t, k] <- Q[i, j, t, k] - Q[i, j, t, k-1];
        }
        p[i, j, t, K[j]] <- 1 - Q[i, j, t, K[j] - 1];

      // incement log probability directly because Y[i, j]
      // has categorical distribution with varying dimension
        increment_log_prob(log(p[i, j, t, Y[i, j, t]]));
      }
    }
  }
}
'

##########################################################################################
# Compound symmetry, population-averaged model
# Delta - treatment effect at each time point
##########################################################################################
bito_stan_long_delta1 <- '
data {
  int n_trt;                     // number of arms
  real m_mu_theta;               // prior for mean of theta
  real var_mu_theta;             // prior for variance of mean of theta
  int<lower=0> T;                // number of categorical time points
  int<lower=0> N;                // number of subjects
  int<lower=0> J;                // number of items
  int<lower=1, upper=7> K[J];    // max level for each item
  int Y[N, J, T];                // data, 3D array of integers
  int trt[N];                    // treatment arm
}

parameters {
    ordered[6] kappa[J];                    // intercept : array of ordered vectors
    vector[T] theta[N];                     // ability : theta[N,T] - array of vectors
    real<lower=0.0, upper=10.0> alpha[J];   // discrimination
    real<lower=0.0, upper=1.0> rho;         // Correlation parameter
    vector[T] mu_theta_raw[n_trt];
}

transformed parameters {
   vector[T] mu_theta[n_trt];
// prior for mu_theta
   mu_theta[1, 1] <- 0.0;  // fix to 0 for identifiability

// set the values of mu_theta before sampling - stan complains otherwise
for (j in 2:n_trt)
   mu_theta[j, 1] <- mu_theta_raw[j, 1];
for (t in 2:T)
   for (j in 1:n_trt)
      mu_theta[j, t] <- mu_theta_raw[j, t];

}

model{
   real sigsq_theta;

   matrix[T, T] var_theta;
   matrix[T, T] L_theta;
   real p[N, J, T, 7];
   real Q[N, J, T, 6];
 
   alpha ~ normal(1.0, 3.0);

   for (t in 1:T)
      for (j in 1:n_trt)
          mu_theta_raw[j, t] ~ normal(m_mu_theta, var_mu_theta);

   for (j in 1:J) {
     for (k in 1:(K[j] - 1)) {
       kappa[j,k] ~ normal(0.0, 2.5);
     }
   }

// Compound symmetry structure for the covariance of thetas
   sigsq_theta <- 1.0;  // fix to 1 for identifiability
   var_theta[1,1] <- sigsq_theta;
   for (t in 2:T) {
      var_theta[t,t] <- sigsq_theta;
      for (j in 1:(t-1)) {
        var_theta[t,j] <- sigsq_theta*rho;
        var_theta[j,t] <- var_theta[t,j];
      }
   }

   L_theta <- cholesky_decompose(var_theta);
for (i in 1:N)
   theta[i] ~ multi_normal_cholesky(mu_theta[trt[i]], L_theta);

  for (t in 1:T){
    for (i in 1:N) {
      for (j in 1:J) {
        for (k in 1:(K[j] - 1)) {
          Q[i, j, t, k] <- inv_logit(kappa[j,k] - alpha[j]*theta[i,t]);
          }
          p[i, j, t, 1] <-  Q[i, j, t, 1];
        for (k in 2:(K[j] - 1)){
          p[i, j, t, k] <- Q[i, j, t, k] - Q[i, j, t, k-1];
        }
        p[i, j, t, K[j]] <- 1 - Q[i, j, t, K[j] - 1];

      // incement log probability directly because Y[i, j]
      // has categorical distribution with varying dimension
        increment_log_prob(log(p[i, j, t, Y[i, j, t]]));
      }
    }
  }
}

generated quantities {
    vector[T-1] deltaout[n_trt-1];             // Treatment effect
// Code the treatment effect - difference in differences
for (t in 1:(T-1))
   for (l in 1:(n_trt - 1))
      deltaout[l, t] <-  mu_theta[l+1, t+1] - mu_theta[1, t+1] - (mu_theta[l+1, 1] - mu_theta[1, 1]);              
}
'


##########################################################################################
# Independent covariance, population-averaged model
##########################################################################################
bito_stan_long3 <- '
data {
  real m_mu_theta;               // hyperprior for mean of theta
  real var_mu_theta;             // hyperprior for variance of theta
  int<lower=0> T;                // number of categorical time points
  int<lower=0> N;                // number of subjects
  int<lower=0> J;                // number of items
  int<lower=1, upper=7> K[J];    // max level for each item
  int Y[N, J, T];                // data, 3D array of integers
}

parameters {
    ordered[6] kappa[J];                    // intercept : array of ordered vectors
    vector[T] theta[N];                     // ability : theta[N,T] - array of vectors
    real<lower=0.0, upper=10.0> alpha[J];   // discrimination
}

model{
   real sigsq_theta;
   vector[T] mu_theta;
   matrix[T, T] var_theta;
   matrix[T, T] L_theta;
   real p[N, J, T, 7];
   real Q[N, J, T, 6];

// Independent covariance structure for the covariance of thetas
   sigsq_theta <- 1.0;   // fix variance of first time point to 1
   var_theta[1,1] <- sigsq_theta;
   for (t in 2:T) {
      var_theta[t,t] <- sigsq_theta;
      for (j in 1:(t-1)) {
        var_theta[t,j] <- 0.0;
        var_theta[j,t] <- var_theta[t,j];
      }
   }

   L_theta <- cholesky_decompose(var_theta);

   theta ~ multi_normal_cholesky(mu_theta, L_theta);

// prior for mu_theta
   mu_theta[1] <- 0.0;  // fix mean of first time point to 0
   for (t in 2:T){
      mu_theta[t] ~ normal(m_mu_theta, var_mu_theta);
   }



//   alpha ~ uniform(0.0, 10.0);
   for (j in 1:J) {
      for (k in 1:(K[j] - 1)) {
         kappa[j,k] ~ normal(0.0, 5.0);
      }
   }
  for (t in 1:T){
    for (i in 1:N) {
      for (j in 1:J) {
        for (k in 1:(K[j] - 1)) {
          Q[i, j, t, k] <- inv_logit(kappa[j,k] - alpha[j]*theta[i,t]);
          }
          p[i, j, t, 1] <-  Q[i, j, t, 1];
        for (k in 2:(K[j] - 1)){
          p[i, j, t, k] <- Q[i, j, t, k] - Q[i, j, t, k-1];
        }
        p[i, j, t, K[j]] <- 1 - Q[i, j, t, K[j] - 1];

      // incement log probability directly because Y[i, j]
      // has categorical distribution with varying dimension
        increment_log_prob(log(p[i, j, t, Y[i, j, t]]));
      }
    }
  }
}
'

#Prepare data
panss.data <- droplevels(na.omit(subset(panss.nsfs, select = c(usubjid, visitnum, qstestcd, qsstresn, trt01a), visitnum <= 12 & visitnum >=4)))
panss.data <- panss.data[!duplicated(panss.data),]
panss.melt <- melt(panss.data, id.vars = c('usubjid', 'visitnum', 'qstestcd', 'trt01a'), measure.vars = 'qsstresn')
panss.wide <- dcast(panss.melt, formula =  usubjid + trt01a + qstestcd ~ visitnum + variable)
panss.wide <- na.omit(panss.wide)
panss.melt2 <- melt(panss.wide)
panss.wide2 <- dcast(panss.melt2, formula = usubjid + trt01a + variable ~ qstestcd)
panss.wide2 <- na.omit(panss.wide2)

X <- abind(split(panss.wide2, panss.wide2$variable), along = 3)
Y <- apply(X[, -c(1,2,3), ], c(2,3), as.numeric) 

# add 1 to make range 1-7 -> for BUGS dcat() function to work
Y <- Y+1


## Y <- aperm(Y, c(1,3,2))

# Test a subset of the data


## Y <- Y[1:50, , 1:2] # note: WinBUGS runs with this subset
## Y <- Y[ , , 1:2]    # note: WinBUGS runs with this subset
## Y <- Y[ , , 1:3]    # ran in BUGS but slow
Y <- Y[ , , 1:4]    # ran in BUGS but slow
attr(Y, "dimnames")[2][[1]] <- c("X1","X2","X3","X4","X5","X6","X7")
attr(Y, "dimnames")[3][[1]] <- c("T1","T2","T3","T4")

# save(Y, file= "simdata1.Rdata")

# Y <- Y[ , , c(1,2,5)]
# Y <- Y[ , , 1] # note: 

## # fake simulated data for error report
## Y <- rmultinom(n = 150, size = 7, prob = c(0.15,0.15,0.15,0.15,0.15,0.15,0.1))
## Y <- t(Y) + 1
## write.csv(Y, file = "simdata.txt", row.names=FALSE) 

# treatment arm
trt <- 4 - as.numeric(panss.wide2[panss.wide2$variable == '4_qsstresn', ]$trt01a)

# time
T <- dim(Y)[3]
# n.variables    
J <- dim(Y)[2]
# sample size
N <- dim(Y)[1]
# K - max level for ordinal variable for each item
K <- apply(apply(Y, c(2,3), max), 1, max)
# number of arms
n.trt <- length(levels(panss.wide2[panss.wide2$variable == '4_qsstresn', ]$trt01a))

## m.alpha <- 1.0
## s.alpha <- 2.5
## m.kappa <- 0.0
## s.kappa <- 2.5
## m.mu.theta <- 0
## s.mu.theta <- 1

m.mu.theta <- 0.0
var.mu.theta <- 1.0

# data for stan()
bito_dat <- list(T = T,
                 N = N,
                 J = J, 
                 Y = Y,
                 K = K,
                 m_mu_theta = m.mu.theta,
                 var_mu_theta = var.mu.theta,
                 trt = trt,
                 n_trt = n.trt)

# set_cppo(mode = "debug")
# set_cppo(mode = "fast")

# initial values
initf3 <- function() {
    list(theta = array(rnorm(N*T), dim=c(N,T)), alpha = runif(J), kappa = t(apply(array(rnorm(J*6), dim=c(J,6)), 1, sort)), rho= runif(1))
}

fitl1 <- stan(model_code = bito_stan_long1, data = bito_dat, 
              iter = 1000, chains = 1, verbose = FALSE)

fitl2 <- stan(model_code = bito_stan_long1, data = bito_dat, 
              iter = 1000, chains = 3, verbose = FALSE, init = initf3)

fitl3 <- stan(model_code = bito_stan_long2, data = bito_dat, 
              iter = 1000, chains = 3, verbose = FALSE, init = initf3)

fitl.trt3b <- stan(model_code = bito_stan_long_trt1, data = bito_dat, 
                   iter = 1000, chains = 3, verbose = FALSE)
# get the initial values used:
fitl.trt3b.inits <- get_inits(fitl.trt3b)


##################################################
fitl.delta1 <- stan(model_code = bito_stan_long_delta1, data = bito_dat, 
                    iter = 1000, chains = 3, verbose = FALSE)
fitl.delta1.inits <- get_inits(fitl.delta1)

# update the model with more iterations:
fitl.delta2 <- stan(fit =fitl.delta1 , data = bito_dat, 
                    iter = 10000, chains = 3, verbose = FALSE)


##########################################################################################
# initial values
initf4 <- function() {
    list(theta = array(rnorm(N*T), dim=c(N,T)), alpha = runif(J), kappa = t(apply(array(rnorm(J*6), dim=c(J,6)), 1, sort)))
}

fitl4 <- stan(model_code = bito_stan_long3, data = bito_dat, 
             iter = 1000, chains = 3, verbose = FALSE, init = initf4)

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

##################################################
print(fitl3b, pars=c("theta", "alpha", "kappa", "rho"), probs=c(.1,.5,.9))
print(fitl3b, pars=c("kappa"), probs=c(.1,.5,.9))
#trace plots
rstan::traceplot(fitl3b, c("alpha"), ncol=2, inc_warmup=F)
rstan::traceplot(fitl3b, c("kappa"), ncol=12, inc_warmup=F)

##################################################
print(fitl.trt3b, pars=c("theta", "alpha", "kappa", "rho"), probs=c(.1,.5,.9))
print(fitl.trt3b, pars=c("kappa"), probs=c(.1,.5,.9))
#trace plots
rstan::traceplot(fitl.trt3b, c("alpha"), ncol=2, inc_warmup=F)
rstan::traceplot(fitl.trt3b, c("kappa"), ncol=12, inc_warmup=F)

##################################################
print(fitl.delta1, pars=c("theta", "alpha", "kappa", "rho", "deltaout"), probs=c(.1,.5,.9))
print(fitl.delta1, pars=c("alpha"), probs=c(.1,.5,.9))
print(fitl.delta1, pars=c("kappa"), probs=c(.1,.5,.9))
print(fitl.delta1, pars=c("deltaout"), probs=c(.1,.5,.9))
#trace plots
rstan::traceplot(fitl.delta1, c("alpha"), ncol=2, inc_warmup=F)
rstan::traceplot(fitl.delta1, c("kappa"), ncol=6, nrow=7, inc_warmup=F)
rstan::traceplot(fitl.delta1, c("deltaout"), ncol=2, inc_warmup=F)

