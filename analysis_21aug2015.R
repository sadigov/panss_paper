require(Hmisc)
require(rstan)
require(rstanmulticore)
# require(eha)
# require(chron)
# require(sqldf)
# require(tseries)
# library(plyr)
require(reshape2)
require(abind)
require(ggplot2)
# require(msm)
# require(equate)

Sys.setlocale("LC_TIME", "English") # with incorrect locale as.Date does not recognize months, months appear in German
setwd("u:/My Documents/Statistics_projects/IRT/bitopertin/")
set.seed(12345)

################################################################################
#                          Import data
################################################################################

# bitodata <- sasxport.get("u:/My Documents/Statistics_projects/IRT/bitopertin/csvdata", method = 'csv')
                                        # names(bitodata)

bitodata_all <- sasxport.get("H:\\cdt4715l\\libraries\\csvdata\\all", method = 'csv')

panssitem_all <- subset(bitodata_all, avisit %in% c('Baseline visit', 'Week 4', 'Week 8', 'Week 12', 'Week 16', 'Week 18', 'Week 20', 'Week 24', 'Week 24 LOCF') & qstestcd %in% c('QS027G01', 'QS027G02', 'QS027G03', 'QS027G04', 'QS027G05', 'QS027G06', 'QS027G07', 'QS027G08', 'QS027G09', 'QS027G10', 'QS027G11', 'QS027G12', 'QS027G13', 'QS027G14', 'QS027G15', 'QS027G16', 'QS027N01', 'QS027N02', 'QS027N03' ,'QS027N04', 'QS027N05', 'QS027N06', 'QS027N07', 'QS027P01', 'QS027P02', 'QS027P03', 'QS027P04', 'QS027P05', 'QS027P06', 'QS027P07') & randfl =='Y')


panssitem <- subset(bitodata, avisit %in% c('Baseline visit', 'Week 4', 'Week 8', 'Week 12', 'Week 16', 'Week 18', 'Week 20', 'Week 24', 'Week 24 LOCF') & qstestcd %in% c('QS027G01', 'QS027G02', 'QS027G03', 'QS027G04', 'QS027G05', 'QS027G06', 'QS027G07', 'QS027G08', 'QS027G09', 'QS027G10', 'QS027G11', 'QS027G12', 'QS027G13', 'QS027G14', 'QS027G15', 'QS027G16', 'QS027N01', 'QS027N02', 'QS027N03' ,'QS027N04', 'QS027N05', 'QS027N06', 'QS027N07', 'QS027P01', 'QS027P02', 'QS027P03', 'QS027P04', 'QS027P05', 'QS027P06', 'QS027P07') & randfl =='Y')

# save(bitodata, file = 'bitodata_21aug2015.RData')
# load(file='bitodata_21aug2015.RData')
# panssitem <- subset(bitodata, avisit %in% c('Baseline visit', 'Week 4', 'Week 8', 'Week 12', 'Week 16', 'Week 18', 'Week 20', 'Week 24', 'Week 24 LOCF') & qstestcd %in% c('QS027G01', 'QS027G02', 'QS027G03', 'QS027G04', 'QS027G05', 'QS027G06', 'QS027G07', 'QS027G08', 'QS027G09', 'QS027G10', 'QS027G11', 'QS027G12', 'QS027G13', 'QS027G14', 'QS027G15', 'QS027G16', 'QS027N01', 'QS027N02', 'QS027N03' ,'QS027N04', 'QS027N05', 'QS027N06', 'QS027N07', 'QS027P01', 'QS027P02', 'QS027P03', 'QS027P04', 'QS027P05', 'QS027P06', 'QS027P07') & randfl =='Y')

# Use armcd variable for study arm, use avisit variable for time

                    
# rm(bitodata)
# save(panssitem, file = 'panssitem_21aug2015.RData')
# save(panssitem_all, file = 'panssitem_all_21aug2015.RData')
# load(file='panssitem_21aug2015.RData')
# load(file='panssitem_all_21aug2015.RData')
panssitem$trt01a <- relevel(panssitem$trt01a, ref = "Placebo" )
panssitem_all$trt01a <- relevel(panssitem_all$trt01a, ref = "Placebo" )

################################################################################
                                        #Prepare data
################################################################################
# N0 LOCF data
# PANSS negative symptom score
panss_ns <- droplevels(na.omit(subset(panssitem, select = c(usubjid, avisitn, avisit, qstestcd, qsstresn, trt01a), qstestcd %in% c('QS027N01', 'QS027N02', 'QS027N03' ,'QS027N04', 'QS027N05', 'QS027N06', 'QS027N07') & avisitn != 10.5 & trt01a %in% c('Bitopertin 10mg', 'Bitopertin 20mg', 'Bitopertin 5mg', 'Placebo'))))
panss_ns <- panss_ns[!duplicated(panss_ns),]
panss.melt <- melt(panss_ns, id.vars = c('usubjid', 'avisitn', 'qstestcd', 'trt01a'), measure.vars = 'qsstresn')

panss.melt <- panss.melt[!duplicated(subset(panss.melt, select = c('usubjid', 'avisitn', 'qstestcd', 'trt01a', 'variable'))), ]
panss.wide <- dcast(panss.melt, formula =  usubjid + trt01a + qstestcd ~ avisitn + variable)

panss.wide <- na.omit(panss.wide)
panss.melt2 <- melt(panss.wide)
panss.wide2 <- dcast(panss.melt2, formula = usubjid + trt01a + variable ~ qstestcd)
panss.wide2 <- na.omit(panss.wide2)

X <- abind(split(panss.wide2, panss.wide2$variable), along = 3)
Y <- apply(X[, -c(1,2,3), ], c(2,3), as.numeric) + 1

# PANSS negative symptom factor score
# Study 310
panss_nsfs <- droplevels(na.omit(subset(panssitem, select = c(usubjid, avisitn, avisit, qstestcd, qsstresn, trt01a), qstestcd %in% c('QS027N01', 'QS027N02', 'QS027N03', 'QS027N04', 'QS027N06', 'QS027G07','QS027G16') & avisitn != 10.5 & trt01a %in% c('Bitopertin 10mg', 'Bitopertin 20mg', 'Bitopertin 5mg', 'Placebo') & studyid == 'NN25310')))
panss_nsfs <- panss_nsfs[!duplicated(panss_nsfs),]
panss_nsfs.melt <- melt(panss_nsfs, id.vars = c('usubjid', 'avisitn', 'qstestcd', 'trt01a'), measure.vars = 'qsstresn')

panss_nsfs.melt <- panss_nsfs.melt[!duplicated(subset(panss_nsfs.melt, select = c('usubjid', 'avisitn', 'qstestcd', 'trt01a', 'variable'))), ]
panss_nsfs.wide <- dcast(panss_nsfs.melt, formula =  usubjid + trt01a + qstestcd ~ avisitn + variable)

panss_nsfs.wide <- na.omit(panss_nsfs.wide)
panss_nsfs.melt2 <- melt(panss_nsfs.wide)
panss_nsfs.wide2 <- dcast(panss_nsfs.melt2, formula = usubjid + trt01a + variable ~ qstestcd)
panss_nsfs.wide2 <- na.omit(panss_nsfs.wide2)

X <- abind(split(panss_nsfs.wide2, panss_nsfs.wide2$variable), along = 3)
Y <- apply(X[, -c(1,2,3), ], c(2,3), as.numeric) + 1

# PANSS negative symptom factor score
# All Studies
panss_nsfs <- droplevels(na.omit(subset(panssitem, select = c(studyid, usubjid, avisitn, avisit, qstestcd, qsstresn, trt01a), qstestcd %in% c('QS027N01', 'QS027N02', 'QS027N03', 'QS027N04', 'QS027N06', 'QS027G07','QS027G16') & avisitn != 10.5 & trt01a %in% c('Bitopertin 10mg', 'Bitopertin 20mg', 'Bitopertin 5mg', 'Placebo'))))
panss_nsfs <- panss_nsfs[!duplicated(panss_nsfs),]
panss_nsfs.melt <- melt(panss_nsfs, id.vars = c('studyid', 'usubjid', 'avisitn', 'qstestcd', 'trt01a'), measure.vars = 'qsstresn')

panss_nsfs.melt <- panss_nsfs.melt[!duplicated(subset(panss_nsfs.melt, select = c('studyid', 'usubjid', 'avisitn', 'qstestcd', 'trt01a', 'variable'))), ]
panss_nsfs.wide <- dcast(panss_nsfs.melt, formula = studyid + usubjid + trt01a + qstestcd ~ avisitn + variable)

panss_nsfs.wide <- na.omit(panss_nsfs.wide)
panss_nsfs.melt2 <- melt(panss_nsfs.wide)
panss_nsfs.wide2 <- dcast(panss_nsfs.melt2, formula = studyid + usubjid + trt01a + variable ~ qstestcd)
panss_nsfs.wide2 <- na.omit(panss_nsfs.wide2)

X <- abind(split(panss_nsfs.wide2, panss_nsfs.wide2$variable), along = 3)
Y <- apply(X[, -c(1:4), ], c(2,3), as.numeric) + 1


################################################################################
# PANSS negative symptom factor score
# All Studies, both indications, data in long format
################################################################################
panss_nsfs <- droplevels((subset(panssitem_all, select = c(studyid, usubjid, avisitn, avisit, qstestcd, qsstresn, trt01a), qstestcd %in% c('QS027N01', 'QS027N02', 'QS027N03', 'QS027N04', 'QS027N06', 'QS027G07','QS027G16') & avisitn != 10.5 & trt01a %in% c('Bitopertin 10mg', 'Bitopertin 20mg', 'Bitopertin 5mg', 'Placebo'))))
panss_nsfs <- panss_nsfs[!duplicated(panss_nsfs),]
panss_nsfs.melt <- melt(panss_nsfs, id.vars = c('studyid', 'usubjid', 'avisitn', 'qstestcd', 'trt01a'), measure.vars = 'qsstresn')

panss_nsfs.melt <- panss_nsfs.melt[!duplicated(subset(panss_nsfs.melt, select = c('studyid', 'usubjid', 'avisitn', 'qstestcd', 'trt01a', 'variable'))), ]


panss_nsfs.melt <- na.omit(panss_nsfs.melt)

ii <- as.numeric(panss_nsfs.melt$usubjid)
jj <- as.numeric(panss_nsfs.melt$qstestcd)
tt <- as.numeric(as.factor(panss_nsfs.melt$avisitn))


n.s <- length(unique(ii))
n.item <- length(unique(jj))
n.times <- length(unique(tt))


Y <- panss_nsfs.melt$value + 1
n.item_cat <- max(unique(Y))
#trt <-  as.numeric(panss_nsfs.melt[panss_nsfs.melt$]$trt01a)
trt <-  as.numeric(panss_nsfs.melt[panss_nsfs.melt$avisitn == 2 & panss_nsfs.melt$qstestcd == 'QS027G07', ]$trt01a)
n_trt <- length(unique(trt))
N_obs <- length(Y)

stan_data <- list(m_mu_theta = 0.0,
                  var_mu_theta = 1.0,
                  N_obs = N_obs,
                  N = n.s,
                  J = n.item,
                  T = n.times,
                  K = n.item_cat,
                  Y = Y,
                  ii = ii,
                  jj = jj,
                  tt = tt,
                  n_trt = n_trt,
                  trt = trt)


##########################################################################################
# AR(1), population-averaged model
# Delta - treatment effect at each time point
##########################################################################################
grad_long_group2 <- '
data {
  int n_trt;                     // number of arms
  real m_mu_theta;               // prior for mean of theta
  real var_mu_theta;             // prior for variance of mean of theta
  int<lower=0> T;                // number of categorical time points
  int<lower=0> N;                // number of subjects
  int<lower=0> J;                // number of items
  int<lower=1, upper=7> K[J];    // max level for each item
  int Y[N, J, T];                // data, 3D array of integers
  int trt[N];                    // treatment trt04a
}

parameters {
    ordered[6] kappa[J];  // intercept : array of ordered vectors
    vector[T] theta[N];                     // ability : theta[N,T] - array of vectors
    real<lower=0.0, upper=10.0> alpha[J];   // discrimination
    real<lower=0.0, upper=1.0> rho;         // Correlation parameter
    vector<lower=-10.0, upper=10.0>[T] mu_theta_raw[n_trt];
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
       kappa[j, k] ~ normal(0.0, 3.0) T[ , 100];
     }
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
        if (Y[i, j, t] <= K[j])
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
###compile the model

# set_cppo(mode = "debug")
# set_cppo(mode = "fast")
c_grad_long_group2 <- stan_model(model_code = 'grad_long_group2')


##########################################################################################
# AR(1), population-averaged model
# Hierarchical priors on item parameters
# Delta - treatment effect at each time point
##########################################################################################

grad_long_group3 <- '
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
    ordered[6] kappa[J];                        // thresholds : array of ordered vectors
    vector<lower=-6, upper = 6>[T] theta[N];    // ability : theta[N,T] - array of vectors
    real<lower=0.0, upper=10.0> alpha[J];       // discrimination
    real<lower=0.0, upper=1.0> rho;             // Correlation parameter
    vector[T] mu_theta_raw[n_trt];
    real <lower=0.0, upper = 5.0> mu_alpha;     // Prior for item discrim.
    real <lower=-20.0, upper = 20.0> mu_kappa;  // prior for item thresholds.
    real <lower=0.0, upper = 5.0> var_alpha;
    real <lower=0.0, upper = 5.0> var_kappa;
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
   real m_mu_alpha;
   real m_mu_kappa;
   real var_mu_alpha;
   real var_mu_kappa;

// hyperpriors for item parameters
   m_mu_alpha <- 1;
   m_mu_kappa <- 0;
   var_mu_alpha <- 100;
   var_mu_kappa <- 100;


// hyperparameters for item parameters;
   mu_alpha ~ normal(m_mu_alpha, sqrt(var_mu_alpha));
   mu_kappa ~ normal(m_mu_kappa, sqrt(var_mu_kappa));
   var_alpha ~ cauchy(0, 3);    //  hyperprior
   var_kappa ~ cauchy(0, 3);

// prior for discrimination, truncated below 0 and above X
   alpha ~ normal(mu_alpha, sqrt(var_alpha));

// prior for item difficulty
  for (j in 1:J) {
     for (k in 1:(K[j] - 1)) {
   kappa[j, k] ~ normal(mu_kappa, sqrt(var_kappa));
    }
   }

// hyperparameters for ability parameters
   for (t in 1:T)
      for (j in 1:n_trt)
          mu_theta_raw[j, t] ~ normal(m_mu_theta, var_mu_theta);


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

#compile the model 
c_grad_long_group3 <- stan_model(model_code = 'grad_long_group3')

# time
T <- dim(Y)[3]
# n.variables    
J <- dim(Y)[2]
# sample size
N <- dim(Y)[1]
# K - max level for ordinal variable for each item
K <- apply(apply(Y, c(2,3), max), 1, max)

trt <- 4- as.numeric(droplevels(panss_nsfs.wide2[panss_nsfs.wide2$variable == '2_qsstresn', ]$trt01a))

n.trt <- length(table(trt))

m.mu.theta <- 0
var.mu.theta <- 1


# data for stan()
stan_data2 <- list(T = T
                   ,N = N
                   ,J = J 
                   ,Y = Y
                   ,K = K                  
                   ,m_mu_theta = m.mu.theta
                   ,var_mu_theta = var.mu.theta
                   ,trt = as.numeric(trt)
                   ,n_trt = n.trt)


      ##   • adapt_engaged (logical) 
      ##   • adapt_gamma (double, positive, defaults to 0.05) 
      ##   • adapt_delta (double, between 0 and 1, defaults to 0.8) 
      ##   • adapt_kappa (double, positive, defaults to 0.75) 
      ##   • adapt_t0 (double, positive, defaults to 10) 
      ##   • adapt_init_buffer (integer, positive, defaults to 75) 
      ##   • adapt_term_buffer (integer, positive, defaults to 50) 
      ##   • adapt_window (integer, positive, defaults to 25) 

 # run chains sequentially
mod_grad_long_group2 <- sampling(c_grad_long_group2, data = stan_data2, iter = 1000, chains = 1)


mod_grad_long_group2b <- pstan(model_code = grad_long_group2, data = stan_data2, warmup = 500, iter = 2000, chains = 3)

mod_grad_long_group2_res <- extract(mod_grad_long_group2b, inc_warmup = FALSE)


# print(fitl1, pars=c("theta", "alpha", "kappastar"), probs=c(.1,.5,.9))
print(mod_grad_long_group2b, pars=c("alpha"), probs=c(.1,.5,.9))
print(mod_grad_long_group2b, pars=c("kappa"), probs=c(.1,.5,.9))
print(mod_grad_long_group2b, pars=c("deltaout"), probs=c(.1,.5,.9))
print(mod_grad_long_group2b, pars=c("mu_theta"), probs=c(.1,.5,.9))

alpha_est <- apply(mod_grad_long_group2_res$alpha, 2, mean)
kappa_est <- apply(mod_grad_long_group2_res$kappa, c(2, 3) , mean)
delta_est <- apply(mod_grad_long_group2_res$deltaout, c(2, 3) , mean)

#trace plots
rstan::traceplot(mod_grad_long_group2, c("alpha"), ncol=2, inc_warmup=F)
rstan::traceplot(mod_grad_long_group2, c("kappa"), ncol=12, inc_warmup=F)
rstan::traceplot(mod_grad_long_group2, c("deltaout"), ncol=12, inc_warmup=F)

mod_grad_long_group2_deltaout <- mod_grad_long_group2_res$deltaout
# mod_grad_long_group2_ci <- coda::HPDinterval(coda::as.mcmc(as.vector(mod_grad_long_group2_deltaout)))
mod_grad_long_group2_ci <- apply(mod_grad_long_group2_deltaout, c(2,3), function(x) coda::HPDinterval(coda::as.mcmc(as.vector(x))))


### Hierarchical model:
 # run chains sequentially
mod_grad_long_group3 <- sampling(c_grad_long_group3, data = stan_data2, iter = 1000, chains = 1)

mod_grad_long_group3 <- sampling(c_grad_long_group3, data = stan_data2, iter = 1000, chains = 1, control = list(adapt_delta=0.95, max_treedepth = 12)) # stepsize=0.001

summary(do.call(rbind, args = get_sampler_params(mod_grad_long_group3, inc_warmup = TRUE)), digits = 2)
lapply(get_sampler_params(mod_grad_long_group3, inc_warmup = TRUE), summary, digits = 2)

# run chains in parallel
mod_grad_long_group3 <- pstan(model_code = grad_long_group3, data = stan_data2, iter = 1000, chains = 3, control = list(adapt_delta=0.8, max_treedepth = 12)) # stepsize=0.001


# mod_grad_long_group2b <- sampling(c_grad_long_group2, data = stan_data2, iter = 1000, chains = 3)

mod_grad_long_group3_res <- extract(mod_grad_long_group3, inc_warmup = FALSE)


# print(fitl1, pars=c("theta", "alpha", "kappastar"), probs=c(.1,.5,.9))
print(mod_grad_long_group3, pars=c("alpha"), probs=c(.1,.5,.9))
print(mod_grad_long_group3, pars=c("deltaout"), probs=c(.1,.5,.9))
print(mod_grad_long_group3, pars=c("mu_theta"), probs=c(.1,.5,.9))

alpha_est <- apply(mod_grad_long_group3_res$alpha, 2, mean)
kappa_est <- apply(mod_grad_long_group3_res$kappa, c(2, 3) , mean)
delta_est <- apply(mod_grad_long_group3_res$deltaout, c(2, 3) , mean)

#trace plots
rstan::traceplot(mod_grad_long_group3, c("alpha"), ncol=2, inc_warmup=F)
rstan::traceplot(mod_grad_long_group3, c("kappa"), ncol=12, inc_warmup=F)


mod_grad_long_group3_deltaout <- mod_grad_long_group3_res$deltaout
mod_grad_long_group3_ci <- coda::HPDinterval(coda::as.mcmc(as.vector(mod_grad_long_group3_deltaout)))

#############################################  
# Model
#############################################
LGRM_Missing_Data <- '
data {
  real m_mu_theta;               // prior mean of the mean of theta
  real var_mu_theta;             // prior variance of the mean of theta
  int<lower=0> N_obs;            // number of observation = N*J*T
  int<lower=0> N;                // number of subjects
  int<lower=0> J;                // number of items
  int<lower=0> T;                // number of categorical time points
  int<lower=1, upper=7> K;       // max level for each item
  int<lower=1, upper=7> Y[N_obs]; // data
  int ii[N_obs];                 // in order to use long format
  int jj[N_obs];
  int tt[N_obs];
  int n_trt;                     // number of arms (number of groups)  
  int trt[N];                    // treatment arm, seems to be a vector of 1, 2, ..., n_trt
}

parameters {
  ordered[K-1] kappa[J];                  // intercept : array of ordered vectors size J,K-1
  vector<lower=-20.0, upper = 20.0>[T] theta[N];   // ability : theta[N,T] - array of vectors
  real<lower=0.0, upper=15.0> alpha[J];   // discrimination
  real<lower=0, upper=1> rho;         // Correlation parameter
  vector[T] mu_theta_raw[n_trt];
  real<lower=0.0, upper = 10.0> mu_alpha; // Hyperparameter for alpha
  real<lower=-20.0, upper = 20.0> mu_kappa;   // Hyperparameter for kappa
  real<lower=0> var_alpha;                // Hyperparameter for alpha
  real<lower=0> var_kappa;                // Hyperparameter for kappa
  vector<lower=-50.0, upper = 50.0>[N_obs] Z; // Latent response
}

transformed parameters {
  vector[T] mu_theta[n_trt];  // prior mean theta, one per group

  mu_theta[1, 1] <- 0.0;      // fix to 0 for identifiability  

  // set the values of mu_theta before sampling - stan complains otherwise
  for (j in 2:n_trt)
    mu_theta[j, 1] <- mu_theta_raw[j, 1];
  for (t in 2:T)
    for (j in 1:n_trt)
      mu_theta[j, t] <- mu_theta_raw[j, t];  
}

model{
  real sigsq_theta;  // variance multiplier in the AR(1) model, drop comment if fixed to 1
  matrix[T, T] var_theta; // covariance matrix of theta
  matrix[T, T] L_theta;   // Cholesky decomposition of theta

  // hyperparameters for item parameters;
  mu_alpha ~ normal(8, 1); // This is close to 8
  mu_kappa ~ normal(0, 1); // This is close to 0
  var_alpha ~ cauchy(0, 5);   // inv gamma hyperprior, this is close to 3
  var_kappa ~ cauchy(0, 5);   // inv gamma hyperprior, this is large
  
  // prior for discrimination, truncated below 0 and above X
  for (j in 1:J){
    alpha[j] ~ normal(mu_alpha, var_alpha); 
    // If mu_alpha=8 and var_alpha=3, this is quite uninformative
  }
  
  // prior for item difficulty
  for (j in 1:J) {
    for (k in 1:(K-1)) {
      kappa[j,k] ~ normal(mu_kappa, var_kappa); 
      // If mu_kappa=0 and var_kappa=large, this is quite uninformative
    }
  }

  // hyperparameters for ability parameters
  for (t in 1:T)
    for (j in 1:n_trt)
      mu_theta_raw[j, t] ~ normal(m_mu_theta, sqrt(var_mu_theta));

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

  L_theta <- cholesky_decompose(var_theta);
  for (i in 1:N){
    theta[i] ~ multi_normal_cholesky(mu_theta[trt[i]], L_theta);
  // print(theta[i]);
  // print(i);
  // print(trt[i]);
  }

  for (n in 1:N_obs){
    Z[n] ~ logistic(alpha[ jj[n] ] * theta[ ii[n] , tt[n] ], 1.0); //  the latent response...
  }

  for (n in 1:N_obs){
    Y[n] ~ ordered_logistic(Z[n], kappa[jj[n]]);
  }

}

generated quantities {
  vector[T-1] deltaout[n_trt-1];  // Treatment effect
  vector[N_obs] RBLR; // Rao-Blackwellized estimates of the latent residual

  // Code the treatment effect - difference in differences
  for (t in 1:(T-1)) {
    for (l in 1:(n_trt - 1)) {
      deltaout[l, t] <-  mu_theta[l+1, t+1] - mu_theta[1, t+1] - (mu_theta[l+1, 1] - mu_theta[1, 1]);
    }
  }

  for (n in 1:N_obs){
    if (Y[n] == 1){
      RBLR[n] <- -exp(- pow( kappa[jj[n],1] - alpha[jj[n]] * theta[ii[n],tt[n]] , 2 )/2.0 )/(sqrt(2 * pi()) * normal_cdf( kappa[jj[n],1] - alpha[jj[n]] * theta[ii[n],tt[n]] , 0.0 , 1.0));
    }
   else if ( (Y[n] != 1) && (Y[n] != K) ) {
      RBLR[n] <- (exp(- pow( kappa[jj[n], Y[n]-1] - alpha[jj[n]] * theta[ii[n],tt[n]] , 2 )/2.0 ) - exp(- pow( kappa[jj[n], Y[n]] - alpha[jj[n]] * theta[ii[n],tt[n]] , 2 )/2.0 ))/(sqrt(2 * pi()) * (normal_cdf( kappa[jj[n], Y[n]] - alpha[jj[n]] * theta[ii[n],tt[n]] , 0.0 , 1.0) - normal_cdf( kappa[jj[n], Y[n]-1] - alpha[jj[n]] * theta[ii[n],tt[n]] , 0.0 , 1.0)));
    }
   else {
          RBLR[n] <- exp(- pow( kappa[jj[n], K-1] - alpha[jj[n]] * theta[ii[n],tt[n]] , 2 )/2.0 )/(sqrt(2 * pi()) * normal_cdf( alpha[jj[n]] * theta[ii[n],tt[n]] - kappa[jj[n], K-1] , 0.0 , 1.0));
    }
  }

}
'
LGRM_Missing_Data.compiled <- stan_model(model_code = LGRM_Missing_Data, verbose = FALSE)

n_chain <-  3
n_iter <- 1000
n_warmup <- 500

# options(mc.cores = parallel::detectCores())
options(mc.cores = 3)

fit_stan <- sampling(LGRM_Missing_Data.compiled,
                     data = stan_data,
                     chains = n_chain,
                     iter = n_iter,
                     thin = 2,
                     warmup = n_warmup,
                  #   init = fit_stan_init,
                     open_progress = TRUE,
                     verbose = TRUE)

fit_stan_init <- get_inits(fit_stan)

mod_grad_long_group3_res <- extract(fit_stan, inc_warmup = FALSE)
mod_grad_long_group3 <- fit_stan

                                        #
# pairs(fit_stan, pars = c("alpha", "kappa", "lp__"))
pairs(fit_stan, pars = c("alpha", "kappa"))

# print(fitl1, pars=c("theta", "alpha", "kappastar"), probs=c(.1,.5,.9))
print(mod_grad_long_group3, pars=c("alpha"), probs=c(.1,.5,.9))
print(mod_grad_long_group3, pars=c("deltaout"), probs=c(.1,.5,.9))
print(mod_grad_long_group3, pars=c("mu_theta"), probs=c(.1,.5,.9))

alpha_est <- apply(mod_grad_long_group3_res$alpha, 2, mean)
kappa_est <- apply(mod_grad_long_group3_res$kappa, c(2, 3) , mean)
delta_est <- apply(mod_grad_long_group3_res$deltaout, c(2, 3) , mean)
theta_est <- apply(mod_grad_long_group3_res$theta, c(2,3), mean)
rho_est   <- mean(mod_grad_long_group3_res$rho)
mu_theta_raw_est <-  apply(mod_grad_long_group3_res$mu_theta_raw, c(2,3), mean)

#trace plots
rstan::traceplot(mod_grad_long_group3, c("alpha"), ncol=2, inc_warmup=F)
rstan::traceplot(mod_grad_long_group3, c("kappa"), ncol=12, inc_warmup=F)
rstan::traceplot(mod_grad_long_group3, c("deltaout"), ncol=2, inc_warmup=F)

mod_grad_long_group3_deltaout <- mod_grad_long_group3_res$deltaout
mod_grad_long_group3_ci <- coda::HPDinterval(coda::as.mcmc(as.vector(mod_grad_long_group3_deltaout)))


save(fit_stan, file = "all_studies_panss_irt_est_12Oct2015.RData", compress = TRUE)


################################################################################
                                        # Residual plots
################################################################################

RBLR_est <- apply(mod_grad_long_group3_res$RBLR, 2, mean)
 RBLR_est <- RBLR_est[-c(which(is.infinite(RBLR_est)))]

# Fitted probabilities: => Also consider doing it in the model! May be more accurate.
# See Albert and Chib 1994, page 7
FP <- rep(0, N_obs)
for (i in 1:N_obs){
  if (Y[i] == 1) { FP[i] <- plogis(kappa_est[jj[i], 1] - alpha_est[jj[i]] * theta_est[ii[i], tt[i]])
  } else if (Y[i] == n.item_cat) {
    FP[i] <- 1 - plogis(kappa_est[jj[i], n.item_cat-1] - alpha_est[jj[i]] * theta_est[ii[i], tt[i]])
  } else {
    FP[i] <- plogis(kappa_est[jj[i], Y[i]] - alpha_est[jj[i]] * theta_est[ii[i], tt[i]]) - 
             plogis(kappa_est[jj[i], Y[i]-1] - alpha_est[jj[i]] * theta_est[ii[i], tt[i]])
  }
}

library(graphics)
library(plot3D)
library(pracma)
library(rgl)

graphics.off() # This closes all of R's graphics window
lines3D(seq(from = FP[1], to = (FP[1]+0.000001), length = length(density(mod_grad_long_group3_res$RBLR[,1])$x)),
      density(mod_grad_long_group3_res$RBLR[,1])$x,
      density(mod_grad_long_group3_res$RBLR[,1])$y,
      phi = 60, theta = -90, ylab = "Residual", xlab = "Fitted prob.",
      xlim = c(0,1), ylim = c(-8,8), col = "white",
      ticktype = "detailed", zlab = "Density") -> res

for (i in 1:N_obs){
  if ((jj[i] == 7) & (Y[i] == 4) & (tt[i] == 1) ){
    lines(trans3D(FP[i],
          density(mod_grad_long_group3_res$RBLR[,i])$x,
          density(mod_grad_long_group3_res$RBLR[,i])$y, pmat = res), col = rgb(0,0,0,alpha=0.2))
  }
}



