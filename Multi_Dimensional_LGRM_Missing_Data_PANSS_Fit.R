#############################################
# Multi-dimensional, longitudinal, graded,  #
# multiple-group, IRT, with missing data    #
#                                           #
#            A PANSS DATA FIT               #
#                                           #
# Authors:      Alexandre Moesching         #
#               Shamil Sadikhov             #
#                                           #
# Version:      V1 (28.10.2015)             #
#############################################
graphics.off() # This closes all of R's graphics windows.
rm(list=ls()) # Careful! This clears all of R's memory!
setwd("~/R/MultiDimensional-PANSS-Fit_V1")
set.seed(1234)
library(rstan) 
library(abind)
library(reshape2)
library(lattice)
library(mvtnorm)
library(rstanmulticore)
library(coda)
library(Hmisc)
library(numbers)
library(graphics)
library(plot3D)
library(pracma)
# library(rgl)
# install.packages("devtools")
# require(devtools)
# install_github('nathanvan/rstanmulticore')

#############################################  
# The model
#############################################
MultiDim_LGRM_Missing_Data <- '
functions {
  // Kronecker product
  matrix kronecker(matrix X, matrix Y){
    int x1; int x2; int y1; int y2;
    matrix[rows(X)*rows(Y), cols(X)*cols(Y)] Z;
    
    x1 <- rows(X); x2 <- cols(X);
    y1 <- rows(Y); y2 <- cols(Y);
    
    for (i in 1:x1){
      for (j in 1:x2){
        for (k in 1:y1){
          for (l in 1:y2){
            Z[k+y1*(i-1), l+y2*(j-1)] <- X[i,j]*Y[k,l];
          }
        }
      }
    }
    return Z;
  }
}

data {
  real m_mu_theta;                  // Prior mean of the mean of theta
  real<lower=0> var_mu_theta;       // Prior variance of the mean of theta
  
  int<lower=1> N_obs;               // Number of observation = N*J*T - #NA
  int<lower=1> N;                   // Number of subjects
  int<lower=1> J;                   // Number of items
  int<lower=1> T;                   // Number of categorical time points
  int<lower=1> A;                   // Number of different abilities
  
  int<lower=1, upper=7> K;          // Max level for each item
  int<lower=1, upper=K> Y[N_obs];   // Data in the long format
  int<lower=1> ii[N_obs];           // Response ID
  int<lower=1> jj[N_obs];           // Response item
  int<lower=1> tt[N_obs];           // Response time
  int<lower=1> FL[J];               // Vector giving the factor loadings of each item
  int<lower=1> n_trt;               // Number of arms/groups
  int<lower=1, upper=n_trt> trt[N]; // Treatment arm, vector of 1, 2, ..., n_trt
}

parameters {
  vector<lower=-20, upper=20>[T*A] theta_long[N];       // Ability in long format: theta_long[N,T*A]. theta_long[i] should give a vector of abilities across time
  vector<lower=-20, upper=20>[T*A] mu_theta_id[n_trt];  // Ability mean : mu_theta_id[n_trt, T*A]
  
  vector<lower=0, upper=15>[A] alpha_id[J];             // Discrimination matrix
  real<lower=0, upper=15> mu_alpha;                     // Hyperparameter for alpha
  real<lower=0> var_alpha;                              // Hyperparameter for alpha
  
  ordered[K-1] kappa[J];                                // Thresholds : array of ordered vectors size J,K-1
  real<lower=-20, upper=20> mu_kappa;                   // Hyperparameter for kappa
  real<lower=0> var_kappa;                              // Hyperparameter for kappa
  
  real<lower=-1, upper=1> rho;                          // Correlation parameter in the AR(1) matrix
  corr_matrix[A] Sigma_theta;                           // Covariance between the abilities  
}

transformed parameters {
  vector<lower=-20, upper=20>[T*A] mu_theta[n_trt];  // Prior mean for theta, one per group: mu_theta[n_trt, T*A]
  vector<lower=-20, upper=20>[A] theta[N, T];        // Ability in array format: theta[N,T,A]. theta[i,t] should give a vector of abilities
  vector<lower=0, upper=15>[A] alpha[J];             // Discrimination matrix, only one non-zero value per row
  
  matrix[A*T, A*T] L_theta;                          // Cholesky decomposition of kronecker(Sigma_theta, Omega_theta)
  matrix[T, T] Omega_theta;                          // Covariance across time, AR(1) process
  
  // Identifiability constraints on mu_theta
  for (t in 1:(T*A)){
    for (j in 1:n_trt){
      mu_theta[j, t] <- mu_theta_id[j, t]; 
    }
  }
  for (i in 1:A){
    mu_theta[1, 1+(i-1)*T] <- 0.0; // Fix to zero for each ability at time point 1 group 1
  }
  
  // Identifiability constraints on alpha
  for (j in 1:J){
    for (l in 1:A){
      if (FL[j] == l){
        alpha[j,l] <- alpha_id[j,l]; // Item j loads into ability l.
      } else {
        alpha[j,l] <- 0.0; // Item j does not load into ability l => fix discrimination to 0
      }
    }
  }
  
  // Format change for theta
  for (i in 1:N){
    for (t in 1:T){
      for (s in 1:A){
        theta[i,t,s] <- theta_long[i, t+((s-1)*T)];
      }
    }
  }
  
  // AR(1) structure for the time-covariance of theta
  Omega_theta[1,1] <- 1.0;
  for (t in 2:T) {
    Omega_theta[t,t] <- 1.0;
    for (j in 1:(t-1)) {
      Omega_theta[t,j] <- pow(rho, t-j);
      Omega_theta[j,t] <- pow(rho, t-j);
    }
  }
  
  // Cholesky decompoistion of the full covariance matrix
  L_theta <- cholesky_decompose(kronecker(Sigma_theta, Omega_theta));
}

model{
  // Hyperparameters for items;
  mu_alpha ~ normal(1, 5);
  mu_kappa ~ normal(0, 5);
  var_alpha ~ normal(1, 1);
  var_kappa ~ normal(1, 1);
  
  // Prior for item discrimination
  for (j in 1:J){
    for (s in 1:A){
      alpha_id[j,s] ~ normal(mu_alpha, var_alpha);
    }
  }
  
  // Prior for items thresholds
  for (j in 1:J) {
    for (k in 1:(K-1)) {
      kappa[j,k] ~ normal(mu_kappa, var_kappa);
    }
  }
  
  // Hyperparameters for ability
  for (t in 1:(T*A)) {
    for (j in 1:n_trt) {
      mu_theta_id[j, t] ~ normal(m_mu_theta, sqrt(var_mu_theta));
    }
  }
  
  // Prior for the covariance between abilities
  Sigma_theta ~ lkj_corr(1);
  
  // Prior for ability
  for (i in 1:N){
    theta_long[i] ~ multi_normal_cholesky(mu_theta[trt[i]], L_theta);
  }
  
  // The likelihood
  for (n in 1:N_obs){
    Y[n] ~ ordered_logistic( dot_product( alpha[ jj[n] ], theta[ ii[n] , tt[n] ]), kappa[jj[n]]);
  }

}

generated quantities {
  real deltaout[n_trt-1, T-1, A];  // Treatment effect
  
  // Code the treatment effect - difference in differences
  for (t in 1:(T-1)) {
    for (l in 1:(n_trt - 1)) {
      for (s in 1:A) {
        deltaout[l, t, s] <-  mu_theta[l+1, t+1+(s-1)*T] - mu_theta[1, t+1+(s-1)*T] - (mu_theta[l+1, 1+(s-1)*T] - mu_theta[1, 1+(s-1)*T]);
      }
    }
  }
}
'

#############################################
# Compilation of the model
#############################################
# MultiDim_LGRM_Missing_Data.compiled <- stan_model(model_code = MultiDim_LGRM_Missing_Data, verbose=FALSE)
# save(list = c("MultiDim_LGRM_Missing_Data.compiled"), file = "CompiledModels.RData")
# rm(MultiDim_LGRM_Missing_Data.compiled)
load("CompiledModels.RData")

#############################################
# PANSS-Data management
#############################################
Sys.setlocale("LC_TIME", "C") # with incorrect locale as.Date does not recognize months, months appear in German

# We load the data
load("panssitem_all_21aug2015.RData")

# This shows the number of observation for each of the 6 studies
table(panssitem_all$studyid)

# We put in the correct order the different groups
panssitem_all$trt01a <- relevel(panssitem_all$trt01a, ref = 'Bitopertin 5mg')
panssitem_all$trt01a <- relevel(panssitem_all$trt01a, ref = 'Placebo')

# This shows that 1st level is Placebo, 2nd level is 5mg, ..., 4th level is 20mg
levels(panssitem_all$trt01a)

# We take a subset of panssitem_all
panss_subset <- subset(# The subset is made on the dataset:
  panssitem_all,
  # Selection rules:
  # - the items we take:
  qstestcd %in% c('QS027G01', 
                  'QS027G02', 
                  'QS027G03',
                  'QS027G04', 
                  'QS027G05', 
                  'QS027G06',
                  'QS027G07', 
                  'QS027G08',
                  'QS027G09', 
                  'QS027G10', 
                  'QS027G11',
                  'QS027G12', 
                  'QS027G13',
                  'QS027G14', 
                  'QS027G15',
                  'QS027G16',
                  'QS027N01', 
                  'QS027N02', 
                  'QS027N03',
                  'QS027N04', 
                  'QS027N05', 
                  'QS027N06', 
                  'QS027N07',
                  'QS027P01', 
                  'QS027P02', 
                  'QS027P03',
                  'QS027P04', 
                  'QS027P05', 
                  'QS027P06', 
                  'QS027P07')
  # - the visit we don't take:
  & avisitn != 10.5
  # - the arms we select:
  & trt01a %in% c('Placebo', 'Bitopertin 5mg', 'Bitopertin 10mg', 'Bitopertin 20mg'),
  # The variables we keep:
  select = c(studyid, usubjid, avisitn, avisit, qstestcd, qsstresn, trt01a))

# If we look at: levels(panss_subset$qstestcd) we see that some items are still there, so we drop them!
panss_droplevels <- droplevels(panss_subset) # Now look at levels(panss_droplevels$qstestcd) and see the difference

# There might be lines that are the same, so we keep those which only appear once
panss_noduplicate <- panss_droplevels[!duplicated(panss_droplevels),]

# Data are melted with measured variable qsstresn that gives the answer from 0 to 6
panss_melt <- melt(panss_noduplicate,
                   id.vars = c('studyid', 'usubjid', 'avisitn', 'qstestcd', 'trt01a'),
                   measure.vars = 'qsstresn')

# Now forget for a moment the observed values, look if there are still duplicated values and drop them
panss_subset2 <- subset(panss_melt, 
                        select = c('studyid', 'usubjid', 'avisitn', 'qstestcd', 'trt01a', 'variable'))
panss_melt2 <- panss_melt[!duplicated(panss_subset2),]

# Now for each respondent and for each item, we look at the observation across the time
panss_wide <- dcast(panss_melt2, 
                    formula = studyid + usubjid + trt01a + qstestcd ~ avisitn + variable)

# At that time, we used to erase rows with NA values. For this we used:
#panss_wide <- na.omit(panss_wide)

# We put them again in the long format (what we gain is that if a respondent did not respond at a certain time point, we get an NA)
panss_melt3 <- melt(panss_wide, 
                    id.vars = c('studyid', 'usubjid', 'trt01a', 'qstestcd'))

# We put them again in a table format, where for each respondent and each time point, we look at the observation across the items
panss_wide2 <- dcast(panss_melt3, 
                     formula = studyid + usubjid + trt01a + variable ~ qstestcd)

# At that time, we used to erase rows with NA values. For this we used:
#panss_wide2 <- na.omit(panss_wide2)

# We put the data into a 3D array
X <- abind(split(panss_wide2, panss_wide2$variable), along = 3)

# The first 4 columns are variable names, we drop them. 
# Then the observation are transformed into numeric values
# Finally we add one so obtain a scale from 1 to 7
Y <- apply(X[, -c(1:4), ], c(2,3), as.numeric) + 1

# Now we take the item names, attention to the order
ItemNames <- dimnames(Y)[[2]]

# Properties of the dataset
n.s        <- dim(Y)[1]
n.item     <- dim(Y)[2]
n.times    <- dim(Y)[3]
n.item_cat <- 7 # I state here that the scale is the same for each item, not a vector


trt        <- as.numeric(droplevels(panss_wide2[panss_wide2$variable == '2_qsstresn', ]$trt01a))
trt_levels <- levels(panss_wide2[panss_wide2$variable == '2_qsstresn', ]$trt01a)  ## That is:
n_trt      <- length(table(trt))                                                  ## trt = 1 <-> Placebo
## trt = 2 <-> 05mg
## trt = 3 <-> 10mg
## trt = 4 <-> 20mg

# Factor loadings: N = Negative, P = Positive, A = Anxiety, D = Disorganised, E = Excited
Factors <- list(N = c("QS027N01", "QS027N02", "QS027N03", "QS027N04", "QS027N06", "QS027G07", "QS027G16"),
                P = c("QS027P01", "QS027P03", "QS027P05", "QS027P06", "QS027G09", "QS027G12"),
                A = c("QS027G01", "QS027G02", "QS027G03", "QS027G04", "QS027G06"),
                D = c("QS027G05", "QS027G11", "QS027G13", "QS027G15", "QS027G10", "QS027N05", "QS027N07", "QS027P02"),
                E = c("QS027G08", "QS027G14", "QS027P04", "QS027P07"))

factor_loadings <- rep(0, length(ItemNames))
for (i in 1:length(ItemNames)){
  for (j in 1:length(Factors)){
    if (ItemNames[i] %in% Factors[[j]]){
      factor_loadings[i] <- j
    }
  }
}
dim <- length(Factors)

# Put the data in a long format
Y <- c(Y)
N_obs <- length(Y)

# Respondent identifier
ii <- rep(c(1:n.s), n.times*n.item)

# Item identifier
jj <- rep(0, n.s*n.item)
for (i in 1:n.item) {
  jj[(n.s*(i-1)+1):(n.s*i)] <- rep(i, n.s)
}
jj <- rep(jj, n.times)

# Time identifier
tt <- rep(0, n.s*n.item*n.times)
for (i in 1:n.times) {
  tt[(n.s*n.item*(i-1)+1):(n.s*n.item*i)] <- rep(i, n.s*n.item)
}

# We drop the NA values
list_indices_with_NA <- c(1:N_obs)
indices_to_drop <- which(is.na(Y))
list_indices_without_NA <- list_indices_with_NA[-indices_to_drop]

Y <- na.omit(Y[list_indices_without_NA])
ii <- na.omit(ii[list_indices_without_NA])
jj <- na.omit(jj[list_indices_without_NA])
tt <- na.omit(tt[list_indices_without_NA])
N_obs <- length(list_indices_without_NA)

#############################################
# Parameters of the model
#############################################
stan_data <- list(m_mu_theta = 0.0,
                  var_mu_theta = 1.0,
                  N_obs = N_obs,
                  N = n.s,
                  J = n.item,
                  T = n.times,
                  A = dim,
                  K = n.item_cat,
                  Y = Y,
                  ii = ii,
                  jj = jj,
                  tt = tt,
                  FL = factor_loadings,
                  n_trt = n_trt,
                  trt = trt)


# init_stan <- function(){
#   a_mod <- avec+0.000001
#   kappa <- bmat
#   theta_long <- cbind(theta[,,1], theta[,,2])
#   alpha_id <- a_mod
#   rho <- rho
#   mu_theta_id <- matrix(rnorm(n_trt*n.times*dim),nrow=n_trt)
#   mu_alpha <- 1
#   mu_kappa <- 0
#   var_alpha <- 10
#   var_kappa <- 10
#   Sigma_theta <- sigma
#   return(list(kappa = kappa,
#               theta_long = theta_long,
#               alpha_id = alpha_id,
#               rho = rho,
#               mu_theta_id = mu_theta_id, 
#               mu_alpha = mu_alpha,
#               mu_kappa = mu_kappa,
#               var_alpha = var_alpha,
#               var_kappa = var_kappa,
#               Sigma_theta = Sigma_theta))
# }

rm(list = c("i", "j", "panss_droplevels", "panss_melt", "panss_melt2", "panss_melt3", 
            "panss_noduplicate", "panss_subset", "panss_subset2", "panss_wide", "panss_wide2"))

#############################################
# Sampling
#############################################
n_chain <- 1
n_iter <- 1000
n_warmup <- 500

# init_stan <- lapply(1:n_chain, function(i) init_stan()) # We need to initialize values for each chains!
fit_stan <- sampling(MultiDim_LGRM_Missing_Data.compiled,
                     data = stan_data,
                     chains = n_chain,
                     iter = n_iter,
                     warmup = n_warmup,
                     #init = init_stan,
                     verbose = FALSE)
save(fit_stan, file = "Fit_StanPANSS.RData")
q(save = "no")

load("Fit_Stan.RData") # Factor loadings

#############################################
# Result analysis
#############################################
fit_stan.res <- rstan::extract(fit_stan)
kappa_est <- apply(fit_stan.res$kappa, c(2,3), mean)
theta_est <- apply(fit_stan.res$theta, c(2,3,4), mean)
alpha_est <- apply(fit_stan.res$alpha, c(2,3), mean)
rho_est <- mean(fit_stan.res$rho)
mu_theta_id_est <- apply(fit_stan.res$mu_theta_id, c(2,3), mean)
mu_alpha_est <- mean(fit_stan.res$mu_alpha)
mu_kappa_est <- mean(fit_stan.res$mu_kappa)
var_alpha_est <- mean(fit_stan.res$var_alpha)
var_kappa_est <- mean(fit_stan.res$var_kappa)
delta_est <- apply(fit_stan.res$deltaout, c(2, 3, 4) , mean)
sigma_theta_est <- apply(fit_stan.res$Sigma_theta, c(2, 3) , mean)
