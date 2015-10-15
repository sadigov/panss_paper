require(Hmisc)
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


################################################################################
#                          Import data
################################################################################

bitodata <- sasxport.get("H:\\cdt4715l\\libraries\\csvdata\\all", method = 'csv')
names(bitodata)

# save(bitodata, file = 'bitodata_21aug2015.RData')
# load(file='bitodata.RData')
panssitem <- subset(bitodata, avisit %in% c('Baseline visit', 'Week 4', 'Week 8', 'Week 12', 'Week 16', 'Week 18', 'Week 20', 'Week 24', 'Week 24 LOCF') & qstestcd %in% c('QS027G01', 'QS027G02', 'QS027G03', 'QS027G04', 'QS027G05', 'QS027G06', 'QS027G07', 'QS027G08', 'QS027G09', 'QS027G10', 'QS027G11', 'QS027G12', 'QS027G13', 'QS027G14', 'QS027G15', 'QS027G16', 'QS027N01', 'QS027N02', 'QS027N03' ,'QS027N04', 'QS027N05', 'QS027N06', 'QS027N07', 'QS027P01', 'QS027P02', 'QS027P03', 'QS027P04', 'QS027P05', 'QS027P06', 'QS027P07') & randfl =='Y')

# Use armcd variable for study arm, use avisit variable for time

                    
rm(bitodata)
save(panssitem, file = 'panssitem_21sep2015.RData')
# load(file='panssitem_21aug2015.RData')


                                        #Prepare data
# N0 LOCF data, PANSS NS
panss_ns <- droplevels(na.omit(subset(panssitem, select = c(studyid, usubjid, avisitn, avisit, qstestcd, qsstresn, arm), qstestcd %in% c('QS027N01', 'QS027N02', 'QS027N03' ,'QS027N04', 'QS027N05', 'QS027N06', 'QS027N07') & !(avisitn %in% c(9, 10.5)) & arm %in% c('Bitopertin 10mg', 'Bitopertin 20mg', 'Bitopertin 5mg', 'Placebo'))))
panss_ns <- panss_ns[!duplicated(panss_ns),]
panss_ns_melt <- melt(panss_ns, id.vars = c("studyid", 'usubjid', 'avisitn', 'qstestcd', 'arm'), measure.vars = 'qsstresn')

panss_ns_melt <- panss_ns_melt[!duplicated(subset(panss_ns_melt, select = c('studyid', 'usubjid', 'avisitn', 'qstestcd', 'arm', 'variable'))), ]
panss_ns_wide <- dcast(panss_ns_melt, formula =  studyid + usubjid + arm + qstestcd ~ avisitn + variable)

panss_ns_wide <- na.omit(panss_ns_wide)
panss_ns_melt2 <- melt(panss_ns_wide)
panss_ns_wide2 <- dcast(panss_ns_melt2, formula = studyid + usubjid + arm + variable ~ qstestcd)
panss_ns_wide2 <- na.omit(panss_ns_wide2)

X <- abind(split(panss_ns_wide2, panss_ns_wide2$variable), along = 3)
Y <- apply(X[, -c(1, 2, 3, 4), ], c(2, 3), as.numeric) + 1

################################################################################
#
#     IRT MODELS
#
################################################################################


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
       kappa[j, k] ~ normal(0.0, 2.5);
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


# time
T <- dim(Y)[3]
# n.variables    
J <- dim(Y)[2]
# sample size
N <- dim(Y)[1]
# K - max level for ordinal variable for each item
K <- apply(apply(Y, c(2,3), max), 1, max)

trt <-  as.numeric(droplevels(panss.wide2[panss.wide2$variable == '2_qsstresn', ]$arm))

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

mod_grad_long_group2 <- sampling(c_grad_long_group2, data = stan_data2, iter = 1000, chains = 3)



mod_grad_long_group2b <- sampling(c_grad_long_group2, data = stan_data2, iter = 10000, chains = 1)

mod_grad_long_group2_res <- extract(mod_grad_long_group2b, inc_warmup = FALSE)


# print(fitl1, pars=c("theta", "alpha", "kappastar"), probs=c(.1,.5,.9))
print(mod_grad_long_group2b, pars=c("alpha"), probs=c(.1,.5,.9))
print(mod_grad_long_group2b, pars=c("deltaout"), probs=c(.1,.5,.9))

alpha_est <- apply(mod_grad_long_group2_res$alpha, 2, mean)
kappa_est <- apply(mod_grad_long_group2_res$kappa, c(2, 3) , mean)
delta_est <- apply(mod_grad_long_group2_res$deltaout, c(2, 3) , mean)

#trace plots
rstan::traceplot(mod_grad_long_group2b, c("alpha"), ncol=2, inc_warmup=F)
rstan::traceplot(mod_grad_long_group2b, c("kappa"), ncol=12, inc_warmup=F)

xyplot(alpha_est ~ a1)
xyplot(kappa_est[, 1] ~ b.mat[1, ])
xyplot(kappa_est[1, ] ~ b.mat[, 1])


mod_grad_long_group3_deltaout <- mod_grad_long_group3_res$deltaout
mod_grad_long_group3_ci <- coda::HPDinterval(coda::as.mcmc(as.vector(mod_grad_long_group3_deltaout)))



################################################################################
#
#                             SEM MODELS
#
################################################################################
require(lavaan)
require(psych)
require(GPArotation)
require(fastICA)
require(nFactors)
require(mirt)

names(panssitem)

                                        #Prepare data
# N0 LOCF data, all PANSS
panss_all <- droplevels(na.omit(subset(panssitem, select = c(studyid, usubjid, avisitn, avisit, qstestcd, qsstresn, arm), qstestcd %in% c('QS027N01', 'QS027N02', 'QS027N03' ,'QS027N04', 'QS027N05', 'QS027N06', 'QS027N07', 'QS027G01', 'QS027G02', 'QS027G03', 'QS027G04', 'QS027G05', 'QS027G06', 'QS027G07', 'QS027G08', 'QS027G09', 'QS027G10', 'QS027G11', 'QS027G12', 'QS027G13', 'QS027G14', 'QS027G15', 'QS027G16', 'QS027P01', 'QS027P02', 'QS027P03', 'QS027P04', 'QS027P05', 'QS027P06', 'QS027P07') & !(avisitn %in% c(9, 10.5)) & arm %in% c('Bitopertin 10mg', 'Bitopertin 20mg', 'Bitopertin 5mg', 'Placebo'))))
panss_all <- panss_all[!duplicated(panss_all),]
panss_all_melt <- melt(panss_all, id.vars = c("studyid", 'usubjid', 'avisitn', 'qstestcd', 'arm'), measure.vars = 'qsstresn')

panss_all_melt <- panss_all_melt[!duplicated(subset(panss_all_melt, select = c('studyid', 'usubjid', 'avisitn', 'qstestcd', 'arm', 'variable'))), ]
panss_all_wide <- dcast(panss_all_melt, formula =  studyid + usubjid + arm + qstestcd ~ avisitn + variable)

panss_all_wide <- na.omit(panss_all_wide)
panss_all_melt2 <- melt(panss_all_wide)
panss_all_wide2 <- dcast(panss_all_melt2, formula = studyid + usubjid + arm + variable ~ qstestcd)
panss_all_wide2 <- na.omit(panss_all_wide2)

X <- abind(split(panss_all_wide2, panss_all_wide2$variable), along = 3)
Y <- apply(X[, -c(1, 2, 3, 4), ], c(2, 3), as.numeric) + 1


                                        # PCA
pca_panss <- principal(Y[, , 1], nfactors = 7, rotate= 'oblimin')
pca_panss <- principal(Y[, , 1], nfactors = 5, rotate= 'oblimin') 
summary(pca_panss)

pca_panss <- princomp(Y[, , 1], cor = TRUE)
loadings(pca_panss) # pc loadings 
plot(pca_panss, type = "lines") # scree plot 
pca_panss$scores # the principal components
biplot(pca_panss)


                                        # Determine Number of Factors to Extract
ev <- eigen(cor(Y[, , 1])) # get eigenvalues
ap <- parallel(subject=nrow(Y[, , 1]),var=ncol(Y[, , 1]), rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS)

                                        # ICA
ica_panss <- fastICA(Y[, , 1], n.comp = 5, alg.typ = "parallel", fun = "logcosh", alpha = 1, method = "C", row.norm = FALSE, maxit = 200, tol = 0.0001, verbose = TRUE)

par(mfrow = c(1, 3))
plot(ica_panss$X, main = "Pre-processed data")
plot(ica_panss$S, main = "ICA components")
plot(ica_panss$X %*% ica_panss$K, main = "PCA components")

summary(ica_panss)

                                        # EFA
# Single factor models

sem_mod1 <- '
N =~ g1*QS027G01 + g2*QS027G02 + g3*QS027G03 + g4*QS027G04 + g5*QS027G05 + g6*QS027G06 + g7*QS027G07 + g8*QS027G08 + g9*QS027G09 + g10*QS027G10 + g11*QS027G11 + g12*QS027G12 + g13*QS027G13 + g14*QS027G14 + g15*QS027G15 + g16*QS027G16 + n1*QS027N01 + n2*QS027N02 + n3*QS027N03 + n4*QS027N04 + n5*QS027N05 + n6*QS027N06 + n7*QS027N07 + p1*QS027P01 + p2*QS027P02 + p3*QS027P03 + p4*QS027P04 + p5*QS027P05 + p6*QS027P06 + p7*QS027P07
'

efa_panss <- sem(sem_mod1, data = Y[ , , 1])

summary(efa_panss, standardized=TRUE)


efa_panss_ord <- sem(sem_mod1, data = Y[ , , 1], ordered = c('QS027N01', 'QS027N02', 'QS027N03' ,'QS027N04', 'QS027N05', 'QS027N06', 'QS027N07', 'QS027G01', 'QS027G02', 'QS027G03', 'QS027G04', 'QS027G05', 'QS027G06', 'QS027G07', 'QS027G08', 'QS027G09', 'QS027G10', 'QS027G11', 'QS027G12', 'QS027G13', 'QS027G14', 'QS027G15', 'QS027G16', 'QS027P01', 'QS027P02', 'QS027P03', 'QS027P04', 'QS027P05', 'QS027P06', 'QS027P07'))

summary(efa_panss_ord, standardized=TRUE)


### Unrestricted/Exploratory MIRT, 5 factor model
# When specifying a single number greater than 1 as the model input to mirt an exploratory IRT
# model will be estimated. Rotation and target matrix options are available if they are passed to
# generic functions such as summary-method and fscores. Factor means and variances are fixed to
# ensure proper identification.
# http://faculty.psy.ohio-state.edu/edwards/documents/IMPS6.30.08.pdf

emirt_panss <- mirt(data = Y[ , , 1], model = 5, itemtype = 'graded')
summary(emirt_panss, rotate = 'varimax' )

summary(emirt_panss, rotate = 'oblimin' )

summary(emirt_panss, rotate = 'simplimax' )

################################################################################
###                   MIRT
################################################################################

                                        # At baseline
mirt_mod1 <- '
N = QS027N01, QS027N02, QS027N03, QS027N04, QS027N06, QS027G07, QS027G16 
P = QS027P01, QS027P03, QS027P05, QS027P06, QS027G09, QS027G12
A = QS027G01, QS027G02, QS027G03, QS027G04, QS027G06
D = QS027G05, QS027G11, QS027G13, QS027G15, QS027G10, QS027N05, QS027N07, QS027P02  
E = QS027G08, QS027G14, QS027P04, QS027P07
COV =  N*P, N*A, N*D, N*E, P*A, P*D, P*E, A*D, A*E, D*E 
'

mirt1 <- mirt.model(mirt_mod1)

mirt1_est <- mirt(data = Y[, , 1], model = mirt1)
summary(mirt1_est)

save.image(file = '21sept2015.RData')


