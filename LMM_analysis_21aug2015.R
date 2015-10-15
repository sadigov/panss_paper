require(Hmisc)
require(rstan)
require(rstanmulticore)
require(reshape2)
require(abind)
require(ggplot2)
require(plyr)

# Sys.setlocale("LC_TIME", "English") # with incorrect locale as.Date does not recognize months, months appear in German
setwd("u:/My Documents/Statistics_projects/IRT/bitopertin/")


################################################################################
#                          Import data
################################################################################

# bitodata <- sasxport.get("u:/My Documents/Statistics_projects/IRT/bitopertin/csvdata", method = 'csv')
# names(bitodata)

# save(bitodata, file = 'bitodata_21aug2015.RData')
# load(file='bitodata_21aug2015.RData')
panssitem2 <- subset(bitodata, avisit %in% c('Baseline visit', 'Week 4', 'Week 8', 'Week 12', 'Week 16', 'Week 18', 'Week 20', 'Week 24', 'Week 24 LOCF') & qstestcd %in% c('PANSFNEG') & randfl =='Y')

# Use armcd variable for study arm, use avisit variable for time

                    
# rm(bitodata)
# save(panssitem2, file = 'panssitem2_06Oct2015.RData')
# load(file='panssitem2_06Oct2015.RData')

load(file='panssitem_all_21aug2015.RData')


panssitem_all2 <- subset(bitodata_all, avisit %in% c('Baseline visit', 'Week 4', 'Week 8', 'Week 12', 'Week 16', 'Week 18', 'Week 20', 'Week 24', 'Week 24 LOCF') & qstestcd %in% c('PANSFNEG') & randfl =='Y')
panssitem_all$trt01a <- relevel(panssitem_all$trt01a, ref = "Placebo" )
################################################################################
#                               Prepare data
################################################################################

panssitem2$trt01a <- relevel(panssitem2$trt01a, ref = "Placebo")

# N0 LOCF data

# PANSS negative symptom factor score
# Study 310
panss_nsfs_tot <- droplevels(na.omit(subset(panssitem2, select = c(studyid, usubjid, avisitn, avisit, qstestcd, qsstresn, trt01a), qstestcd %in% c('PANSFNEG') & avisitn != 10.5 & trt01a %in% c('Bitopertin 10mg', 'Bitopertin 20mg', 'Bitopertin 5mg', 'Placebo'))))
panss_nsfs_tot <- panss_nsfs_tot[!duplicated(panss_nsfs_tot), ]
panss_nsfs_tot.melt <- melt(panss_nsfs_tot, id.vars = c('studyid', 'usubjid', 'avisitn', 'qstestcd', 'trt01a'), measure.vars = 'qsstresn')

panss_nsfs_tot.melt <- panss_nsfs_tot.melt[!duplicated(subset(panss_nsfs_tot.melt, select = c('studyid', 'usubjid', 'avisitn', 'qstestcd', 'trt01a', 'variable'))), ]

panss_nsfs_tot.wide <- dcast(panss_nsfs_tot.melt, formula = studyid + usubjid + trt01a + qstestcd ~ avisitn + variable)

panss_nsfs_tot.wide <- na.omit(panss_nsfs_tot.wide)
panss_nsfs_tot.melt2 <- melt(panss_nsfs_tot.wide)
panss_nsfs_tot.wide2 <- dcast(panss_nsfs_tot.melt2, formula = studyid + usubjid + trt01a + variable ~ qstestcd)
panss_nsfs_tot.wide2 <- na.omit(panss_nsfs_tot.wide2)

# N0 LOCF data

# PANSS negative symptom factor score
# All six studies
panss_nsfs_tot <- droplevels(na.omit(subset(panssitem_all2, select = c(studyid, usubjid, avisitn, avisit, qstestcd, qsstresn, trt01a), qstestcd %in% c('PANSFNEG') & avisitn != 10.5 & trt01a %in% c('Bitopertin 10mg', 'Bitopertin 20mg', 'Bitopertin 5mg', 'Placebo'))))
panss_nsfs_tot <- panss_nsfs_tot[!duplicated(panss_nsfs_tot), ]
panss_nsfs_tot.melt <- melt(panss_nsfs_tot, id.vars = c('studyid', 'usubjid', 'avisitn', 'qstestcd', 'trt01a'), measure.vars = 'qsstresn')

panss_nsfs_tot.melt <- panss_nsfs_tot.melt[!duplicated(subset(panss_nsfs_tot.melt, select = c('studyid', 'usubjid', 'avisitn', 'qstestcd', 'trt01a', 'variable'))), ]

panss_nsfs_tot.wide <- dcast(panss_nsfs_tot.melt, formula = studyid + usubjid + trt01a + qstestcd ~ avisitn + variable)

panss_nsfs_tot.wide <- na.omit(panss_nsfs_tot.wide)
panss_nsfs_tot.melt2 <- melt(panss_nsfs_tot.wide)
panss_nsfs_tot.wide2 <- dcast(panss_nsfs_tot.melt2, formula = studyid + usubjid + trt01a + variable ~ qstestcd)
panss_nsfs_tot.wide2 <- na.omit(panss_nsfs_tot.wide2)


################################################################################
                                        #
################################################################################
require(nlme)

stdize <- function(x) (x-mean(x))/sd(x)

panss_nsfs_tot.melt$avisitn <- as.factor(panss_nsfs_tot.melt$avisitn)
panss_nsfs_tot.melt$value_std <- stdize(panss_nsfs_tot.melt$value)

mod1 <- lme(value_std ~ trt01a + avisitn + avisitn:trt01a, random = ~1|usubjid,  correlation = corAR1(), data = panss_nsfs_tot.melt)
summary(mod1)

mod1.fixef <- fixef(mod1)

# confint(mod1, parm = c("trt01aBitopertin 10mg:avisitn4",  "trt01aBitopertin 20mg:avisitn4", "trt01aBitopertin 10mg:avisitn6", "trt01aBitopertin 20mg:avisitn6", "trt01aBitopertin 10mg:avisitn7", "trt01aBitopertin 20mg:avisitn7",  "trt01aBitopertin 10mg:avisitn8", "trt01aBitopertin 20mg:avisitn8",  "trt01aBitopertin 10mg:avisitn9", "trt01aBitopertin 20mg:avisitn9",  "trt01aBitopertin 10mg:avisitn10", "trt01aBitopertin 20mg:avisitn10"))

mod1.ci <- intervals(mod1, level = 0.95, which = "fixed")[[1]]

mod1.ci.res <- as.data.frame(mod1.ci[11:28, ])
mod1.ci.res$avisitn <- rep(c(4,6,7,8,9,10), each = 3)
mod1.ci.res$trt <- c(rep(c("10mg", "20mg", "5mg"), 6))

mod1.ci.res <- mod1.ci.res[order(mod1.ci.res$trt, mod1.ci.res$avisitn), ]

# panel plot
# xyplot(value ~ ES | model*method, data = two_pl_results_melt[c(1:16), ], type='b', ylab = "Power")

                                        # Overlaid
pdf(file = "panss_nsfs_all_mmrm.pdf")
xYplot(Cbind(est., lower, upper) ~ avisitn, type='b', ylab = "PANSS NSFS", group = trt, data = mod1.ci.res, main = "Std. total Score, MMRM", ylim = c(-0.18, 0.18))
  ##      panel = function(x, y, ...) { 
  ##        panel.xYplot(x, y, ...)
  ##        panel.xyplot(mod1.ci.res$avisitn, mod1.ci.res$value,  jitter.data = TRUE, cex=0.5, horizontal=FALSE,  ...)
  ## })
dev.off()

################################################################################
                                        # Analysis of the MLIRT results

load('panss_est07oct2015.RData')

mod_grad_long_group2_res <- extract(mod_grad_long_group2c, inc_warmup = FALSE)


# print(fitl1, pars=c("theta", "alpha", "kappastar"), probs=c(.1,.5,.9))
print(mod_grad_long_group2c, pars=c("alpha"), probs=c(.1,.5,.9))
print(mod_grad_long_group2c, pars=c("kappa"), probs=c(.1,.5,.9))
print(mod_grad_long_group2c, pars=c("deltaout"), probs=c(.1,.5,.9))
print(mod_grad_long_group2c, pars=c("mu_theta"), probs=c(.1,.5,.9))

alpha_est <- apply(mod_grad_long_group2_res$alpha, 2, mean)
kappa_est <- apply(mod_grad_long_group2_res$kappa, c(2, 3) , mean)
delta_est <- apply(mod_grad_long_group2_res$deltaout, c(2, 3) , mean)

#trace plots
rstan::traceplot(mod_grad_long_group2c, c("alpha"), ncol=2, inc_warmup=F)
rstan::traceplot(mod_grad_long_group2c, c("kappa"), ncol=12, inc_warmup=F)
rstan::traceplot(mod_grad_long_group2c, c("deltaout"), ncol=12, inc_warmup=F)

mod_grad_long_group3_deltaout <- mod_grad_long_group3_res$deltaout
# mod_grad_long_group2_ci <- coda::HPDinterval(coda::as.mcmc(as.vector(mod_grad_long_group2_deltaout)))
mod_grad_long_group3_ci <- apply(mod_grad_long_group3_deltaout, c(2,3), function(x) coda::HPDinterval(coda::as.mcmc(as.vector(x))))

deltaout <- adply(mod_grad_long_group3_ci, 1:3)
deltaout$trt <- rep(rep(c("10mg", "20mg", "5mg"), 6), each = 2)
deltaout$avisitn <- rep(c(4,6,7,8,9,10) , each = 6)
deltaout$ci <- rep(c( "Lower", "Upper"),  18)
deltaout2 <- deltaout[, c('V1', 'trt', 'avisitn', 'ci')]
names(deltaout2)[1] <- 'confint'

deltaout_melt <- melt(deltaout2, id.vars = c('trt', 'avisitn', 'ci'), measure.vars = 'confint')
deltaout_wide <- dcast(deltaout_melt,  formula = trt + avisitn ~ ci, value.var = "value" )

dimnames(delta_est)[[1]] <- c("10mg", "20mg", "5mg")
dimnames(delta_est)[[2]] <- c(4,6,7,8,9,10)

delta_est2 <- as.data.frame(delta_est)
delta_est2_melt <- melt(delta_est2 )
names(delta_est2_melt) <- c('avisitn', 'estimate')
delta_est2_melt$trt <- rep(c("10mg", "20mg", "5mg"), 6)
delta_est2_melt$avisitn <- as.numeric(as.character(delta_est2_melt$avisitn))


deltaout_final <- merge(delta_est2_melt, deltaout_wide, by = c('avisitn', 'trt'))
deltaout_final <- deltaout_final[order(deltaout_final$trt, deltaout_final$avisitn), ]


                                        # Overlaid
pdf(file = "panss_nsfs_all_lirt.pdf")
Hmisc::xYplot(Cbind(estimate, Lower, Upper) ~ avisitn, type='b', ylab = "PANSS NSFS", group = trt, data = deltaout_final, main = "LIRT", ylim = c(-0.18, 0.18))
dev.off()
