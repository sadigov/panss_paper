################################################################################
#
#   Analysis of PANSS from Bitopertin program
#
#
################################################################################

setwd("u:/My Documents/Statistics_projects/IRT/bitopertin/")
set.seed(1234)

require(mirt)
require(reshape2)
require(abind)
require(psych)

load(file='panssitem_all_21aug2015.RData')

panssitem_all$trt01a <- relevel(panssitem_all$trt01a, ref = "Placebo" )

################################################################################
# PANSS negative symptom factor score
# All three Studies in the NS program
################################################################################
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

trt <- panss_nsfs.wide2[panss_nsfs.wide2$variable == '2_qsstresn', ]$trt01a

data1 <- Y[ , , 1]
# data1 <- Y[ , , 7]
data_trt <- cbind(trt, data1)

itemnames <- colnames(data1) 

################################################################################
# PANSS negative symptom factor score
# All Studies
################################################################################
panss_nsfs <- droplevels(subset(panssitem_all, select = c(studyid, usubjid, avisitn, avisit, qstestcd, qsstresn, trt01a), qstestcd %in% c('QS027N01', 'QS027N02', 'QS027N03', 'QS027N04', 'QS027N06', 'QS027G07','QS027G16') & avisitn != 10.5 & trt01a %in% c('Bitopertin 10mg', 'Bitopertin 20mg', 'Bitopertin 5mg', 'Placebo')))
panss_nsfs <- panss_nsfs[!duplicated(panss_nsfs),]
panss_nsfs.melt <- melt(panss_nsfs, id.vars = c('studyid', 'usubjid', 'avisitn', 'qstestcd', 'trt01a'), measure.vars = 'qsstresn')

panss_nsfs.melt <- panss_nsfs.melt[!duplicated(subset(panss_nsfs.melt, select = c('studyid', 'usubjid', 'avisitn', 'qstestcd', 'trt01a', 'variable'))), ]
panss_nsfs.wide <- dcast(panss_nsfs.melt, formula = studyid + usubjid + trt01a + qstestcd ~ avisitn + variable)

# panss_nsfs.wide <- na.omit(panss_nsfs.wide)
panss_nsfs.melt2 <- melt(panss_nsfs.wide)
panss_nsfs.wide2 <- dcast(panss_nsfs.melt2, formula = studyid + usubjid + trt01a + variable ~ qstestcd)
panss_nsfs.wide2 <- na.omit(panss_nsfs.wide2)

X <- abind(split(panss_nsfs.wide2, panss_nsfs.wide2$variable), along = 3)
Y <- apply(X[, -c(1:4), ], c(2,3), as.numeric) + 1

trt <- panss_nsfs.wide2[panss_nsfs.wide2$variable == '2_qsstresn', ]$trt01a

data1 <- Y[ , , 1] # pick baseline visit
# data1 <- Y[ , , 7]
data_trt <- cbind(trt, data1)

itemnames <- colnames(data1) 

################################################################################
# polychoric correlations
################################################################################
psych::polychoric(data1 )


################################################################################
## Test for DIF between groups at baseline

# completely independent model
# help(multipleGroup)
(mod <- multipleGroup(data1, model = 1, group = trt, invariance  = c("free_means")))
coef(mod, simplify = TRUE, IRTpars = TRUE)
coef(mod, simplify = TRUE, IRTpars = FALSE)

mirtCluster() ## to run in parallel

#test whether adding slopes and intercepts constraints results in DIF
mod_constr1 <- DIF(mod, which.par = c('a1', 'd1', 'd2', 'd3', 'd4', 'd5', 'd6'), p.adjust = 'fdr', method = 'QMCEM') # , technical = list(NCYCLES = 1000)



(mod2a <- multipleGroup(data1, model = 1, group = trt, invariance  = c(itemnames, "free_means",  'free_var'), technical = list(NCYCLES = 1000)))
coef(mod2a, simplify = TRUE, IRTpars = TRUE)
coef(mod2a, simplify = TRUE, IRTpars = FALSE)

(mod2b <- multipleGroup(data1, model = 1, group = trt, invariance  = c( "free_means",  'free_var'), technical = list(NCYCLES = 1000)))
coef(mod2b, simplify = TRUE, IRTpars = TRUE)
coef(mod2b, simplify = TRUE, IRTpars = FALSE)

mirtCluster() ## to run in parallel
mod_constr2a <- DIF(mod2a, which.par = c('a1', 'd1', 'd2', 'd3', 'd4', 'd5', 'd6'), p.adjust = 'fdr', method = 'QMCEM') # , technical = list(NCYCLES = 1000)
mod_constr2a

mod_constr2b <- DIF(mod2b, which.par = c('a1', 'd1', 'd2', 'd3', 'd4', 'd5', 'd6'), p.adjust = 'fdr', method = 'QMCEM') # , technical = list(NCYCLES = 1000)
mod_constr2b 

#plot curves between groups
plot(mod)
plot(mod, type = 'trace')
plot(mod, type = 'info')
plot(mod, type = 'RE')

itemplot(mod, item = 2)
itemplot(mod, item = 1)
itemplot(mod, 2, type = 'info')
itemplot(mod, 1, type = 'info')
itemplot(mod, 2, type = 'RE')

## model/item fit also applicable here
itemfit(mod)
M2(mod)

#### Using anchor items, letting means and variances be free except one group
#### using items 3,4,6  as anchors
model_anchor <- multipleGroup(data1, model = 1, group = trt, invariance = c(itemnames[c(3,4,6)], 'free_means', 'free_var'))
coef(model_anchor, simplify = TRUE, IRTpars = TRUE)
anchor <- DIF(model_anchor, which.par = c('a1', 'd1', 'd2', 'd3', 'd4', 'd5', 'd6'), items2test = c(1,2,5,7))
anchor


mirtCluster() ## to run in parallel
mod_constr1 <- DIF(mod, which.par = c('a1', 'd1', 'd2', 'd3', 'd4', 'd5', 'd6')) #test whether adding slopes and intercepts constraints results in DIF

mod_constr2 <- DIF(mod2, which.par = c('a1', 'd1', 'd2', 'd3', 'd4', 'd5', 'd6'))
# DIF(mod.freegroup, which.par = c('a1', 'd'), items2test = 3:5, scheme = 'drop') #drop constraints


#step down procedure for highly constrained model
model_sd <- multipleGroup(data1, 1, group = trt, invariance = c(itemnames, 'free_means', 'free_var'))
model_sd <- multipleGroup(data1, 1, group = trt, invariance = c(itemnames, 'free_means'))
coefs_model_sd <- coef(model_sd, simplify = TRUE, IRTpars = TRUE)
means <- sapply(coefs_model_sd, function(x) do.call("rbind", x[c(2)]))
stepdown <- DIF(model_sd, c('a'), scheme = 'drop_sequential')
stepdown

################################################################################
## Test for DIF between groups at baseline

data1 <- Y[ , , 7] # pick last visit

data_trt <- cbind(trt, data1)

itemnames <- colnames(data1) 

(mod2b <- multipleGroup(data2, model = 1, group = trt, invariance  = c( "free_means",  'free_var'), technical = list(NCYCLES = 1000)))
coef(mod2b, simplify = TRUE, IRTpars = TRUE)
coef(mod2b, simplify = TRUE, IRTpars = FALSE)

mod_constr2a <- DIF(mod2a, which.par = c('a1', 'd1', 'd2', 'd3', 'd4', 'd5', 'd6'), p.adjust = 'fdr', method = 'QMCEM') # , technical = list(NCYCLES = 1000)
mod_constr2a

mod_constr2b <- DIF(mod2b, which.par = c('a1', 'd1', 'd2', 'd3', 'd4', 'd5', 'd6'), p.adjust = 'fdr', method = 'QMCEM') # , technical = list(NCYCLES = 1000)
mod_constr2b 

################################################################################
#
#            From DIF function manual in mirt package
#
#
#
################################################################################

#simulate data where group 2 has a smaller slopes and more extreme intercepts
set.seed(12345)
a1 <- a2 <- matrix(abs(rnorm(15,1,.3)), ncol=1)
d1 <- d2 <- matrix(rnorm(15,0,.7),ncol=1)
a2[1:2, ] <- a1[1:2, ]/3
d1[c(1,3), ] <- d2[c(1,3), ]/4

head(data.frame(a.group1 = a1, a.group2 = a2, d.group1 = d1, d.group2 = d2))
itemtype <- rep('dich', nrow(a1))
N <- 1000
dataset1 <- simdata(a1, d1, N, itemtype)
dataset2 <- simdata(a2, d2, N, itemtype, mu = .1, sigma = matrix(1.5))
dat <- rbind(dataset1, dataset2)
group <- c(rep('D1', N), rep('D2', N))

#### no anchors, all items tested for DIF by adding item constrains one item at a time.
# define a parallel cluster (optional) to help speed up internal functions
mirtCluster()

#  Information matrix with crossprod (not controlling for latent group differences)
model <- multipleGroup(dat, 1, group, SE = TRUE)

#test whether adding slopes and intercepts constraints results in DIF. Plot items showing DIF
resulta1d <- DIF(model, c('a1', 'd'), plotdif = TRUE)
resulta1d



