
################################################################################
#
#                             SEM MODELS
#
################################################################################

setwd("u:/My Documents/Statistics_projects/IRT/bitopertin/")
set.seed(1234)

require(lavaan)
require(psych)
require(GPArotation)
## require(fastICA)
require(nFactors)
require(mirt)


load(file='panssitem_all_21aug2015.RData')

names(panssitem_all)

                                        #Prepare data
# N0 LOCF data, all PANSS
panss_all <- droplevels(na.omit(subset(panssitem_all, select = c(studyid, usubjid, avisitn, avisit, qstestcd, qsstresn, arm), qstestcd %in% c('QS027N01', 'QS027N02', 'QS027N03' ,'QS027N04', 'QS027N05', 'QS027N06', 'QS027N07', 'QS027G01', 'QS027G02', 'QS027G03', 'QS027G04', 'QS027G05', 'QS027G06', 'QS027G07', 'QS027G08', 'QS027G09', 'QS027G10', 'QS027G11', 'QS027G12', 'QS027G13', 'QS027G14', 'QS027G15', 'QS027G16', 'QS027P01', 'QS027P02', 'QS027P03', 'QS027P04', 'QS027P05', 'QS027P06', 'QS027P07') & !(avisitn %in% c(9, 10.5)) & arm %in% c('Bitopertin 10mg', 'Bitopertin 20mg', 'Bitopertin 5mg', 'Placebo'))))
panss_all <- panss_all[!duplicated(panss_all),]
panss_all_melt <- melt(panss_all, id.vars = c("studyid", 'usubjid', 'avisitn', 'qstestcd', 'arm'), measure.vars = 'qsstresn')

panss_all_melt <- panss_all_melt[!duplicated(subset(panss_all_melt, select = c('studyid', 'usubjid', 'avisitn', 'qstestcd', 'arm', 'variable'))), ]
panss_all_wide <- dcast(panss_all_melt, formula =  studyid + usubjid + arm + qstestcd ~ avisitn + variable)

# panss_all_wide <- na.omit(panss_all_wide)
panss_all_melt2 <- melt(panss_all_wide)
panss_all_wide2 <- dcast(panss_all_melt2, formula = studyid + usubjid + arm + variable ~ qstestcd)
panss_all_wide2 <- na.omit(panss_all_wide2)

X <- abind(split(panss_all_wide2, panss_all_wide2$variable), along = 3)
Y <- apply(X[, -c(1, 2, 3, 4), ], c(2, 3), as.numeric) + 1


                                        # PCA
pca_panss <- principal(Y[, , 1], nfactors = 7, rotate= 'oblimin')
pca_panss <- principal(Y[, , 1], nfactors = 5, rotate= 'oblimin') 
summary(pca_panss)

# optimal number of principal components
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

##                                         # ICA
## ica_panss <- fastICA(Y[, , 1], n.comp = 5, alg.typ = "parallel", fun = "logcosh", alpha = 1, method = "C", row.norm = FALSE, maxit = 200, tol = 0.0001, verbose = TRUE)

## par(mfrow = c(1, 3))
## plot(ica_panss$X, main = "Pre-processed data")
## plot(ica_panss$S, main = "ICA components")
## plot(ica_panss$X %*% ica_panss$K, main = "PCA components")

## summary(ica_panss)

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
###                Exploratory   MIRT
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
