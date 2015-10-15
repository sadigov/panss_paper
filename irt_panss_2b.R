##########################################################################################
#
#   IRT Analysis of Bitopertin PANSS NSFS
#   Trial n25310
#   Author: Shamil Sadikhov
#   Date 10 July 2014
#
##########################################################################################

rm(list = ls())


# library("R2WinBUGS")
require(Hmisc)
library("R2OpenBUGS")
library("ltm")
library("mcmcplots")
require(plyr)
require(reshape2)
require(BRugs)
require(abind)

setwd('u:/My Documents/Statistics_projects/IRT/bitopertin')
Sys.setlocale("LC_TIME", "English") # with incorrect locale as.Date does not recognize months, months appear in German

# BUGs directory:
obd <- "C:/Users/sadikhos/Programs/OpenBUGS322/OpenBUGS.exe"
wbd <- "C:/Users/sadikhos/Programs/WinBUGS14"

# bitodata <- sasxport.get("u:/My Documents/Bitopertin/MAAP/PPNS/Analysis/n25310a/data_21112013/csv_21112013", method='csv')
# names(bitodata)
# panss <- bitodata$panss
load('panss.Rdata')
panss.nsfs <- subset(panss, qsscat %in% c('NEGATIVE SCALE', 'GENERAL PSYCHOPATHOLOGY SCALE') & qstestcd %in% c('QS027N01', 'QS027N02', 'QS027N03', 'QS027N04', 'QS027N06', 'QS027G07', 'QS027G16'))

# rm(bitodata)
# save(panss, file = 'panss.Rdata')

lgrmar2b <- function(){
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

 shift <- mean(theta[ , ])
 scale <- sd(theta[ , ])

for (i in 1:n){
    theta.raw[i, 1:T] ~ dmnorm(mu.theta[], Pr.theta[ , ])
    # theta[i, 1:T] ~ dmnorm(mu.theta[], Pr.theta[ , ])
    for (t in 1:T){

    theta[i, t] <- (theta.raw[i, t] - shift) / scale
}
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
    for (k in 1:K[j] - 1){
        kappa.star [j, k] ~ dnorm(m.kappa, pr.kappa)
        kappa[j, k] <- ranked(kappa.star[j, 1:(K[j] - 1)], k)
    }
}
pr.kappa <- pow(s.kappa, -2)       
}

## some temporary filename:
filename <- file.path(getwd(), "lgrmar2b.txt")
## write model file:
write.model(lgrmar2b, filename)
file.show(filename)

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

# Y <- Y[1:50, , 1:2] # note: WinBUGS runs with this subset
# Y <- Y[ , , 1:2]    # note: WinBUGS runs with this subset
# Y <- Y[ , , 1:3]    # ran in BUGS but slow
Y <- Y[ , , 1:4]    # ran in BUGS but slow


# Y <- Y[ , , c(1,2,5)]

# treatment arm
trt <- 3 - as.numeric(panss.wide2[panss.wide2$variable == '4_qsstresn', ]$trt01a)


# time
T <- dim(Y)[3]
# n.variables    
p <- dim(Y)[2]
# sample size
n <- dim(Y)[1]
# K - max level for ordinal variable for each item
K <- apply(apply(Y, c(2,3), max), 1, max)
# number of arms
n.trt <- length(levels(panss.wide2[panss.wide2$variable == '4_qsstresn', ]$trt01a))

m.alpha <- 1.0
s.alpha <- 2.5
m.kappa <- 0.0
s.kappa <- 2.5
m.mu.theta <- 0
s.mu.theta <- 1


data <- list("Y", "n", "p", "T", "K", "m.kappa", "s.kappa", "m.alpha", "s.alpha", "m.mu.theta", "s.mu.theta")

monitor <- c("alpha", "kappa", "theta", "rho")
n.burn <- 4000
n.thin <- 10
n.sim <- 500 * n.thin + n.burn

lgrmar2b.out <- R2OpenBUGS::bugs(data = data, inits = NULL,
                                 parameters.to.save = monitor,
                                 model.file = file.path(getwd(), "lgrmar2b.txt"),
                                 n.iter = n.sim, n.thin = n.thin, n.burnin = n.burn,
                                 OpenBUGS.pgm = obd,  working.directory = getwd(),
                                 debug = TRUE,
                                 codaPkg = FALSE)

lgrmar2.out <- R2WinBUGS::bugs(data = data, inits = NULL,
                                 parameters.to.save = monitor,
                                 model.file = file.path(getwd(), "lgrmar2.txt"),
                                 n.iter = n.sim, n.thin = n.thin, n.burnin = n.burn,
                                 bugs.directory = wbd,  working.directory = getwd(),
                                 debug = FALSE,
                                 codaPkg = FALSE)

plot(lgrmar2.out)

lgrmar2.out2 <- R2WinBUGS::bugs(data = data, inits = NULL,
                                 parameters.to.save = monitor,
                                 model.file = file.path(getwd(), "lgrmar2.txt"),
                                 n.iter = n.sim, n.thin = n.thin, n.burnin = n.burn,
                                 bugs.directory = wbd,  working.directory = getwd(),
                                 debug = FALSE,
                                 codaPkg = FALSE)

plot(lgrmar2.out2)
