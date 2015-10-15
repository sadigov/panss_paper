################################################################################
################################################################################
##                                                                            ##
##   Title:  Generation of longitudinal multiple-group graded response data   ##
##   Author: Shamil Sadikhov                                                  ##
##   Date:   August 20, 2015                                                  ##
##                                                                            ##
################################################################################
################################################################################

require(mvtnorm)
################################################################################
# Vector of probabilities where each row gives the probability of answering in #
# the corresponding category                                                   #  
################################################################################
irf <- function(theta = NULL                     # ability parameter
                 ,a = NULL                       # item discrimination
                 ,b = NULL                       # item difficulty (if dichotomous data) or vector of thresholds (if graded response data)
                 ,cc = 0                         # guessing parameter
                 ,type = "3pl"){      # 3pl use the logistic function as the link function, norm is the normal standard CDF
  
  # Dichotomous data
  if (length(b) == 1) {
    if (type == "3pl")  prob <- cc+(1-cc)/(1+exp(-a*theta + b))
    if (type == "norm") prob <- cc+(1-cc)*pnorm(a*theta - b)
    return(prob)
  }
  
  # Graded response data (length(b)>1)
  else {
    ncat = length(b) + 1 # number of categories
    prob <- numeric(ncat)
    probc <- numeric(ncat) # cumulative probabilities
    if (type == "3pl"){
      probc[1] <- 1
      for (i in 2:ncat) {
        probc[i] <- cc + (1-cc)/(1+exp(-a*theta + b[i-1])) # probc[c] = P(y > c-1|theta,a,b) (=P(y >= c|theta,a,b))
      } 
    }
    if (type == "norm"){
      probc[1] <- 1
      for (i in 2:ncat) {
        probc[i] <- cc + (1-cc)*pnorm(a*theta - b[i-1]) # probc[c] = P(y > c-1|theta,a,b)
      }
    }
    prob <- probc - c(probc[-1], 0) # prob[c] = P(y=c|theta,a,b) = P(y > c-1|theta,a,b) - P(y > c|theta,a,b)
    return(prob)
    }
}

################################################################################
# Generates item responses for graded/dichotomous, multiple items, multiple    #
# persons, multiple groups (one or two) BUT single time                        #
################################################################################
irtgen <- function(n.s = 10                 # number of subjects
                    ,avec = 1               # discrimination parameter vector
                    ,bmat = 0               # threshold/difficulty parameter matrix
                    ,cvec = c(0)            # guessing parameter
                    ,type = "3pl"
                    ,mtheta = c(-0.5, 0.5)  # ability parameter mean
                    ,sdtheta = 1){          # ability parameter sd
  
  ncat <- ncol(bmat) + 1               # number of categories, starting with 1,...
  nitem <-  length(avec)
  data <- matrix(0, nrow = n.s, ncol = nitem)
  
  cut.p <- floor(n.s/2) # Cut the group in two half
  if (length(mtheta) > 1) {
    theta1 <- rnorm(cut.p, mean = mtheta[1], sd = sdtheta)
    theta2 <- rnorm(n.s-cut.p, mean = mtheta[2], sd = sdtheta)
    theta <- c(theta1, theta2) # Generate a vector of abilities
  } 
  else {
    theta <- rnorm(n.s, mean = mtheta, sd = sdtheta)
  } # Generate a vector of abilities
  
  if (ncat>2){
    for (i in 1:nitem) {
      p <- t(apply(as.array(theta), 1, FUN = function(x) irf(theta = x, a = avec[i], b = bmat[i, ], c = cvec[i], type = type)))
      data[, i] <- t(apply(p, 1, function(x) sample(1:ncat, size = 1, replace = TRUE, prob = x)))
    }
  } else {
    for (i in 1:nitem) {
      p <- t(apply(as.array(theta), 1, FUN = function(x) irf(theta = x, a = avec[i], b = bmat[i, ], c = cvec[i], type = type)))
      data[, i] <- t(apply(p, 1, function(x) rbinom(n.s, 1, prob = x)))
    }
  }
  
  list <- list(data, theta)
  return(list)
}

################################################################################
# Generates longitudinal item responses for graded/dichotomous, multiple       #
# items, multiple persons, multiple groups (one or two)                        #
################################################################################
lirtgen <- function( n.s = 10               # number of subjects
                     ,n.times = 7
                     ,avec = c(1)             # discrimination parameter
                     ,bmat = c(0)             # difficulty parameter
                     ,cvec = c(0)             # guessing parameter
                     ,type = "3pl"
                     ,mtheta = c(-0.5, 0.5)   # ability parameter mean vector
                     ,trend = c("linear")    # Either linear or quadratic
                     ,cov.struct = "AR1"    # Can be 1D, unif, Markov or AR1
                     ,rho = 0.5
                     ,sigma = 1
                     ,t.obs = c(1:n.times)
                     ,format = "array"){
  n.item <- nrow(avec) # number of items
  n.cat <- ncol(bmat)+1 # Number of categories
  cut.n <- floor(n.s/2)
  
  ## Covariance structure
  if (cov.struct == "1D") { # 1-dependent correlation matrix
    vcov <- diag(sigma^2, n.times)
    for (t in 2:n.times) {
      vcov[t-1,t] <- sigma^2*rho
      vcov[t,t-1] <- sigma^2*rho
    }
  } else if (cov.struct == "unif") { # Uniform
    vcov <- matrix(rho*sigma^2,n.times,n.times)
    diag(vcov) <- sigma^2
  } else if (cov.struct == "Markov") { # Markov
    vcov <- matrix(0, n.times, n.times)
    vcov[1,1] <- sigma^2
    for (t in 2:n.times) {
      vcov[t,t] <- sigma^2
      for (j in 1:(t-1)) {
        d_tj <- abs(t.obs[t]-t.obs[j])
        vcov[t,j] <- sigma^2*(rho^d_tj)
        vcov[j,t] <- vcov[t,j]
      }
    }
  } else if (cov.struct == "AR1"){ # AR1
    vcov <- matrix(0, n.times, n.times)
    vcov[1,1] <- sigma^2
    for (t in 2:n.times) {
      vcov[t,t] <- sigma^2
      for (j in 1:(t-1)) {
        vcov[t,j] <- sigma^2*(rho^(t-j))
        vcov[j,t] <- vcov[t,j]
      }
    }
  } else {
    stop('Covariance structure not found. Can be either "1D", "unif", "Markov" or "AR1". Default is "AR1".')
  }
  
  ## Time-dependent mean of the ability parameter
  if (length(mtheta)>1){ # If more than one group
    if (trend=="linear") {
      mean.vec1 <- c(0, 1:(n.times-1)*mtheta[1]/(n.times-1))
      mean.vec2 <- c(0, 1:(n.times-1)*mtheta[2]/(n.times-1))
    } else if (trend == "quadratic") {
      x.mat <- rbind(rep(0, n.times),
                     0:(n.times-1),
                     (0:(n.times-1))^2)
      beta1.1 <- 2*mtheta[1]/(n.times-1)
      beta2.1 <- (mtheta[1]-(n.times-1)*beta1.1)/(n.times-1)^2
      mean.vec1 <- t(c(0, beta1.1, beta2.1))%*%x.mat
      beta1.2 <- 2*mtheta[2]/(n.times-1)
      beta2.2 <- (mtheta[2]-(n.times-1)*beta1.2)/(n.times-1)^2
      mean.vec2 <- t(c(0, beta1.2, beta2.2))%*%x.mat
    }
    
    theta1 <- mvtnorm::rmvnorm(n= cut.n, mean = mean.vec1, sigma = vcov)
    theta2 <- mvtnorm::rmvnorm(n= n.s - cut.n, mean = mean.vec2, sigma = vcov)
    theta <- rbind(theta1, theta2)
  } else { # Only one group
    if (trend=="linear") {
      mean.vec <- c(0, 1:(n.times-1)*mtheta/(n.times-1))
    } else if (trend == "quadratic") {
      x.mat <- rbind(rep(0, n.times),
                     0:(n.times-1),
                     (0:(n.times-1))^2)
      beta1 <- 2*mtheta/(n.times-1)
      beta2 <- (mtheta-(n.times-1)*beta1)/(n.times-1)^2
      mean.vec <- t(c(0, beta1, beta2))%*%x.mat
    }
    
    theta <- mvtnorm::rmvnorm(n= n.s, mean = mean.vec, sigma = vcov)
  }
  # Theta has n.s rows and n.times columns
  
  # Pre-allocate:
  data <- array(0,dim=c(n.s,n.item,n.times))
  
  # We now generate the data
  if (n.cat>2){ # Graded
    for (t in 1:n.times) {
      for (i in 1:n.item) {
        p <- t(apply(as.array(theta[,t]), 1, FUN = function(x) irf(theta = x, a = avec[i], b = bmat[i, ], c = cvec[i], type = type)))
        data[,i,t] <- t(apply(p, 1, function(x) sample(1:n.cat, size = 1, replace = TRUE, prob = x)))
      }
    }
  } else { # Dichotomous
    for (t in 1:n.times) {
      for (i in 1:n.item) {
        p <- t(apply(as.array(theta[,t]), 1, FUN = function(x) irf(theta = x, a = avec[i], b = bmat[i, ], c = cvec[i], type = type)))
        data[,i,t] <- t(apply(p, 1, function(x) rbinom(n.s, 1, prob = x)))
      }
    }
  }
  
  if (format == "long") {
    data <- c(data)
    
    ii <- rep(c(1:n.s), n.times*n.item)
    
    jj <- rep(0, n.s*n.item)
    for (i in 1:n.item) {
      jj[(n.s*(i-1)+1):(n.s*i)] <- rep(i, n.s)
    }
    jj <- rep(jj, n.times)
    
    tt <- rep(0, n.s*n.item*n.times)
    for (i in 1:n.times) {
      tt[(n.s*n.item*(i-1)+1):(n.s*n.item*i)] <- rep(i, n.s*n.item)
    }
    
    list <- list(data, theta, ii = ii, jj = jj, tt = tt)
    
    return(list)
    
  } else {
    list <- list(data, theta)
    
    return(list)
    
  }
}

# # TEST
# n.item <- 5
# n.s <- 6
# n.item_cat <- 7
# n.times <- 7
# avec <- matrix(rlnorm(n.item, meanlog=.2, sdlog=.2))
# avec <- matrix(3,n.item)
# bmat <- matrix(rnorm(n.item*(n.item_cat - 1)), n.item)
# bmat <- t(apply(bmat, 1, sort, decreasing = FALSE)) # Thresholds at each time point -> matrix of size (n.item)x(n.itemLevels-1)
# if(n.item_cat == 2){
#  bmat <- t(bmat)
# }
# cvec <- matrix(rep(0, n.item))
# mtheta = c(-0.5, 0.5)
# mtheta = 0.5
# type = "norm"
# trend <- "quadratic"
# cov.struct <- "AR1"
# rho <- 0#.5
# sigma <- 1
# format <- "long"
# 
# result <- lirtgen(n.s = n.s, n.times = n.times, avec = avec, bmat = bmat, cvec = cvec, type = type,
#                   mtheta = mtheta, trend = trend, cov.struct = cov.struct, rho = rho, sigma = sigma, format = format)
# theta_mean_1 <- apply(result[[2]][1:10000,], 2, mean)
# theta_mean_2 <- apply(result[[2]][10001:20000,], 2, mean)
# str(result)
# str(result[[1]])
# result[[1]][1,,]