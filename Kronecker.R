#######################################################
### Kronecker product                               ###
### Stan code to be put inside the function section ###
#######################################################
'
functions{
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
'

#######################################################
### Testing the kronecker product function          ###
#######################################################
rm(list = ls())
library(rstan)
# Set Seed
set.seed(1234)

model <- '
functions{
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
  int<lower=0> J; // number of schools 
  real y[J]; // estimated treatment effects
  real<lower=0> sigma[J]; // s.e. of effect estimates
  matrix[2,2] X;
  matrix[3,3] Y;
}

parameters {
  real mu; 
  real<lower=0> tau;
  real eta[J];
}

transformed parameters {
  real theta[J];
  for (j in 1:J)
  theta[j] <- mu + tau * eta[j];
}

model {
  eta ~ normal(0, 1);
  y ~ normal(theta, sigma);
}

generated quantities {
  matrix[rows(X)*rows(Y), cols(X)*cols(Y)] Z;
  Z <- kronecker(X, Y);
}
'

model_comp <- stan_model(model_code = model, verbose=FALSE)

X <- matrix(c(1,2,3,5), nrow = 2)
Y <- matrix(c(7,11,13,17,19,23,29,31,37), nrow = 3)
schools_dat <- list(J = 8, 
                    y = c(28,  8, -3,  7, -1,  1, 18, 12),
                    sigma = c(15, 10, 16, 11,  9, 11, 10, 18),
                    X = X,
                    Y = Y)


fit <- sampling(model_comp,
                data = schools_dat, 
                iter = 1000, 
                chains = 1)

# print(fit)
la <- extract(fit, permuted = TRUE) # return a list of arrays
Z <- la$Z
Z <- apply(Z, c(2,3), mean)

TrueZ <- kronecker(X, Y) # Done by an R function
TrueZ - Z