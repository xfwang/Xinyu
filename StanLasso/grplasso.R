rm(list=ls())
library(rstan)
library(glmnet)
library(grplasso)

## Simulation - Group LASSO 

##  parameters
N = 100 # number of samples
K = 100 # number of parameters
K1= 10
lambda0 = 0.1
# group indicator of length K: 1, 2, 3... represent group No.; 0 represents inactive group
ind <- c(rep(1,K1), rep(2,K1), rep(0,(K-K1*2)))
n.levels <- length(unique(ind))

beta <- rep(0, K); beta[ind == 1] <- 1; beta[ind == 2] <- -2;
ind.sparse <- c(4:6, 11:13); beta[ind.sparse] <- 0

x = matrix(rnorm(1:(N*K)), nrow=N)
y = as.vector(x %*% beta + rnorm(N, sd = 1))

# write model
# cat("
# data {
#     int<lower=0> N;
#     int<lower=1> K;
#     vector[N] y;
#     matrix[N,K] x;
#     vector[K] ind; // group indicator
#     real<lower=0> lambda;
#     }
#     
#     parameters {
#     vector[K] beta;
#     vector[", sum(ind==1), "] ind1;
#     vector[", sum(ind==2), "] ind2;
#     vector[", sum(ind==0), "] ind3;
# 
#     ind1 = ", paste("{",paste(which(ind==1),collapse=","),"};"),"
#     ind2 = ", paste("{",paste(which(ind==2),collapse=","),"};"),"
#     ind3 = ", paste("{",paste(which(ind==0),collapse=","),"};"),"
#     }
#     
#     transformed parameters {
#     real<lower=0> squared_error;
#     real<lower=0> S1;
#     real<lower=0> S2;
#     real<lower=0> S3;
#     
#     squared_error = dot_self(y - x * beta);
#     S1 = dot_self(beta[ind1]);
#     S2 = dot_self(beta[ind2]);
#     S3 = dot_self(beta[ind3]);",
#     "
#     }
# 
# model {
#   // beta ~ double_exponential(mu, sigma);
#   target += -squared_error;
#   target += - lambda * (S1+S2+S3); 
# }
# 
# generated quantities {
#   real<lower=0> sigma_squared;
#   sigma_squared = squared_error / N;
# }
# ",
# file = "test3_LASSO/grplasso.stan");


## group lasso, without any prior
fit <- stan(file = "test3_LASSO/grplasso.stan", data = list(N,K, y,x,lambda=lambda0), pars=c("beta"), chains = 1)
crossprod(summary(fit)$summary[1:100,1]-beta)

## group lasso, Laplace prior
fit2 <- stan(file = "test3_LASSO/grplasso2.stan", data = list(N,K, y,x,lambda=lambda0), pars=c("beta"), chains = 1)
crossprod(summary(fit2)$summary[1:100,1]-beta)

## group lasso, Laplace prior, lambda as a PARAMETER
fit3 <- stan(file = "test3_LASSO/grplasso3.stan", data = list(N,K, y,x), pars=c("beta","lambda"), chains = 1)
crossprod(summary(fit3)$summary[1:100,1]-beta)

## Sparse group lasso, without any prior
fit.sparse1 <- stan(file = "test3_LASSO/sparseGrpLasso.stan", data = list(N,K, y,x,lambda=lambda0), pars=c("beta"), chains = 1)
crossprod(summary(fit.sparse1)$summary[1:100,1]-beta)

## Sparse group lasso, Laplace prior
fit.sparse2 <- stan(file = "test3_LASSO/sparseGrpLasso2.stan", data = list(N,K, y,x,lambda=lambda0), pars=c("beta"), chains = 1)
crossprod(summary(fit.sparse2)$summary[1:100,1]-beta)

## Sparse group lasso, Laplace prior, lambda as a PARAMETER
fit.sparse3 <- stan(file = "test3_LASSO/sparseGrpLasso3.stan", data = list(N,K, y,x,lambda=lambda0), pars=c("beta"), chains = 1)
crossprod(summary(fit.sparse3)$summary[1:100,1]-beta)


fit.glasso = grplasso(x=cbind(1,x), y, index = c(NA,ind), lambda=3, model = LinReg())
crossprod(fit.glasso$coefficients-beta)





