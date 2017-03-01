rm(list=ls())
library(rstan)
library(glmnet)
library(grplasso)

## Simulation - LASSO 

##  parameters
N = 30 # number of samples
K = 10 # number of parameters

x = matrix(rnorm(1:(N*K)), nrow=N)
beta = c(rep(1,4), rep(0,6))
y = as.vector(x %*% beta + rnorm(N, sd = 0.1))



fit <- stan(file = "test3_LASSO/lasso.stan", data = list(N,K, y,x, lambda=0.1), pars=c("beta"), chains = 1)
print(fit, digits=3)

fit2=glmnet(x,y,lambda=1)
fit2$beta


## Simulation - Group LASSO 
x = matrix(rnorm(1:(N*K)), nrow=N)
beta2 = c(rep(1,4), rep(-1,3), rep(0,3))
y = as.vector(x %*% beta2 + rnorm(N, sd = 0.1))



fit <- stan(file = "test3_LASSO/grplasso.stan", data = list(N,K, y,x, lambda=0.1), pars=c("beta"), chains = 1)
print(fit, digits=3)

fit3 = grplasso(x=cbind(1,x), y, index = c(NA,rep(1,4), rep(2,3), rep(3,3)), lambda=3, model = LinReg())
fit3$coefficients
