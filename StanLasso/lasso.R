rm(list=ls())
library(rstan)
library(glmnet)
library(grplasso)

## Simulation - LASSO 

##  parameters
N = 100 # number of samples
K = 100 # number of parameters
K1= 10
lambda0 = 0.1

beta = c(rep(1,K1), rep(0,(K-K1)))
x = matrix(rnorm(1:(N*K)), nrow=N)
y = as.vector(x %*% beta + rnorm(N, sd = 0.1))



## non-informative prior 
fit <- stan(file = "test3_LASSO/lasso.stan", data = list(N,K, y,x, lambda=lambda0), pars=c("beta"), chains = 1)
# print(fit, digits=3)
crossprod(summary(fit)$summary[1:100,1]-beta)

## Laplace prior, lambda as a CONSTANT
fit2 <- stan(file = "test3_LASSO/lasso2.stan", data = list(N,K, y,x, lambda=lambda0), pars=c("beta","mu", "sigma"), chains = 1)
crossprod(summary(fit2)$summary[1:100,1]-beta)

## Laplace prior, lambda as a PARAMETER
fit3 <- stan(file = "test3_LASSO/lasso3.stan", data = list(N,K, y,x), pars=c("beta","mu", "sigma","lambda"), chains = 1)
crossprod(summary(fit3)$summary[1:100,1]-beta)

## non-informative prior  
fit.linear <- stan(file = "test3_LASSO/linear.stan", data = list(N,K, y,x), pars=c("beta"), chains = 1)
crossprod(summary(fit.linear)$summary[1:100,1]-beta)

## weakly informative normal prior distribution
fit.linear2 <- stan(file = "test3_LASSO/linear2.stan", data = list(N,K, y,x), pars=c("beta"), chains = 1)
crossprod(summary(fit.linear2)$summary[1:100,1]-beta)

## normal prior distribution
fit.linear3 <- stan(file = "test3_LASSO/linear3.stan", data = list(N,K, y,x), pars=c("beta"), chains = 1)
crossprod(summary(fit.linear3)$summary[1:100,1]-beta)

## ridge regression, lambda as a PARAMETER
fit.ridge <- stan(file = "test3_LASSO/ridge.stan", data = list(N,K, y,x), pars=c("beta","lambda"), chains = 1)
crossprod(summary(fit.ridge)$summary[1:100,1]-beta)

## weak normal prior ridge regression, lambda as a PARAMETER
fit.ridge2 <- stan(file = "test3_LASSO/ridge2.stan", data = list(N,K, y,x), pars=c("beta","lambda"), chains = 1)
crossprod(summary(fit.ridge2)$summary[1:100,1]-beta)

## weak normal prior ridge regression, lambda as a PARAMETER
fit.ridge3 <- stan(file = "test3_LASSO/ridge3.stan", data = list(N,K, y,x), pars=c("beta","lambda"), chains = 1)
crossprod(summary(fit.ridge3)$summary[1:100,1]-beta)



## GLM
fit2=glmnet(x,y,lambda = lambda0)
fit2$beta
crossprod(fit2$beta-beta)


