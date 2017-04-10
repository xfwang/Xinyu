rm(list=ls())
library(rstan)
library(glmnet)
library(grplasso)
library(SGL)
library(magrittr)
library(mvtnorm)
sapply(paste("./functions", dir(paste0("./functions")),sep="\\"),source)

## Simulation - Group LASSO 

##  parameters
N = 100 # number of samples
G = 10 # number of grouops
G1 = 3 # number of active grouops
K = N/G # number of covariates in each group
lambda0 = 0.1

# group indicator of length K: 1, 2, 3... represent group No.; 0 represents inactive group
ind <- rep(1:G, each = K)

beta1 <- c(-2:2, rep(0,5))
beta <- rep(beta1, G1) %>% 
  c(rep(0, (G-G1)*K))

# x = matrix(rnorm(1:(N*K*G)), nrow=N)
x = X_groupLasso(N, K*G)
y = as.vector(x %*% beta + rnorm(N, sd = 2))


res <- list()
## group lasso, without any prior
res$fit <- stan(file = "test3_LASSO/grplasso.stan", data = list(N,K,G, y,x,lambda=lambda0), pars=c("beta"), chains = 1)
## group lasso, Laplace prior
res$fit2 <- stan(file = "test3_LASSO/grplasso2.stan", data = list(N,K,G, y,x,lambda=lambda0), pars=c("beta"), chains = 1)
## group lasso, Laplace prior, lambda as a PARAMETER
res$fit3 <- stan(file = "test3_LASSO/grplasso3.stan", data = list(N,K,G, y,x), pars=c("beta","lambda"), chains = 1)
## group lasso, Laplace prior, lambda as a PARAMETER, with hyperparameters
res$fit4 <- stan(file = "test3_LASSO/grplasso4.stan", data = list(N,K,G, y,x), pars=c("beta","lambda", "mu"), chains = 1)
## group lasso, Gasussian prior, lambda as a PARAMETER, with hyperparameter
res$fit5 <- stan(file = "test3_LASSO/grplasso4_normal.stan", data = list(N,K,G, y,x), pars=c("beta","lambda", "mu"), chains = 1)
## group lasso, Laplace (normal_gamma) prior, lambda as a PARAMETER, with hyperparameter
res$fit6 <- stan(file = "test3_LASSO/grplasso4_normal_gamma.stan", data = list(N,K,G, y,x), pars=c("beta","lambda", "mu"), chains = 1)
## group lasso, t prior 
res$fit7 <- stan(file = "test3_LASSO/grplasso4_t.stan", data = list(N,K,G, y,x,nu=3.0), pars=c("beta","lambda"), chains = 1)

## Sparse group lasso, without any prior
res$fit.sparse1 <- stan(file = "test3_LASSO/sparseGrpLasso.stan", data = list(N,K,G, y,x,lambda=lambda0), pars=c("beta"), chains = 1)
## Sparse group lasso, Laplace prior
res$fit.sparse2 <- stan(file = "test3_LASSO/sparseGrpLasso2.stan", data = list(N,K,G, y,x,lambda=lambda0), pars=c("beta"), chains = 1)
## Sparse group lasso, Laplace prior, lambda as a PARAMETER
res$fit.sparse3 <- stan(file = "test3_LASSO/sparseGrpLasso3.stan", data = list(N,K,G, y,x,lambda=lambda0), pars=c("beta","lambda"), chains = 1)
## Sparse group lasso, Laplace prior, lambda as a PARAMETER
res$fit.sparse4 <- stan(file = "test3_LASSO/sparseGrpLasso4.stan", data = list(N,K,G, y,x,lambda=lambda0), pars=c("beta","lambda"), chains = 1)


## PARAMETER ESTIMATION 
lapply(res, function(fit) crossprod(shrink.stan(fit, confidence.level = 0.99)-beta))
lapply(res, function(fit) which(shrink.stan(fit, confidence.level = 0.99) != 0))

lapply(res, function(fit) crossprod(summary(fit)$summary[1:100,1]-beta))
lapply(res, function(fit) table((param.stan(fit, confidence.level = 0.99)==0), (beta==0)))

## PACKAGE: grplasso
fit.glasso = grplasso(x=cbind(1,x), y, index = c(NA,ind), lambda=100, model = LinReg())
crossprod(fit.glasso$coefficients[-1]-beta)

## PACKAGE: SGL
cv.sgl <- list(x = x, y = matrix(y, ncol=1)) %>% 
  cvSGL(ind, type = "linear")
# largest value of lambda such that error is within 1 standard error of the minimum
lambda.1se <- cv.sgl$lambdas[which(cv.sgl$lldiff == max(cv.sgl$lldiff[cv.sgl$lldiff < (cv.sgl$llSD[which.min(cv.sgl$lldiff)] + min(cv.sgl$lldiff))]))]
fit.sgl <- list(x = x, y = matrix(y, ncol=1)) %>% 
  SGL(ind, type = "linear", lambda = lambda.1se)
crossprod(fit.sgl$beta-beta)

