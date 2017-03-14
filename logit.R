rm(list=ls())
library(rstan)
library(glmnet)
library(rmutil)

## Simulation - LASSO 

##  parameters
N = 100 # number of samples
K = 100 # number of covariates
K1= 10
lambda0 = 0.1

# beta = c(rep(1,K1), rep(0,(K-K1)))
beta <- rlaplace(100)
x = matrix(rnorm(1:(N*K)), nrow=N)
p = 1 / (1 + exp( - x %*% beta))
y = sapply(p, rbinom, n = 1, size = 1)



## non-informative prior 
fit <- stan(file = "test5_logistic/logit.stan", data = list(N,K, y,x, lambda=lambda0), pars=c("beta"), chains = 1)
# print(fit, digits=3)
crossprod(summary(fit)$summary[1:100,1]-beta)


## GLM
cvfit <- cv.glmnet(x,y, family = "binomial")
plot(cvfit) # lambda.min & lambda.1se
#final model
fit2 <- glmnet(x,y, family = "binomial",lambda = cvfit$lambda.1se)
fit2$beta
crossprod(fit2$beta-beta)


