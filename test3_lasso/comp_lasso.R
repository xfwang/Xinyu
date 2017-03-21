library(rstan); library(mvtnorm); library(dplyr)
sapply(paste("./functions", dir(paste0("./functions")),sep="\\"),source)

## COMPARE weak prior and strong prior in LASSO
## Hypothesis: weak prior -> high variance in posterior draws -> higher probability in covering 0 -> sparser model

## simulation
N = 100 # number of samples
K = 100 # number of covariates
K1= 10

beta = c(rep(1,K1), rep(0,(K-K1)))
x = cov_groupLasso(N, K)
y = as.vector(x %*% beta + rnorm(N, sd = 1))
res <- list()

res$fit1 <- stan(file = "test3_LASSO/lasso3.stan", data = list(N,K, y,x), pars=c("beta","mu", "sigma","lambda"), chains = 1)
res$fit2 <- stan(file = "test3_LASSO/lasso3_normal.stan", data = list(N,K, y,x), pars=c("beta","mu", "sigma","lambda"), chains = 1)
res$fit3 <- stan(file = "test3_LASSO/lasso3_normal_gamma.stan", data = list(N,K, y,x), 
                 pars=c("beta","mu", "sigma","lambda","tau"), chains = 1)
res$fit4 <- stan(file = "test3_LASSO/lasso3_normal_gamma2.stan", data = list(N,K, y,x), 
                 pars=c("beta","mu", "sigma","lambda","tau"), chains = 1)

lapply(res, function(fit) crossprod(param.stan(fit, confidence.level = 0.99)-beta))
lapply(res, function(fit) which(param.stan(fit, confidence.level = 0.99) != 0))


##lambda
lapply(res, function(fit) summary(fit)$summary["sigma", ])
lapply(res, function(fit) summary(fit)$summary["lambda", ])


