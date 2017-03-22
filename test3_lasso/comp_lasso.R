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
                 pars=c("beta","mu", "sigma","lambda"), chains = 1)
res$fit4 <- stan(file = "test3_LASSO/lasso3_normal_gamma2.stan", data = list(N,K, y,x), 
                 pars=c("beta","mu", "sigma","lambda","tau"), chains = 1)
res$fit5 <- stan(file = "test3_LASSO/lasso3_normal_gamma3.stan", data = list(N,K, y,x), 
                 pars=c("beta","mu", "sigma","lambda","tau"), chains = 1)

lapply(res, function(fit) crossprod(param.stan(fit, confidence.level = 0.99)-beta))
lapply(res, function(fit) which(param.stan(fit, confidence.level = 0.99) != 0))


##lambda
lapply(res, function(fit) summary(fit)$summary["sigma", ])
lapply(res, function(fit) summary(fit)$summary["lambda", ])


## -------------------------------------------------------------------------------------------------------------- ##
## Group LASSO simulation
N = 100 # number of samples
G = 10 # number of grouops
G1 = 3 # number of active grouops
K = 10 # number of covariates in each group

# group indicator of length K: 1, 2, 3... represent group No.; 0 represents inactive group
ind <- rep(1:K, each = G)

beta <- rep(c(1,-1,0.5), each = K) %>% 
  c(rep(0, (G-G1)*K))            # for Group LASSO

x = cov_groupLasso(N, K*G)
y = as.vector(x %*% beta + rnorm(N, sd = 2))

res <- list()
res$fit1 <- stan(file = "test3_LASSO/grplasso4.stan", data = list(N,K,G, y,x), pars=c("beta","mu","lambda"), chains = 1)
res$fit2 <- stan(file = "test3_LASSO/grplasso4_normal.stan", data = list(N,K,G, y,x), pars=c("beta","lambda"), chains = 1)
res$fit3 <- stan(file = "test3_LASSO/grplasso4_normal_gamma.stan", data = list(N,K,G, y,x), 
                 pars=c("beta","mu","lambda"), chains = 1)
res$fit4 <- stan(file = "test3_LASSO/grplasso4_normal_gamma2.stan", data = list(N,K,G, y,x), 
                 pars=c("beta","mu","lambda"), chains = 1)

## -------------------------------------------------------------------------------------------------------------- ##
## sparse Group LASSO simulation
beta1 <- c(-3:3, rep(0,3)) / 3 # for sparse Group LASSO
beta <- rep(beta1, G1) %>%
  c(rep(0, (G-G1)*K))

# x = matrix(rnorm(1:(N*K*G)), nrow=N)
x = cov_groupLasso(N, K*G)
y = as.vector(x %*% beta + rnorm(N, sd = 2))

res <- list()
res$fit.sparse1 <- stan(file = "test3_LASSO/sparseGrpLasso3.stan", data = list(N,K,G, y,x), pars=c("beta","lambda"), chains = 1)
res$fit.sparse2 <- stan(file = "test3_LASSO/sparseGrpLasso4.stan", data = list(N,K,G, y,x), pars=c("beta","lambda"), chains = 1)
res$fit.sparse3 <- stan(file = "test3_LASSO/sparseGrpLasso4_normal.stan", data = list(N,K,G, y,x), pars=c("beta","lambda"), chains = 1)
res$fit.sparse4 <- 
