rm(list=ls())
library(dplyr)
library(survival)
library(ggplot2)
library(Rcpp)
library(rstan)
source("test4_survival/fnt_simulation.R")

alpha <- 100
## sample sizes from TCGA blca data
nobs <- 200 
ncen <- 2

K = 100 # number of covariates
K1= 10

beta = c(rep(1,K1), rep(0,(K-K1)))
x = matrix(rnorm(1:((nobs + ncen)*K), sd = 0.1), nrow=(nobs + ncen))
mu = x %*% beta
idx.event <- c(rep(TRUE, nobs), rep(FALSE, ncen))

## test these inputs for arbitrary values of alpha & mu
simulated_data <- sim_ph(alpha = alpha,
                         mu = mu,
                         idx.event = idx.event) 

stan_data <- list(
  Nobs = nobs,
  Ncen = ncen,
  yobs = simulated_data[idx.event, "os_months"],
  ycen = simulated_data[!idx.event, "os_months"],
  K = K,
  xobs = x[idx.event, ],
  xcen = x[!idx.event, ]
)

stan_file <- "test4_survival/coxph.stan"
## MODEL1 random initial values
recover_simulated <- stan(stan_file,
                          data = stan_data,
                          chains = 1,
                          pars = c("beta","alpha","lambda"), 
                          iter = 1000,
                          seed = 1328025050
)
xx=summary(recover_simulated)$summary
head(xx)
crossprod(summary(recover_simulated)$summary[1:100,1]-beta)
hist(xx[, "Rhat"])
plot(summary(recover_simulated)$summary[1:100,1])


