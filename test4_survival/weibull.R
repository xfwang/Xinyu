rm(list=ls())
library(dplyr)
library(survival)
library(ggplot2)
library(Rcpp)
source("test4_survival/functions.R")

#####  -----    COMMENTS   ------  #####
## mu = -3 => large variance => small effect size in estimation
## mu = -1 => small variance => relatively easy to estimate alpha
test_alpha <- 0.8
test_mu <- -3

## sample sizes from TCGA blca data
test_nobs <- 179 
test_ncen <- 230

## test these inputs for arbitrary values of alpha & mu
simulated_data <- sim_data(alpha = test_alpha,
                           mu = test_mu,
                           Nobs = test_nobs,
                           Ncen = test_ncen) 

## plot KM curve from simulated data
simulated_data <- 
  simulated_data %>%            ##  x %>% f(y) turns into f(x, y)
  dplyr::mutate(os_deceased = os_status == 'DECEASED')


observed_data <- simulated_data %>%
  dplyr::filter(os_status == 'DECEASED')

censored_data <- simulated_data %>%
  dplyr::filter(os_status != 'DECEASED')

stan_data <- list(
  Nobs = nrow(observed_data),
  Ncen = nrow(censored_data),
  yobs = observed_data$os_months,
  ycen = censored_data$os_months
)
rm(censored_data)
rm(observed_data)
str(stan_data)


## MODEL1 random initial values
# strong prior on alpha
stan_file <- "test4_survival/weibull.stan"
recover_simulated <- 
  rstan::stan(stan_file,
              data = stan_data,
              chains = 1,
              iter = 1000,
              seed = 1328025050
  )
print(recover_simulated)

# weak prior on alpha
stan_file2 <- "test4_survival/weibull2.stan"
recover_simulated2 <- 
  rstan::stan(stan_file2,
              data = stan_data,
              chains = 1,
              iter = 1000,
              seed = 1328025050
  )
print(recover_simulated2)

# strong prior on both alpha and mu
stan_file3 <- "test4_survival/weibull3.stan"
recover_simulated3 <- 
  rstan::stan(stan_file3,
              data = stan_data,
              chains = 1,
              iter = 1000,
              seed = 1328025050
  )
print(recover_simulated3)

## MODEL2 set initial values
# recover_simulated2 <- 
#   rstan::stan(stan_file,
#               data = stan_data,
#               chains = 4,
#               iter = 1000,
#               init = gen_inits
#   )
# print(recover_simulated2)

## CONVERGENCE
rstan::traceplot(recover_simulated2, 'lp__')
rstan::traceplot(recover_simulated2, c('alpha','mu'), ncol = 1)

## DIAGNOSIS 
if (interactive())
  shinystan::launch_shinystan(recover_simulated2)


