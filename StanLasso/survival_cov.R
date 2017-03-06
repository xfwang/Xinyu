library(dplyr)
library(survival)
library(ggplot2)
library(Rcpp)

sim_data <- function(alpha, mu, Nobs, Ncen) {
  observed_data <- data.frame(os_status = rep_len('DECEASED', Nobs),
                              os_months = rweibull(n = Nobs, alpha, exp(-(mu)/alpha)),
                              stringsAsFactors = F
  )
  
  censored_data <- data.frame(os_status = rep_len('LIVING', Ncen),
                              os_months = runif(Ncen) * rweibull(Ncen, alpha, exp(-(mu)/alpha)),
                              stringsAsFactors = F
  )
  
  return(observed_data %>% bind_rows(censored_data))
}

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
