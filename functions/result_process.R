library(dplyr);library(rstan)

#####---------------------------------------------------------------------------------------------------------#####
##### R script to shrink parameters using 95% CI calculated in stan()                                         #####
#####                                                                                                         #####
##### Input:                                                                                                  #####
#####       stan.ext          :=  the matrix of parameter(s) iteration record                                 #####
#####       confidence.level  :=  confidence level of confidence interval, default is 0.95                    #####
##### Output:                                                                                                 #####
#####       idx <= 0          :=  a vector of indicators on if a parameter is 0 or not                        #####
#####---------------------------------------------------------------------------------------------------------#####
shrink_ind.stan <- function(stan.ext, confidence.level = 0.95) {
  idx <- stan.ext %>% 
    apply(2, quantile, c((1 - confidence.level) / 2, (1 + confidence.level) / 2)) %>%
    apply(2, prod)
  return((idx <= 0))
}

shrink.stan <- function(stan.obj, param = "beta", confidence.level = 0.95) {
  beta_all <- param.stan(stan.obj, param = param)
  beta <- colMeans(beta_all) * as.numeric(!shrink_ind.stan(beta_all, confidence.level = confidence.level))
  return(beta)
}

param.stan <- function(stan.obj, param = "beta") {
  la <- rstan::extract(stan.obj, permuted = T)
  beta <- la[[param]]
  return(beta)
}

#####---------------------------------------------------------------------------------------------------------#####
##### R script to make plots on a list of stan objects                                                        #####
#####                                                                                                         #####
##### Input:                                                                                                  #####
#####       stan.list         :=  a list of stan objects                                                      #####
#####       type              :=  the type of plots                                                           #####
##### Output:                                                                                                 #####
#####       output            :=                                                                              #####
#####---------------------------------------------------------------------------------------------------------#####
stan_plot <- function(stan.list, type = c("Rhat")) {
  pp = lapply(stan.list, stan_rhat)
  return(pp)
}

#####---------------------------------------------------------------------------------------------------------#####
##### R script to find out the rank of shrinkage in the coefficients in stan                                  #####
##### ASSUMPTION: normal distribution of beta                                                                 #####
#####                                                                                                         #####
##### Input:                                                                                                  #####
#####       beta.list         :=  the series of sampled posterior coefficients for one parameter              #####
#####       beta.array        :=  the array of sampled posterior coefficients for several parameters (niter*P)#####
##### Output:                                                                                                 #####
#####       credible_level    :=  how likely the true parameter is 0, which is the cumulative probability of 0#####
#####---------------------------------------------------------------------------------------------------------#####
cred_prob <- function(beta.list) {
  cum_prob <- pnorm(0, mean = mean(beta.list), sd = sd(beta.list))
  return(min(cum_prob, 1-cum_prob)*2)
}

cred_quantile <- function(beta.list) {
  cut_perc <- mean(beta.list < 0)
  return(min(cut_perc, 1-cut_perc)*2)
}

rank_select <- function(stan.obj, param = "beta") {
  beta.array <- param.stan(stan.obj, param = param)
  return(apply(beta.array, 2, cred_quantile))
}
