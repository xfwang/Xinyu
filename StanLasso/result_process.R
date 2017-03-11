#####---------------------------------------------------------------------------------------------------------#####
##### R script to shrink parameters using 95% CI calculated in stan()                                         #####
#####                                                                                                         #####
##### Input:                                                                                                  #####
#####       stan.ext          :=  the matrix of parameter(s) iteration record                                 #####
#####       confidence.level  :=  confidence level of confidence interval, default is 0.95                    #####
##### Output:                                                                                                 #####
#####       output            :=  a vector of indicators on if a parameter is 0 or not                        #####
#####---------------------------------------------------------------------------------------------------------#####
shrink.stan <- function(stan.ext, confidence.level = 0.95) {
  idx <- stan.ext %>% 
    apply(2, quantile, c((1 - confidence.level) / 2, (1 + confidence.level) / 2)) %>%
    apply(2, prod)
  return((idx <= 0))
}

param.stan <- function(stan.obj, param = "beta") {
  fit.beta <- summary(stan.obj)$summary[grep(param, rownames(summary(stan.obj)$summary)), "mean"]
  la <- rstan::extract(stan.obj, permuted = T)
  beta <- fit.beta * as.numeric(!shrink.stan(la[[param]]))
  return(beta)
}
